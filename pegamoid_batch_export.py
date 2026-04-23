#!/usr/bin/env python3
"""
Batch orbital export using Pegamoid's Orbitals class for grid computation
and VTK offscreen rendering for isosurface visualization.

Usage:
  python pegamoid_batch_export.py <file.h5> --active --neighbors 2
  python pegamoid_batch_export.py <file.h5> --orbitals 41,42,43,44,45,46
  python pegamoid_batch_export.py <file.molden> --orbitals 5,6,7,8

Requires: Pegamoid==2.12.4, numpy, vtk, h5py, Pillow (for montage)
License: GPL v3.0 (derives from Pegamoid by Ignacio Fdez. Galvan)
"""

import argparse
import os
import shutil
import sys
import types

import numpy as np

# ---------------------------------------------------------------------------
# Pegamoid import
# ---------------------------------------------------------------------------
# Pegamoid installs as a single executable script (pegamoid.py), not as an
# importable Python module.  Its module-level code unconditionally creates a
# QApplication and starts the Qt event loop, so a normal import would launch
# the GUI and block.  We work around this by loading the source, truncating
# it before the QApplication line, and exec-ing the result.
#
# The Pegamoid version is pinned to 2.12.x because the Orbitals class API
# and the module-level code layout may change in future releases.
# ---------------------------------------------------------------------------

_PEGAMOID_SUPPORTED_VERSIONS = ('2.12',)


def _find_pegamoid_script():
    """Locate pegamoid.py installed by pip."""
    import importlib.metadata

    # 1. Check that Pegamoid is installed and the version is supported
    try:
        version = importlib.metadata.version('Pegamoid')
    except importlib.metadata.PackageNotFoundError:
        sys.exit(
            "Error: Pegamoid is not installed.\n"
            "Install it with:  pip install Pegamoid==2.12.4"
        )

    if not any(version.startswith(v) for v in _PEGAMOID_SUPPORTED_VERSIONS):
        print(
            f"Warning: Pegamoid {version} found, but this script was tested "
            f"with 2.12.x.  Proceed with caution.",
            file=sys.stderr,
        )

    # 2. Find pegamoid.py on PATH
    path = shutil.which('pegamoid.py')
    if path and os.path.isfile(path):
        return path

    # 3. Fall back to the pip scripts directory
    import sysconfig
    scripts_dir = sysconfig.get_path('scripts')
    if scripts_dir:
        candidate = os.path.join(scripts_dir, 'pegamoid.py')
        if os.path.isfile(candidate):
            return candidate

    sys.exit(
        "Error: pegamoid.py script not found.\n"
        "Make sure Pegamoid is installed:  pip install Pegamoid==2.12.4"
    )


def _import_pegamoid_orbitals():
    """Import only the Orbitals class from pegamoid.py, skipping the GUI."""
    pegamoid_path = _find_pegamoid_script()

    with open(pegamoid_path, 'r') as f:
        source = f.read()

    # Truncate before the module-level QApplication line
    marker = '\napp = QApplication(sys.argv)'
    idx = source.find(marker)
    if idx > 0:
        source = source[:idx]

    code = compile(source, pegamoid_path, 'exec')
    mod = types.ModuleType('pegamoid')
    mod.__file__ = pegamoid_path
    exec(code, mod.__dict__)
    return mod


_pegamoid = _import_pegamoid_orbitals()
Orbitals = _pegamoid.Orbitals

# On headless systems (no DISPLAY/WAYLAND_DISPLAY) force VTK to use EGL
# offscreen rendering before the window factory is initialised.
# Must come before "import vtk" — without this, VTK defaults to
# vtkXOpenGLRenderWindow which requires a DISPLAY and crashes on HPC.
_headless = not (os.environ.get('DISPLAY') or os.environ.get('WAYLAND_DISPLAY'))
if _headless:
    os.environ.setdefault('VTK_DEFAULT_RENDER_WINDOW_OFFSCREEN', 'true')

import vtk
from vtk.util import numpy_support


# ---------------------------------------------------------------------------
# Orbital info & selection
# ---------------------------------------------------------------------------

def get_orbital_info(filepath):
    """Read orbital type indices, energies, and occupations from HDF5."""
    import h5py
    info = []
    with h5py.File(filepath, 'r') as f:
        if 'MO_TYPEINDICES' in f:
            types = [x.decode() for x in f['MO_TYPEINDICES'][:]]
        else:
            types = ['?' for _ in f['MO_ENERGIES'][:]]
        energies = f['MO_ENERGIES'][:] if 'MO_ENERGIES' in f else []
        occupations = f['MO_OCCUPATIONS'][:] if 'MO_OCCUPATIONS' in f else []
        for i, (t, e, o) in enumerate(zip(types, energies, occupations)):
            info.append({'index': i, 'type': t, 'energy': e, 'occup': o})
    return info


def select_orbitals(orbital_info, active=False, neighbors=0, manual=None):
    """Determine which orbital indices to export."""
    if manual is not None:
        return sorted(manual)

    if not active:
        # Default: all occupied + first few virtuals
        occ = [o['index'] for o in orbital_info if o['occup'] > 0.01]
        if occ:
            return occ
        return list(range(min(10, len(orbital_info))))

    # Find active space orbitals (type '2' = RAS2, '1' = RAS1, '3' = RAS3)
    active_types = {'1', '2', '3'}
    active_indices = [o['index'] for o in orbital_info if o['type'] in active_types]

    if not active_indices:
        print("Warning: No active orbitals found in type indices. Falling back to occupied orbitals.")
        return [o['index'] for o in orbital_info if o['occup'] > 0.01]

    if neighbors > 0:
        # Group active orbitals into contiguous blocks (one per irrep)
        blocks = []
        current_block = [active_indices[0]]
        for i in range(1, len(active_indices)):
            if active_indices[i] - active_indices[i - 1] <= 1:
                current_block.append(active_indices[i])
            else:
                blocks.append(current_block)
                current_block = [active_indices[i]]
        blocks.append(current_block)

        # For each block, add N neighbors before and after
        selected = []
        for block in blocks:
            first = min(block)
            last = max(block)
            start = max(0, first - neighbors)
            end = min(len(orbital_info) - 1, last + neighbors)
            selected.extend(range(start, end + 1))
        return sorted(set(selected))

    return active_indices


# ---------------------------------------------------------------------------
# Grid & rendering
# ---------------------------------------------------------------------------

def build_grid(orbitals_obj, resolution, padding):
    """Build a 3D grid of coordinates covering the molecule."""
    coords = np.array([c['xyz'] for c in orbitals_obj.centers])
    # Coordinates are in Bohr in the HDF5 file
    lo = coords.min(axis=0) - padding
    hi = coords.max(axis=0) + padding

    x = np.linspace(lo[0], hi[0], resolution)
    y = np.linspace(lo[1], hi[1], resolution)
    z = np.linspace(lo[2], hi[2], resolution)

    spacing = [(hi[i] - lo[i]) / (resolution - 1) for i in range(3)]
    origin = lo.tolist()

    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    # VTK vtkImageData expects x-fastest (Fortran) ordering
    return (X.flatten(order='F'), Y.flatten(order='F'), Z.flatten(order='F'),
            origin, spacing, (resolution, resolution, resolution))


def create_molecule_actors(orbitals_obj, atom_radius=0.5, atom_color=(0.85, 0.85, 0.85)):
    """Create VTK actors for atom nuclei (spheres in Bohr coordinates)."""
    actors = []
    for c in orbitals_obj.centers:
        if c['Z'] == 0:
            continue
        sphere = vtk.vtkSphereSource()
        sphere.SetCenter(*c['xyz'])  # Already in Bohr, same as grid
        sphere.SetRadius(atom_radius)
        sphere.SetThetaResolution(24)
        sphere.SetPhiResolution(24)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(sphere.GetOutputPort())

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*atom_color)
        actor.GetProperty().SetAmbient(0.3)
        actor.GetProperty().SetDiffuse(0.7)
        actor.GetProperty().SetSpecular(0.2)
        actors.append(actor)

    return actors


def render_orbital(orbitals_obj, orb_index, xf, yf, zf, origin, spacing, dims,
                   isovalue=0.03, width=800, height=600,
                   pos_color=(0.12, 0.35, 0.78), neg_color=(0.78, 0.24, 0.12),
                   bg_color=(1, 1, 1), opacity=0.7, camera='x'):
    """Compute orbital on grid, create isosurface, render offscreen, return image data."""
    # Compute orbital values on grid
    print(f"  Computing orbital {orb_index + 1} on grid...", end=' ', flush=True)
    mo_values = orbitals_obj.mo(orb_index, xf, yf, zf)
    print("done.")

    # Create VTK image data (all coordinates in Bohr)
    image_data = vtk.vtkImageData()
    image_data.SetDimensions(dims)
    image_data.SetOrigin(*origin)
    image_data.SetSpacing(*spacing)

    # mo_values is already 1D in x-fastest order (matching VTK), no reshape needed
    vtk_data = numpy_support.numpy_to_vtk(mo_values, deep=True, array_type=vtk.VTK_DOUBLE)
    vtk_data.SetName('orbital')
    image_data.GetPointData().SetScalars(vtk_data)

    # Positive isosurface
    contour_pos = vtk.vtkContourFilter()
    contour_pos.SetInputData(image_data)
    contour_pos.SetValue(0, isovalue)

    mapper_pos = vtk.vtkPolyDataMapper()
    mapper_pos.SetInputConnection(contour_pos.GetOutputPort())
    mapper_pos.ScalarVisibilityOff()

    actor_pos = vtk.vtkActor()
    actor_pos.SetMapper(mapper_pos)
    actor_pos.GetProperty().SetColor(*pos_color)
    actor_pos.GetProperty().SetOpacity(opacity)
    actor_pos.GetProperty().SetAmbient(0.2)
    actor_pos.GetProperty().SetDiffuse(0.7)
    actor_pos.GetProperty().SetSpecular(0.3)

    # Negative isosurface
    contour_neg = vtk.vtkContourFilter()
    contour_neg.SetInputData(image_data)
    contour_neg.SetValue(0, -isovalue)

    mapper_neg = vtk.vtkPolyDataMapper()
    mapper_neg.SetInputConnection(contour_neg.GetOutputPort())
    mapper_neg.ScalarVisibilityOff()

    actor_neg = vtk.vtkActor()
    actor_neg.SetMapper(mapper_neg)
    actor_neg.GetProperty().SetColor(*neg_color)
    actor_neg.GetProperty().SetOpacity(opacity)
    actor_neg.GetProperty().SetAmbient(0.2)
    actor_neg.GetProperty().SetDiffuse(0.7)
    actor_neg.GetProperty().SetSpecular(0.3)

    # Set up renderer
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(*bg_color)
    renderer.AddActor(actor_pos)
    renderer.AddActor(actor_neg)

    # Add molecule
    for atom_actor in create_molecule_actors(orbitals_obj):
        renderer.AddActor(atom_actor)

    # Add light
    light = vtk.vtkLight()
    light.SetPosition(1, 1, 1)
    light.SetFocalPoint(0, 0, 0)
    renderer.AddLight(light)

    # Offscreen render window.
    # On headless nodes, explicitly use EGL so VTK never tries to connect
    # to a DISPLAY. vtkXOpenGLRenderWindow ignores SetOffScreenRendering
    # and bus-errors when no X server is present.
    if _headless and hasattr(vtk, 'vtkEGLRenderWindow'):
        render_window = vtk.vtkEGLRenderWindow()
    else:
        render_window = vtk.vtkRenderWindow()
    render_window.SetOffScreenRendering(True)
    render_window.SetSize(width, height)
    render_window.AddRenderer(renderer)

    renderer.ResetCamera()
    cam = renderer.GetActiveCamera()
    fp = cam.GetFocalPoint()
    dist = cam.GetDistance()

    # Parse camera direction: 'x', 'y', 'z', or 'azimuth,elevation' (degrees)
    if camera in ('x', 'X'):
        cam.SetPosition(fp[0] + dist, fp[1], fp[2])
        cam.SetViewUp(0, 0, 1)
    elif camera in ('y', 'Y'):
        cam.SetPosition(fp[0], fp[1] + dist, fp[2])
        cam.SetViewUp(0, 0, 1)
    elif camera in ('z', 'Z'):
        cam.SetPosition(fp[0], fp[1], fp[2] + dist)
        cam.SetViewUp(0, 1, 0)
    else:
        # azimuth,elevation in degrees
        parts = [float(x) for x in camera.split(',')]
        az = np.radians(parts[0])
        el = np.radians(parts[1]) if len(parts) > 1 else 0.0
        cam.SetPosition(
            fp[0] + dist * np.cos(el) * np.cos(az),
            fp[1] + dist * np.cos(el) * np.sin(az),
            fp[2] + dist * np.sin(el))
        cam.SetViewUp(0, 0, 1)

    renderer.ResetCameraClippingRange()
    render_window.Render()

    # Capture to image
    w2i = vtk.vtkWindowToImageFilter()
    w2i.SetInput(render_window)
    w2i.SetInputBufferTypeToRGB()
    w2i.Update()

    return w2i.GetOutput()


def save_vtk_image_as_png(vtk_image, filepath):
    """Save a VTK image data object as PNG."""
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(filepath)
    writer.SetInputData(vtk_image)
    writer.SetCompressionLevel(9)
    writer.Write()


# ---------------------------------------------------------------------------
# Montage
# ---------------------------------------------------------------------------

def create_montage(image_paths, labels, output_path, cols=4, img_size=(400, 300),
                   label_height=40, font_size=14):
    """Create a montage grid image from individual PNGs using Pillow."""
    from PIL import Image, ImageDraw, ImageFont

    n = len(image_paths)
    rows = (n + cols - 1) // cols

    cell_w = img_size[0]
    cell_h = img_size[1] + label_height
    montage_w = cols * cell_w
    montage_h = rows * cell_h

    montage = Image.new('RGB', (montage_w, montage_h), 'white')
    draw = ImageDraw.Draw(montage)

    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf", font_size)
    except (IOError, OSError):
        try:
            font = ImageFont.truetype("C:/Windows/Fonts/consola.ttf", font_size)
        except (IOError, OSError):
            font = ImageFont.load_default()

    for i, (img_path, label) in enumerate(zip(image_paths, labels)):
        row = i // cols
        col = i % cols
        x = col * cell_w
        y = row * cell_h

        img = Image.open(img_path)
        img = img.resize(img_size, Image.LANCZOS)
        montage.paste(img, (x, y))

        # Draw label below image
        draw.text((x + 5, y + img_size[1] + 2), label, fill='black', font=font)

    montage.save(output_path)
    print(f"Montage saved to: {output_path}")


# ---------------------------------------------------------------------------
# File format detection
# ---------------------------------------------------------------------------

def detect_format(filepath):
    """Detect file format from extension."""
    ext = os.path.splitext(filepath)[1].lower()
    if ext in ('.h5', '.hdf5'):
        return 'hdf5'
    elif ext in ('.molden',):
        return 'molden'
    elif ext in ('.scforb', '.rassorb', '.pt2orb', '.gssorb', '.lprorb'):
        return 'inporb'
    else:
        # Try HDF5
        import h5py
        if h5py.is_hdf5(filepath):
            return 'hdf5'
        return 'molden'


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Batch export orbital images from OpenMolcas files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  %(prog)s calculation.rasscf.h5 --active --neighbors 2
  %(prog)s calculation.rasscf.h5 --orbitals 41,42,43,44,45,46 --camera 20,15
  %(prog)s calculation.scf.molden --orbitals 1,2,3,4,5 --isovalue 0.05
""")
    parser.add_argument('file', help='Input file (HDF5 or Molden format)')
    parser.add_argument('--output-dir', '-o', default='./orbitals',
                        help='Output directory for PNG files (default: ./orbitals)')
    parser.add_argument('--orbitals', type=str, default=None,
                        help='Comma-separated list of orbital indices (1-based)')
    parser.add_argument('--active', action=argparse.BooleanOptionalAction, default=True,
                        help='Auto-select active space orbitals from MO_TYPEINDICES (default: on)')
    parser.add_argument('--neighbors', '-n', type=int, default=10,
                        help='Include N neighbor orbitals on each side of active space (default: 10)')
    parser.add_argument('--isovalue', type=float, default=0.03,
                        help='Isosurface cutoff value (default: 0.03)')
    parser.add_argument('--resolution', '-r', type=int, default=60,
                        help='Grid points per axis (default: 60)')
    parser.add_argument('--padding', type=float, default=10.0,
                        help='Extra space around molecule bounding box in Bohr (default: 10.0)')
    parser.add_argument('--width', type=int, default=800, help='Image width (default: 800)')
    parser.add_argument('--height', type=int, default=600, help='Image height (default: 600)')
    parser.add_argument('--montage-cols', type=int, default=4,
                        help='Number of columns in montage (default: 4)')
    parser.add_argument('--no-montage', action='store_true', help='Skip montage creation')
    parser.add_argument('--opacity', type=float, default=0.7, help='Orbital opacity (default: 0.7)')
    parser.add_argument('--camera', type=str, default='45,30',
                        help='Camera direction: x, y, z, or azimuth,elevation in degrees '
                             '(default: 45,30). Examples: "x", "z", "20,15"')
    parser.add_argument('--sort', type=str, default='energy', choices=['energy', 'index'],
                        help='Sort orbitals by energy or index (default: energy)')
    parser.add_argument('--colorful-background', action=argparse.BooleanOptionalAction, default=True,
                        help='Color background by orbital type: pastel green for active '
                             '(RAS1/2/3), light blue for core/virtual (default: on)')

    args = parser.parse_args()

    filepath = os.path.abspath(args.file)
    if not os.path.exists(filepath):
        print(f"Error: File not found: {filepath}")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    # Detect format and load orbitals
    fmt = detect_format(filepath)
    print(f"Loading {fmt} file: {filepath}")
    orb = Orbitals(filepath, fmt)

    # Determine orbital selection
    manual_indices = None
    if args.orbitals:
        manual_indices = [int(x) - 1 for x in args.orbitals.split(',')]  # Convert to 0-based

    if fmt == 'hdf5':
        orbital_info = get_orbital_info(filepath)
        selected = select_orbitals(orbital_info, active=args.active,
                                   neighbors=args.neighbors, manual=manual_indices)
    else:
        # For molden files, no type indices available
        if manual_indices is not None:
            selected = manual_indices
        else:
            n_orbs = len(orb.MO)
            selected = list(range(min(n_orbs, 10)))

    # Sort selected orbitals
    if fmt == 'hdf5' and args.sort == 'energy':
        selected = sorted(selected, key=lambda i: orbital_info[i]['energy'])
        print(f"Selected orbitals (1-based, sorted by energy): {[i + 1 for i in selected]}")
    else:
        print(f"Selected orbitals (1-based): {[i + 1 for i in selected]}")

    # Build grid once (reuse for all orbitals)
    print(f"Building grid: {args.resolution}^3 points, padding={args.padding} Bohr")
    xf, yf, zf, origin, spacing, dims = build_grid(orb, args.resolution, args.padding)

    # Render each orbital
    image_paths = []
    labels = []
    basename = os.path.splitext(os.path.basename(filepath))[0]

    for idx in selected:
        orb_num = idx + 1  # 1-based for display

        # Build label
        if fmt == 'hdf5' and idx < len(orbital_info):
            info = orbital_info[idx]
            label = f"MO {orb_num} | {info['type']} | E={info['energy']:.4f} | occ={info['occup']:.2f}"
        else:
            label = f"MO {orb_num}"

        print(f"Rendering {label}")

        if args.colorful_background and fmt == 'hdf5' and idx < len(orbital_info):
            orb_type = orbital_info[idx]['type']
            if orb_type in ('1', '2', '3'):
                bg = (0.78, 0.94, 0.78)   # pastel green — active (RAS1/2/3)
            else:
                bg = (0.80, 0.88, 0.96)   # pastel blue — core / virtual
        else:
            bg = (1.0, 1.0, 1.0)

        vtk_image = render_orbital(orb, idx, xf, yf, zf, origin, spacing, dims,
                                   isovalue=args.isovalue,
                                   width=args.width, height=args.height,
                                   bg_color=bg,
                                   opacity=args.opacity, camera=args.camera)

        png_path = os.path.join(args.output_dir, f"{basename}_mo{orb_num:03d}.png")
        save_vtk_image_as_png(vtk_image, png_path)
        print(f"  Saved: {png_path}")

        image_paths.append(png_path)
        labels.append(label)

    # Create montage
    if not args.no_montage and len(image_paths) > 1:
        montage_path = os.path.join(args.output_dir, f"{basename}_montage.png")
        try:
            create_montage(image_paths, labels, montage_path,
                           cols=args.montage_cols,
                           img_size=(args.width // 2, args.height // 2))
        except ImportError:
            print("Warning: Pillow not installed, skipping montage. Install with: pip install Pillow")

    print(f"\nDone! {len(image_paths)} orbital images saved to {args.output_dir}/")


if __name__ == '__main__':
    main()
