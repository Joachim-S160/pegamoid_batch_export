# Tutorial 5 — Interactive dissociation dashboard

This tutorial shows how to build a localhost orbital inspection dashboard for a
CASSCF dissociation curve. Starting from a set of `rasscf.h5` files, you will render
orbital images for every geometry, generate the data file the dashboard reads, and
serve an interactive web page where you can click any geometry to inspect its active
space orbitals.

## What the dashboard shows

- Dissociation curve: relative energy in eV vs. bond distance (Å)
- Per-point orbital panel: active space orbitals + ±5 neighbours, with natural
  orbital occupation numbers
- Geometry navigation with ← → buttons (and keyboard arrow keys)
- Optional two-stage mode: toggle between two sets of h5 files (e.g. a per-geometry
  active space vs. a union active space) to see how orbital character changes

## When to use this

- Checking whether your active space follows the correct orbitals as the bond stretches
- Spotting sudden orbital rotations or active space discontinuities along a curve
- Comparing two pipeline stages (autoCAS per-geometry vs. union CAS)
- Generating a visual supplement for a dissociation curve publication or presentation

---

## Prerequisites

```bash
pip install -r requirements.txt     # numpy, vtk, h5py, Pillow
```

You need a set of `*.rasscf.h5` files from a CASSCF dissociation curve calculation,
one file per geometry.  They can come from OpenMolcas, or any program that writes
the same HDF5 schema (`MO_TYPEINDICES`, `MO_OCCUPATIONS`, `ROOT_ENERGIES`,
`CENTER_COORDINATES`).

---

## Step 1 — Render orbital images

Use `pegamoid_batch_export.py` on each h5 file. The `run_pegamoid_parallel.py`
helper script (in `autoCAS4HE_priv/scripts/`) runs multiple systems concurrently.

For a manual loop (single CPU):

```bash
for i in $(seq 0 24); do
    python pegamoid_batch_export.py \
        results/system_${i}.rasscf.h5 \
        --active \
        --neighbors 5 \
        --camera 45,30 \
        --jobs 4 \
        --output-dir dashboard/static/images/main/system_${i}
done
```

For parallel rendering across systems:

```bash
python autoCAS4HE_priv/scripts/run_pegamoid_parallel.py \
    --results-dir results \
    --workers 4
```

This produces one PNG per orbital per geometry in
`dashboard/static/images/main/system_N/`.

### White-background fix

VTK renders on a black background. The parallel runner applies a PIL fix
automatically. If you use a manual loop, apply it after each render:

```python
from PIL import Image
import numpy as np, pathlib

for p in pathlib.Path("dashboard/static/images/main/system_0").glob("*.png"):
    img = Image.open(p).convert("RGB")
    arr = np.array(img)
    arr[np.all(arr < 15, axis=-1)] = 255
    Image.fromarray(arr).save(p)
```

---

## Step 2 — Generate `data.json`

### Single-stage mode

Place all h5 files in one directory (or use the directory directly):

```
results/
    system_0.rasscf.h5
    system_1.rasscf.h5
    ...
```

```bash
python generate_dashboard_data.py \
    --h5-dir results/ \
    --molecule N2 \
    --basis ANO-RCC-VDZP \
    --method CASSCF \
    --out dashboard/data.json
```

Bond distances are extracted automatically from `CENTER_COORDINATES` for diatomics.
For polyatomics, supply a CSV:

```bash
python generate_dashboard_data.py \
    --h5-dir results/ \
    --geometry-csv geometries.csv \
    --molecule H2O ...
```

`geometries.csv` format:
```
id,distance_ang
system_0,2.00
system_1,1.80
system_2,1.60
```

### Two-stage mode

If your results directory has `pergeom/` and `final/` subdirectories, the dashboard
automatically shows a stage toggle:

```
results/
    pergeom/
        system_0.rasscf.h5  ...
    final/
        system_0.rasscf.h5  ...
```

```bash
python generate_dashboard_data.py \
    --h5-dir results/ \
    --molecule PbO \
    --basis ANO-RCC-VQZP \
    --method CASSCF \
    --exp-de-ev 3.87 \
    --exp-ref "Drowart 1965" \
    --stage-label-pergeom "Per-geometry CAS" \
    --stage-label-final "Union CAS" \
    --images-dir dashboard/static/images \
    --out dashboard/data.json
```

For two-stage mode, render images into `dashboard/static/images/pergeom/system_N/`
and `dashboard/static/images/final/system_N/` (matching the stage keys).

### Per-geometry active space override (autoCAS pipelines)

When the `pergeom/` h5 files come from an autoCAS pipeline, the CASSCF is typically
run with the **union** orbital window (e.g. CAS(10,8) for all geometries), even though
autoCAS selected a smaller per-geometry active space at some geometries (e.g. CAS(6,6)
at the dissociation limit). The h5 `MO_TYPEINDICES` reflects the union window, so the
script would wrongly report CAS(10,8) everywhere.

To get the correct per-geometry sizes, provide the `per_geometry_cas_spaces` file
written by the autoCAS pipeline:

```bash
python generate_dashboard_data.py \
    --h5-dir results/ \
    --molecule PbO ...
    # per_geometry_cas_spaces is auto-detected if it exists in --h5-dir
```

The file is auto-detected when it exists at `<h5-dir>/per_geometry_cas_spaces`.
You can also point to it explicitly:

```bash
python generate_dashboard_data.py \
    --h5-dir results/ \
    --pergeom-cas path/to/per_geometry_cas_spaces \
    ...
```

`per_geometry_cas_spaces` format (written by autoCAS):
```
system system_0: CAS(6,6) indices: [41, 42, 44, 45, 46, 47]
system system_0: occupation: [2, 2, 2, 0, 0, 0]
system system_1: CAS(6,6) indices: [41, 42, 44, 45, 46, 47]
...
```

Indices are 0-based Python indices; the script converts them to 1-based MO numbers.
The `CAS(n_e, n_o)` values are used directly for the dashboard labels.
If no such file is found, the script falls back to reading active MOs from h5
`MO_TYPEINDICES` — appropriate when each h5 was computed with its own active space.

---

## Step 3 — Serve the dashboard

```bash
python dashboard/app.py
```

Open http://localhost:8080 in your browser.

```
Click any point → orbital panel opens on the right
← / →  (or keyboard arrow keys) → navigate between geometries
Click an orbital image → full-size lightbox
Click "Neighbouring orbitals" → expand ±5 neighbours
```

---

## Worked example — PbO dissociation, two-stage

This example uses the PbO data from the autoCAS4HE benchmarks
(`tests/autocas/external_scf/Pb_benchmark/PbO/results_N25_pergeom/`).

```bash
# Render images for both stages (25 geometries × 2 stages = 50 renders)
python autoCAS4HE_priv/scripts/run_pegamoid_parallel.py \
    --results-dir tests/autocas/external_scf/Pb_benchmark/PbO/results_N25_pergeom \
    --workers 4

# Generate data.json
python generate_dashboard_data.py \
    --h5-dir tests/autocas/external_scf/Pb_benchmark/PbO/results_N25_pergeom \
    --molecule PbO --basis ANO-RCC-VQZP --method CASSCF \
    --exp-de-ev 3.87 --exp-ref "Drowart 1965" \
    --stage-label-pergeom "Per-geometry CAS" \
    --stage-label-final "Union CAS" \
    --out dashboard/data.json

# Serve
python dashboard/app.py
```

Expected output:
- 25 geometries, R = 1.55–8.00 Å
- Equilibrium near R ≈ 1.90 Å (ΔE ≈ 0, minimum of curve)
- Dissociation limit at top of curve (ΔE ≈ 3.5–4.0 eV)
- Per-geometry stage: active space grows from CAS(6,8) at large R → CAS(10,8) at short R
- Union stage: uniform CAS(10,8) across all geometries

---

## data.json schema reference

The dashboard reads `data.json` from the same directory as `index.html`. The schema:

```json
{
  "molecule": "PbO",
  "basis": "ANO-RCC-VQZP",
  "method": "CASSCF",
  "stage_labels": {
    "pergeom": "Per-geometry CAS",
    "final": "Union CAS"
  },
  "exp_de_ev": 3.87,
  "exp_ref": "Drowart 1965",
  "points": [
    {
      "id": "system_0",
      "distance": 8.0,
      "pergeom": {
        "energy": -20933.981,
        "cas_electrons": 6,
        "cas_orbitals": 8,
        "active_mo": [41, 42, 43, 44, 45, 46, 47, 48],
        "neighbor_mo_below": [36, 37, 38, 39, 40],
        "neighbor_mo_above": [49, 50, 51, 52, 53],
        "mo_occ": {"41": 1.9948, "42": 1.0, ...},
        "mo_en":  {"41": -2.185, "42": -1.037, ...}
      },
      "final": { "..." }
    }
  ]
}
```

`stage_labels` keys must match the stage keys used in `points[*]` and in the image
path `dashboard/static/images/<stage_key>/system_N/`.

For single-stage mode, `stage_labels` has one entry (key `"main"`) and the stage
toggle is hidden automatically.

---

## Notes

- Energies in `data.json` are stored as absolute values in Eh. The dashboard converts
  them to relative eV (relative to each stage's minimum) for display.
- `generate_dashboard_data.py` reads the active space from `MO_TYPEINDICES` in the h5
  file. This reflects what was actually computed in the CASSCF, not a post-hoc autoCAS
  selection. For autoCAS-specific per-geometry selections, use the PbO-specific
  `generate_website_data.py` instead.
- PNG filenames must match the pattern `<system_id>.rasscf_mo<NNN>.png` (3-digit
  zero-padded 1-based MO index). This is the default output format of
  `pegamoid_batch_export.py`.
