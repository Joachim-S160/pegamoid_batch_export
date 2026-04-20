# Tutorial 1 — Visualize the CAS active space with neighbours

This tutorial shows how to export the active space orbitals from a RASSCF or CASSCF
calculation, together with a window of neighbouring orbitals on each side, to get
context on what surrounds the active space.

## When to use this

- Inspecting CASSCF/RASSCF natural orbital character after a calculation
- Checking whether the right orbitals are in the active space
- Identifying orbital contamination (d/f orbitals in a valence active space)
- Comparing active space content across a dissociation curve

## Prerequisites

```bash
pip install -r requirements.txt
```

## Basic usage

Export the active space (detected automatically from `MO_TYPEINDICES`) plus 10
neighbours on each side:

```bash
python pegamoid_batch_export.py calculation.rasscf.h5 \
    --active \
    --neighbors 10 \
    --camera 45,30 \
    --output-dir ./orbital_images
```

This produces:
- `./orbital_images/calculation.rasscf_mo<NNN>.png` — one PNG per orbital
- `./orbital_images/calculation.rasscf_montage.png` — all orbitals in a grid

## Camera angle

For diatomic molecules along the x-axis, `--camera 45,30` gives a useful 3D perspective.
For molecules in the xy plane, try `--camera 0,90` (top-down). The default `--camera x`
looks straight along x.

```bash
# Diatomic (bond along x-axis)
python pegamoid_batch_export.py PbO.rasscf.h5 --active --neighbors 10 --camera 45,30

# Planar molecule (in xy plane)
python pegamoid_batch_export.py benzene.rasscf.h5 --active --neighbors 5 --camera 0,90

# Default side view
python pegamoid_batch_export.py N2.rasscf.h5 --active --neighbors 10
```

## Reading the montage

Each orbital panel in the montage is labelled:

```
MO 42  type=2  E=-0.421  occ=1.847
```

| Field | Meaning |
|-------|---------|
| `MO 42` | 1-based global orbital index |
| `type=2` | orbital class: I=inactive, 1/2/3=RAS1/2/3 active, S=secondary |
| `E=...` | orbital energy (Hartree) |
| `occ=...` | occupation number (RASSCF natural orbital) |

Orbitals with `type=2` (or `1`/`3`) are in the active space. The neighbours have
`type=I` (occupied inactive) or `type=S` (virtual secondary).

## Adjusting the isovalue

The default isovalue (0.03) works well for most valence orbitals. For diffuse or Rydberg
orbitals, lower it:

```bash
python pegamoid_batch_export.py calculation.rasscf.h5 --active --neighbors 5 --isovalue 0.015
```

## Full example: PbO CAS(8,7) at equilibrium geometry

```bash
python pegamoid_batch_export.py \
    PbO_cas87_mertens.rasscf.h5 \
    --active \
    --neighbors 10 \
    --camera 45,30 \
    --output-dir ./pbo_cas87_active
```

The montage will show 7 active orbitals flanked by 10 inactive and 10 virtual neighbours
(27 orbitals total), allowing you to assess the quality of the active space selection.
