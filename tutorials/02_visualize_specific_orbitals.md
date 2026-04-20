# Tutorial 2 — Visualize a chosen set of orbitals by index

This tutorial shows how to render a specific set of orbitals identified by their
1-based global index, independent of their type (inactive, active, or virtual).

## When to use this

- Inspecting a range of orbitals to find a specific MO character (e.g. searching for
  π bonding orbitals buried in the inactive or virtual space)
- Checking candidate orbitals before using the OpenMolcas `ALTER` keyword to swap them
  into the active space
- Reproducing a specific set of orbitals from a literature reference

## Basic usage

Pass a comma-separated list of 1-based orbital indices:

```bash
python pegamoid_batch_export.py calculation.rasscf.h5 \
    --orbitals 34,35,41,42,43,44,45,46,47,48 \
    --camera 45,30 \
    --output-dir ./orbital_images
```

## Exploring a range

To scan a range of orbitals (e.g. inactive orbitals 19–31) use shell expansion or a
Python one-liner to build the list:

```bash
# Shell brace expansion
python pegamoid_batch_export.py calculation.rasscf.h5 \
    --orbitals $(seq -s, 19 31) \
    --camera 45,30

# Or write it out explicitly
python pegamoid_batch_export.py calculation.rasscf.h5 \
    --orbitals 19,20,21,22,23,24,25,26,27,28,29,30,31 \
    --camera 45,30 \
    --output-dir ./scan_19_to_31
```

## Sorting

By default orbitals are sorted by energy. To preserve the order you specified:

```bash
python pegamoid_batch_export.py calculation.rasscf.h5 \
    --orbitals 46,47,34,35 \
    --sort index \
    --camera 45,30
```

## Combining with custom layout

If you are rendering many orbitals at once, adjust the number of columns in the montage:

```bash
python pegamoid_batch_export.py calculation.rasscf.h5 \
    --orbitals $(seq -s, 30 50) \
    --montage-cols 7 \
    --camera 45,30 \
    --output-dir ./wide_scan
```

## Full example: searching for π bonding MOs in PbO

The active Mertens CAS(8,7) calculation contains π* antibonding orbitals at positions
46 and 47, but the π bonding pair may be at lower energy in the inactive space. To search:

```bash
# Step 1 — scan the inactive region
python pegamoid_batch_export.py PbO_cas87_mertens.rasscf.h5 \
    --orbitals $(seq -s, 19 41) \
    --camera 45,30 \
    --output-dir ./pbo_inactive_scan

# Step 2 — once located, render the candidates together with the active space
python pegamoid_batch_export.py PbO_cas87_mertens.rasscf.h5 \
    --orbitals 34,35,42,43,44,45,46,47,48 \
    --camera 45,30 \
    --output-dir ./pbo_candidates
```

Inspect the montage: orbitals with π character (two lobes, one nodal plane along the bond)
at low energy are good swap-in candidates for the `ALTER` keyword.
