# Tutorial 4 — Parallel orbital export with multiple CPU workers

This tutorial explains how to speed up batch orbital rendering using the `--jobs` flag,
which runs multiple orbitals simultaneously across CPU cores.

## When to use this

- Exporting a large active space (e.g. CAS(12,10) with 10 neighbors = ~30 orbitals)
- Generating images for many geometries along a dissociation curve
- Any situation where sequential rendering feels slow

## How it works

Each orbital is completely independent: the script builds the 3D grid once, then sends
one orbital per worker process. Each worker:

1. Opens its own copy of the HDF5 file
2. Computes the orbital on the pre-built grid (pure CPU/numpy work)
3. Runs VTK offscreen rendering to produce the isosurface image
4. Writes the PNG and exits

Workers are OS processes — one per CPU core. No MPI, no GPU required.

## Basic usage

```bash
# Sequential (default — safe, predictable output order)
python pegamoid_batch_export.py calculation.rasscf.h5 --active --neighbors 10

# 4 parallel workers (one per core)
python pegamoid_batch_export.py calculation.rasscf.h5 --active --neighbors 10 --jobs 4
```

Output files and ordering are identical between `--jobs 1` and `--jobs N` — the montage
comes out in the same orbital order regardless of which worker finished first.

## Benchmark (N2, CAS(10,8), ±5 neighbors = 15 orbitals, resolution 50³)

| `--jobs` | Wall time | Speedup | Efficiency |
|----------|-----------|---------|------------|
| 1        | 13.7 s    | 1.00×   | 100%       |
| 4        | 4.6 s     | 3.01×   | 75%        |

The ~25% efficiency loss is one-time process spawn + HDF5 open overhead per worker.
For larger active spaces with many orbitals, efficiency approaches 100% of linear speedup.

## How many workers to use

- **Local workstation**: `--jobs` = number of physical cores (check with `nproc`)
- **HPC (PBS/Torque)**: request `ppn=N` and pass `--jobs N`
- **Rule of thumb**: `--jobs 4` is a safe default that works well on most machines

```bash
# Check your local core count
nproc

# Use all available cores
python pegamoid_batch_export.py calculation.rasscf.h5 --active --neighbors 10 \
    --jobs $(nproc)
```

## HPC example (PBS job script)

```bash
#!/bin/bash
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=8gb
#PBS -A starting_2025_097

source /path/to/setup_hortense.sh

python /path/to/pegamoid_batch_export.py \
    autocas_project/final/system_0.rasscf.h5 \
    --active \
    --neighbors 10 \
    --camera 45,30 \
    --resolution 60 \
    --jobs 4 \
    --output-dir ./orbitals/system_0
```

## Benchmark script

A ready-made benchmark is included in the repository:

```bash
python benchmark_parallel.py calculation.rasscf.h5 \
    --resolution 50 \
    --neighbors 5 \
    --jobs-list 1,2,4
```

This runs the export at each job count, reports wall time, speedup, and parallel efficiency.

## Notes

- Print output from workers will be interleaved — this is normal. Each line is prefixed
  with the worker's process ID (`[pid NNNNN]`) so you can trace which orbital is running.
- The grid arrays (~5 MB at resolution=60) are sent to each worker via pickle once per
  orbital. At resolution=80 this grows to ~12 MB — still negligible overhead.
- Multiple workers reading the same HDF5 file simultaneously is safe; h5py opens files
  read-only and each worker has its own file handle.
