#!/usr/bin/env python3
"""
Benchmark --jobs 1 vs --jobs 4 for pegamoid_batch_export.py.

Usage:
  source /path/to/setup_autocas_env.sh
  python benchmark_parallel.py <file.rasscf.h5> [--resolution 60] [--neighbors 5]
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time


SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pegamoid_batch_export.py')


def run_export(h5file, outdir, jobs, resolution, neighbors):
    cmd = [
        sys.executable, SCRIPT,
        h5file,
        '--output-dir', outdir,
        '--active',
        '--neighbors', str(neighbors),
        '--resolution', str(resolution),
        '--jobs', str(jobs),
        '--no-montage',
        '--camera', '45,30',
    ]
    t0 = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        print(f"ERROR (jobs={jobs}):\n{result.stderr[-2000:]}", file=sys.stderr)
        sys.exit(1)
    return elapsed


def count_pngs(directory):
    return len([f for f in os.listdir(directory) if f.endswith('.png')])


def main():
    parser = argparse.ArgumentParser(description='Benchmark parallel vs sequential orbital export')
    parser.add_argument('file', help='Input HDF5 file (e.g. system_0.rasscf.h5)')
    parser.add_argument('--resolution', '-r', type=int, default=50,
                        help='Grid resolution (default: 50 — lower = faster benchmark)')
    parser.add_argument('--neighbors', '-n', type=int, default=5,
                        help='Neighbor orbitals on each side of active space (default: 5)')
    parser.add_argument('--jobs-list', type=str, default='1,4',
                        help='Comma-separated list of job counts to test (default: 1,4)')
    args = parser.parse_args()

    jobs_list = [int(j) for j in args.jobs_list.split(',')]
    h5file = os.path.abspath(args.file)

    print(f"File       : {h5file}")
    print(f"Resolution : {args.resolution}^3")
    print(f"Neighbors  : ±{args.neighbors}")
    print(f"Jobs tested: {jobs_list}")
    print()

    timings = {}
    with tempfile.TemporaryDirectory(prefix='pegamoid_bench_') as tmproot:
        for jobs in jobs_list:
            outdir = os.path.join(tmproot, f'jobs_{jobs}')
            os.makedirs(outdir)

            print(f"--- Running with --jobs {jobs} ---")
            elapsed = run_export(h5file, outdir, jobs, args.resolution, args.neighbors)
            n_imgs = count_pngs(outdir)
            timings[jobs] = elapsed
            print(f"  Orbitals rendered : {n_imgs}")
            print(f"  Wall time         : {elapsed:.1f} s")
            print()

    print("=" * 50)
    print("BENCHMARK SUMMARY")
    print("=" * 50)
    baseline = timings[jobs_list[0]]
    for jobs in jobs_list:
        t = timings[jobs]
        speedup = baseline / t
        efficiency = speedup / jobs * 100
        print(f"  jobs={jobs:2d}  {t:6.1f} s   speedup={speedup:.2f}x   efficiency={efficiency:.0f}%")
    print()
    if len(jobs_list) > 1:
        best_jobs = min(timings, key=timings.get)
        print(f"Fastest: --jobs {best_jobs} ({timings[best_jobs]:.1f} s)")


if __name__ == '__main__':
    main()
