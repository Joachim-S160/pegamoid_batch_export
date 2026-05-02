#!/usr/bin/env python3
"""
Generate dashboard/data.json for the pegamoid dissociation dashboard.

Works directly from OpenMolcas rasscf.h5 files — no autoCAS metadata required.

Two-stage mode (pergeom / final):
  Pass a root directory that contains `pergeom/` and `final/` subdirectories,
  each holding one rasscf.h5 per geometry:
    results/
      pergeom/system_0.rasscf.h5  system_1.rasscf.h5  ...
      final/  system_0.rasscf.h5  system_1.rasscf.h5  ...

Single-stage mode:
  Pass a directory containing rasscf.h5 files directly:
    results/system_0.rasscf.h5  system_1.rasscf.h5  ...

Bond distances are extracted automatically from CENTER_COORDINATES in the h5
(works for any diatomic molecule). For polyatomics, supply --geometry-csv.

Usage:
  python generate_dashboard_data.py --h5-dir results/ \\
      --molecule PbO --basis ANO-RCC-VQZP --method CASSCF \\
      --exp-de-ev 3.87 --exp-ref "Drowart 1965" \\
      --images-dir dashboard/static/images \\
      --out dashboard/data.json
"""

import argparse
import ast
import csv
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np

BOHR_TO_ANG = 0.529177


def h5_sort_key(p: Path) -> int:
    """Sort h5 files by the integer index in their stem (system_N.rasscf.h5 → N)."""
    m = re.search(r"(\d+)", p.stem)
    return int(m.group(1)) if m else 0


def read_distance_from_h5(h5_path: Path) -> float | None:
    """Auto-detect bond distance (Å) for a diatomic from CENTER_COORDINATES."""
    try:
        with h5py.File(h5_path, "r") as f:
            if "CENTER_COORDINATES" not in f:
                return None
            coords = np.array(f["CENTER_COORDINATES"])
            if len(coords) != 2:
                return None
            return float(np.linalg.norm(coords[1] - coords[0])) * BOHR_TO_ANG
    except Exception:
        return None


def infer_molecule_name(h5_path: Path) -> str:
    """Infer molecule name from CENTER_LABELS in h5."""
    try:
        with h5py.File(h5_path, "r") as f:
            if "CENTER_LABELS" not in f:
                return ""
            labels = [lbl.decode().strip().rstrip("0123456789") for lbl in f["CENTER_LABELS"]]
            return "".join(labels)
    except Exception:
        return ""


def parse_per_geometry_cas(path: Path) -> dict[str, tuple[int, int, list[int]]]:
    """Parse a per_geometry_cas_spaces file.

    Returns {system_id: (n_electrons, n_orbitals, [1-based MO numbers])}.
    Indices in the file are 0-based; +1 converts to 1-based MO numbers.
    """
    pattern = re.compile(r"system (system_\d+): CAS\((\d+),(\d+)\) indices: (\[.*?\])")
    result: dict[str, tuple[int, int, list[int]]] = {}
    for line in path.read_text().splitlines():
        m = pattern.match(line)
        if m:
            sid = m.group(1)
            n_elec = int(m.group(2))
            n_orbs = int(m.group(3))
            indices_0based = ast.literal_eval(m.group(4))
            result[sid] = (n_elec, n_orbs, sorted(i + 1 for i in indices_0based))
    return result


def read_h5_point(h5_path: Path, n_neighbors: int,
                  override: tuple[int, int, list[int]] | None = None) -> dict:
    """Read energy, active MOs, occupations, and orbital energies from a rasscf.h5 file."""
    with h5py.File(h5_path, "r") as f:
        # Total energy
        for key in ("ROOT_ENERGIES", "RASSCF_ENERGY", "SCF_ENERGY", "CASPT2_ENERGY"):
            if key in f:
                energy = float(np.array(f[key]).flat[0])
                break
        else:
            raise KeyError(f"No energy key found in {h5_path}")

        # MO types: b'I'=inactive, b'2'=RAS2/active, b'S'=secondary
        mo_types = np.array(f["MO_TYPEINDICES"]).tolist()
        n_mos = len(mo_types)

        occ_arr = np.array(f["MO_OCCUPATIONS"]) if "MO_OCCUPATIONS" in f else None
        en_arr  = np.array(f["MO_ENERGIES"])    if "MO_ENERGIES"    in f else None

    if override is not None:
        n_elec_override, n_orbs_override, active_mo = override
    else:
        n_elec_override = None
        inactive_types = {b'I', b'S', 'I', 'S'}
        active_mo = sorted([i + 1 for i, t in enumerate(mo_types) if t not in inactive_types])

    if not active_mo:
        raise ValueError(f"No active MOs found in {h5_path}")

    first_active = active_mo[0]
    last_active  = active_mo[-1]

    below_start = max(1, first_active - n_neighbors)
    neighbor_mo_below = list(range(below_start, first_active))

    above_end = min(n_mos, last_active + n_neighbors)
    neighbor_mo_above = list(range(last_active + 1, above_end + 1))

    mo_occ, mo_en = {}, {}
    for mo in active_mo + neighbor_mo_below + neighbor_mo_above:
        if occ_arr is not None and mo - 1 < len(occ_arr):
            mo_occ[str(mo)] = round(float(occ_arr[mo - 1]), 4)
        if en_arr is not None and mo - 1 < len(en_arr):
            mo_en[str(mo)] = round(float(en_arr[mo - 1]), 4)

    # Use override electron count when available; otherwise derive from occupation sum.
    # (Override is needed when the h5 was run with a larger union active space: the
    # occupation sum over the per-geometry subset would miss electrons in excluded MOs.)
    if n_elec_override is not None:
        n_elec = n_elec_override
    else:
        n_elec = round(sum(mo_occ.get(str(mo), 0.0) for mo in active_mo))

    return {
        "energy": energy,
        "cas_electrons": n_elec,
        "cas_orbitals": len(active_mo),
        "active_mo": active_mo,
        "neighbor_mo_below": neighbor_mo_below,
        "neighbor_mo_above": neighbor_mo_above,
        "mo_occ": mo_occ,
        "mo_en": mo_en,
    }


def collect_stage(stage_dir: Path, geometry_map: dict[str, float] | None,
                  n_neighbors: int,
                  active_override: dict[str, tuple[int, int, list[int]]] | None = None,
                  ) -> tuple[list[dict], list[str]]:
    """Return (points_list, errors) for all h5 files in stage_dir.

    active_override maps system_id → (n_elec, n_orbs, 1-based MO list); when present,
    overrides the h5-derived active space and electron count (needed when h5 files were
    computed with a larger union active space but the true per-geometry selection is
    stored externally in a per_geometry_cas_spaces file).
    """
    h5_files = sorted(stage_dir.glob("*.rasscf.h5"), key=h5_sort_key)
    if not h5_files:
        return [], [f"No rasscf.h5 files found in {stage_dir}"]

    points = []
    errors = []
    for h5 in h5_files:
        stem = h5.stem.replace(".rasscf", "")   # e.g. "system_0"
        try:
            if geometry_map and stem in geometry_map:
                distance = geometry_map[stem]
            else:
                distance = read_distance_from_h5(h5)
                if distance is None:
                    errors.append(f"{stem}: could not auto-detect distance; skipping")
                    continue
            ov = active_override.get(stem) if active_override else None
            info = read_h5_point(h5, n_neighbors, override=ov)
            points.append({"id": stem, "distance": distance, "info": info})
        except Exception as e:
            errors.append(f"{stem}: {e}")

    return points, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--h5-dir", required=True,
                        help="Root directory: contains h5 files directly (single-stage) "
                             "or pergeom/ and final/ subdirs (two-stage)")
    parser.add_argument("--geometry-csv",
                        help="CSV with columns 'id,distance_ang' (e.g. system_0,1.92). "
                             "Overrides auto-detection for named systems.")
    parser.add_argument("--molecule", default="",
                        help="Molecule name for dashboard header (default: inferred from h5)")
    parser.add_argument("--basis",  default="", help="Basis set label")
    parser.add_argument("--method", default="CASSCF", help="Method label (default: CASSCF)")
    parser.add_argument("--exp-de-ev", type=float, default=None,
                        help="Experimental Dₑ in eV (shown as reference on chart)")
    parser.add_argument("--exp-ref", default="",
                        help="Citation for experimental Dₑ (e.g. 'Drowart 1965')")
    parser.add_argument("--neighbors", type=int, default=5,
                        help="Number of neighbour MOs above/below active space (default: 5)")
    parser.add_argument("--images-dir", default="dashboard/static/images",
                        help="Path to rendered PNG image tree (default: dashboard/static/images)")
    parser.add_argument("--out", default="dashboard/data.json",
                        help="Output path for data.json (default: dashboard/data.json)")
    parser.add_argument("--stage-label-pergeom", default="Per-geometry CAS",
                        help="Label for the pergeom stage in two-stage mode")
    parser.add_argument("--stage-label-final", default="Union CAS",
                        help="Label for the final stage in two-stage mode")
    parser.add_argument("--pergeom-cas",
                        help="Path to per_geometry_cas_spaces file (or equivalent). "
                             "Overrides h5 MO_TYPEINDICES for the pergeom stage. "
                             "Auto-detected if a file named 'per_geometry_cas_spaces' "
                             "exists directly in --h5-dir.")
    args = parser.parse_args()

    h5_dir = Path(args.h5_dir)
    out_path = Path(args.out)

    # Load optional geometry CSV
    geometry_map: dict[str, float] = {}
    if args.geometry_csv:
        with open(args.geometry_csv) as f:
            for row in csv.DictReader(f):
                geometry_map[row["id"]] = float(row["distance_ang"])

    # Detect mode: two-stage if pergeom/ and final/ subdirs exist
    pergeom_dir = h5_dir / "pergeom"
    final_dir   = h5_dir / "final"
    two_stage = pergeom_dir.is_dir() and final_dir.is_dir()

    # Load optional per-geometry active space override (for pergeom stage)
    pergeom_active_override: dict[str, list[int]] | None = None
    pergeom_cas_path = Path(args.pergeom_cas) if args.pergeom_cas else h5_dir / "per_geometry_cas_spaces"
    if pergeom_cas_path.is_file():
        pergeom_active_override = parse_per_geometry_cas(pergeom_cas_path)
        print(f"Per-geometry CAS override loaded: {pergeom_cas_path} ({len(pergeom_active_override)} systems)")

    if two_stage:
        stage_labels = {
            "pergeom": args.stage_label_pergeom,
            "final":   args.stage_label_final,
        }
        print(f"Two-stage mode: {pergeom_dir}  /  {final_dir}")
        pergeom_pts, pergeom_errs = collect_stage(pergeom_dir, geometry_map or None, args.neighbors,
                                                   active_override=pergeom_active_override)
        final_pts,   final_errs   = collect_stage(final_dir,   geometry_map or None, args.neighbors)
        all_errors = pergeom_errs + final_errs

        # Merge by system id
        pergeom_by_id = {p["id"]: p for p in pergeom_pts}
        final_by_id   = {p["id"]: p for p in final_pts}
        all_ids = sorted(set(pergeom_by_id) | set(final_by_id),
                         key=lambda x: h5_sort_key(Path(x)))

        points = []
        for sid in all_ids:
            pg = pergeom_by_id.get(sid)
            fn = final_by_id.get(sid)
            distance = (pg or fn)["distance"]
            pt = {"id": sid, "distance": round(distance, 6)}
            pt["pergeom"] = pg["info"] if pg else None
            pt["final"]   = fn["info"] if fn else None
            points.append(pt)
            pg_cas = f"CAS({pg['info']['cas_electrons']},{pg['info']['cas_orbitals']})" if pg else "—"
            fn_ok  = "ok" if fn else "—"
            print(f"  {sid}  R={distance:.3f} Å  pergeom={pg_cas}  final={fn_ok}")
    else:
        stage_labels = {"main": args.method}
        print(f"Single-stage mode: {h5_dir}")
        pts, all_errors = collect_stage(h5_dir, geometry_map or None, args.neighbors)
        points = []
        for p in pts:
            pt = {"id": p["id"], "distance": round(p["distance"], 6), "main": p["info"]}
            cas = p["info"]
            print(f"  {p['id']}  R={p['distance']:.3f} Å  CAS({cas['cas_electrons']},{cas['cas_orbitals']})")
            points.append(pt)

    # Infer molecule name from first available h5 if not supplied
    molecule = args.molecule
    if not molecule:
        first_h5 = next((h5_dir / "pergeom" if two_stage else h5_dir).glob("*.rasscf.h5"), None)
        if first_h5:
            molecule = infer_molecule_name(first_h5)

    data = {
        "molecule": molecule,
        "basis":  args.basis,
        "method": args.method,
        "stage_labels": stage_labels,
        "points": points,
    }
    if args.exp_de_ev is not None:
        data["exp_de_ev"] = args.exp_de_ev
        data["exp_ref"]   = args.exp_ref

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2))
    print(f"\nWrote {out_path}  ({len(points)} points, {len(stage_labels)} stage(s))")

    if all_errors:
        print(f"\nWarnings ({len(all_errors)}):")
        for e in all_errors:
            print(f"  {e}")


if __name__ == "__main__":
    main()
