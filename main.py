#!/usr/bin/env python3
"""Run PLIP and generate per-complex interaction plots for one PDB or a folder containing several PDBs."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from scripts.plip_analysis import run_plip
from scripts.plip_interaction_map import parse_plip_file, plot_interactions
from scripts.plip_interaction_counts import parse_plip_counts, plot_single_report


def find_pdb_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        if input_path.suffix.lower() not in {".pdb", ".ent"}:
            raise ValueError(f"Input file must be a PDB: {input_path}")
        return [input_path]
    if input_path.is_dir():
        pdbs = sorted(
            p for p in input_path.iterdir()
            if p.is_file() and p.suffix.lower() in {".pdb", ".ent"})
        
        if not pdbs:
            raise FileNotFoundError(f"No .pdb or .ent files found in {input_path}")
        return pdbs
    raise FileNotFoundError(f"Input path does not exist: {input_path}")

# def write_summary(rows: list[dict], output_dir: Path) -> None:  # this could be useful if folder contains a lot of pdbs
#     if not rows:
#         return
#     path = output_dir / "plip_summary.csv"
#     fieldnames = list(rows[0].keys())
#     with path.open("w", newline="", encoding="utf-8") as handle:
#         writer = csv.DictWriter(handle, fieldnames=fieldnames)
#         writer.writeheader()
#         writer.writerows(rows)
#     print(f"[INFO] Summary saved: {path}")

def main() -> int:
    parser = argparse.ArgumentParser(description="Run PLIP and generate interaction maps/count plots.")

    parser.add_argument(
        "--input", "-i", required=True, type=Path, 
        help="One PDB file or a folder containing PDB files.")
    
    parser.add_argument(
        "--output", "-o", type=Path, default=Path("output"), 
        help="Output directory (default: ./output).")
    
    parser.add_argument(
        "--plip-executable", default="plip", 
        help="PLIP executable name or path (default: plip).")
    
    parser.add_argument(
        "--chain-a", default="A", 
        help="Protein chain to count as chain A (default: A).")
    
    parser.add_argument(
        "--chain-b", default="B", 
        help="Protein chain to count as chain B (default: B).")
    
    parser.add_argument(   # this is useful for receptors with different chain IDs, e.g. receptor-ligand from a larger complex
        "--chain-aliases", default="D:A,E:B", 
        help="Comma-separated PLIP-chain:analysis-chain aliases (default: D:A,E:B).")
    
    parser.add_argument(
        "--skip-plip", action="store_true", 
        help="Do not run PLIP, use existing report.txt files in output folders.")
    
    args = parser.parse_args()

    aliases = {}
    try:
        for pair in args.chain_aliases.split(","):
            if pair.strip():
                source, target = pair.split(":")
                aliases[source.strip()] = target.strip()
    except ValueError:
        parser.error("--chain-aliases must look like D:A,E:B")

    pdb_files = find_pdb_files(args.input)
    args.output.mkdir(parents=True, exist_ok=True)  # Create output directory if it doesn't exist
    failures = []

    for pdb in pdb_files:
        sample_name = pdb.stem
        sample_out = args.output / sample_name
        plip_out = sample_out / "plip"
        plots_out = sample_out / "plots"
        plip_out.mkdir(parents=True, exist_ok=True)
        plots_out.mkdir(parents=True, exist_ok=True)

        print(f"\n{'=' * 68}\n[INFO] Processing {pdb.name}\n{'=' * 68}")

        try:
            if not args.skip_plip:
                report = run_plip(
                    pdb_path=pdb, output_dir=plip_out, executable=args.plip_executable)
                
            else:
                report = plip_out / "report.txt"
                if not report.is_file():
                    raise FileNotFoundError(
                        f"--skip-plip was used, but report is missing: {report}")

            interactions = parse_plip_file(report, chain_aliases=aliases)
            map_path = plots_out / f"{sample_name}_interaction_map.png"
            plot_interactions(
                interactions, report, map_path, chain_a=args.chain_a, chain_b=args.chain_b)

            counts = parse_plip_counts(
                report, chain_a=args.chain_a, chain_b=args.chain_b, chain_aliases=aliases)
            
            plot_single_report(
                counts, sample_name, plots_out,
                chain_a=args.chain_a, chain_b=args.chain_b,)

        except Exception as exc:
            failures.append((sample_name, str(exc)))
            print(f"[ERROR] {sample_name}: {exc}", file=sys.stderr)


    # write_summary(summary_rows, args.output)                             # useful if folder contains a lot of pdbs, might use later
    # alascan_summary = write_alascan_summary(alascan_rows, args.output)
    # if alascan_summary:
    #     print(f"[INFO] BAlaS summary saved: {alascan_summary}")

    if failures:
        print("\nSome complexes could not be processed:")
        for name, error in failures:
            print(f" - {name}: {error}")
        return 1

    print(f"\n[INFO] Finished :) Results are in: {args.output.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
