"""Run PLIP on a PDB and locate its text report."""

from __future__ import annotations

import subprocess
from pathlib import Path


def run_plip(
    pdb_path: Path,
    output_dir: Path,
    executable: str = "plip",
) -> Path:
    """Run PLIP in text-report mode and return the resulting report.txt path.

    PLIP CLI options used here are: -f input PDB, -o output directory, -t text.
    If your installed PLIP version uses different flags, adjust this command.
    """
    pdb_path = Path(pdb_path).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not pdb_path.is_file():
        raise FileNotFoundError(f"PDB not found: {pdb_path}")

    command = [
        executable,
        "-f", str(pdb_path),
        "-o", str(output_dir),
        "-t",
        "--peptides", "X"
    ]
    # plip -f input.pdb -o output -t
    
    print("[PLIP]", " ".join(command))
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        check=False,
    )

    if completed.stdout:
        print(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError(
            f"PLIP exited with code {completed.returncode}.\n"
            f"STDERR:\n{completed.stderr}"
        )

    # PLIP versions/builds can place report.txt in a nested output directory.
    candidates = sorted(output_dir.rglob("report.txt"))
    if not candidates:
        candidates = sorted(output_dir.rglob("*.txt"))
    if not candidates:
        raise FileNotFoundError(
            f"PLIP completed but no text report was found under {output_dir}. "
            "Check the PLIP version/output options."
        )

    # Prefer a file explicitly named report.txt.
    report = next((p for p in candidates if p.name.lower() == "report.txt"), candidates[0])
    print(f"[INFO] PLIP report: {report}")
    return report
