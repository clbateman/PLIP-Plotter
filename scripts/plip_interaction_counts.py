"""Count and plot PLIP interactions for one report."""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


INTERACTION_HEADERS = {
    "Hydrogen Bonds": "Hydrogen Bonds",
    "Hydrophobic Interactions": "Hydrophobic Interactions",
    "Salt Bridges": "Salt Bridges",
    "Pi Stacking": "Pi Stacking",
    "Pi-Stacking": "Pi Stacking",
    "Water Bridges": "Water Bridges",
    "Halogen Bonds": "Halogen Bonds",
    "Metal Complexes": "Metal Complexes",
}


def parse_plip_counts(
    report_path: Path,
    chain_a: str = "A",
    chain_b: str = "B",
    chain_aliases: dict[str, str] | None = None,
) -> dict:
    """Count interaction rows, ligand positions and unique protein residues."""
    report_path = Path(report_path)
    if not report_path.is_file():
        raise FileNotFoundError(f"PLIP report not found: {report_path}")

    aliases = {"D": "A", "E": "B"}
    if chain_aliases:
        aliases.update(chain_aliases)

    lines = report_path.read_text(encoding="utf-8", errors="replace").splitlines()
    total = 0
    chain_counts = {chain_a: 0, chain_b: 0}
    peptide_positions = set()
    interacting_residues = {chain_a: set(), chain_b: set()}
    interaction_type_counts = {}
    current_type = None

    for line in lines:
        heading_match = re.match(r"\s*\*\*(.+?)\*\*\s*$", line)
        if heading_match:
            heading = heading_match.group(1).strip()
            current_type = INTERACTION_HEADERS.get(heading, heading)
            continue

        if not line.strip().startswith("|") or current_type is None:
            continue

        parts = [part.strip() for part in line.strip().strip("|").split("|")]
        if len(parts) < 4:
            continue
        try:
            residue_number = int(parts[0])
            ligand_position = int(parts[3])
        except (ValueError, IndexError):
            continue

        chain_raw = parts[2].upper()
        chain = aliases.get(chain_raw, chain_raw)
        total += 1
        interaction_type_counts[current_type] = interaction_type_counts.get(current_type, 0) + 1

        if chain in chain_counts:
            chain_counts[chain] += 1
            peptide_positions.add(ligand_position)
            interacting_residues[chain].add((residue_number, parts[1].upper()))

    unique_a = len(interacting_residues[chain_a])
    unique_b = len(interacting_residues[chain_b])
    return {
        "total_interactions": total,
        "chain_A_interactions": chain_counts[chain_a],
        "chain_B_interactions": chain_counts[chain_b],
        "peptide_positions_involved": len(peptide_positions),
        "chain_A_unique_residues": unique_a,
        "chain_B_unique_residues": unique_b,
        "total_unique_residues": unique_a + unique_b,
        "interaction_types": interaction_type_counts,
    }


def plot_single_report(
    counts: dict,
    label: str,
    output_dir: Path,
    chain_a: str = "A",
    chain_b: str = "B",
) -> list[Path]:
    """Create the three single-complex count plots, adapted from the WT/mutant script."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    overall = pd.Series({
        "Total interactions": counts["total_interactions"],
        f"Chain {chain_a}": counts["chain_A_interactions"],
        f"Chain {chain_b}": counts["chain_B_interactions"],
        "Ligand positions involved": counts["peptide_positions_involved"],
    })
    fig, ax = plt.subplots(figsize=(9, 5))
    overall.plot(kind="bar", ax=ax, edgecolor="black")
    ax.set_ylabel("Count")
    ax.set_title(f"PLIP interaction counts — {label}")
    ax.tick_params(axis="x", rotation=20)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=9)
    fig.tight_layout()
    path = output_dir / f"{label}_interaction_counts.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    type_counts = counts["interaction_types"]
    if type_counts:
        series = pd.Series(type_counts).sort_index()
        fig, ax = plt.subplots(figsize=(max(8, len(series) * 1.4), 5))
        series.plot(kind="bar", ax=ax, edgecolor="black")
        ax.set_ylabel("Number of interaction rows")
        ax.set_title(f"PLIP interaction types — {label}")
        ax.tick_params(axis="x", rotation=30)
        ax.grid(axis="y", linestyle="--", alpha=0.35)
        ax.set_axisbelow(True)
        for container in ax.containers:
            ax.bar_label(container, padding=3, fontsize=9)
        fig.tight_layout()
        path = output_dir / f"{label}_interaction_types.png"
        fig.savefig(path, dpi=250, bbox_inches="tight")
        plt.close(fig)
        outputs.append(path)

    residue_counts = pd.Series({
        "All protein chains": counts["total_unique_residues"],
        f"Chain {chain_a}": counts["chain_A_unique_residues"],
        f"Chain {chain_b}": counts["chain_B_unique_residues"],
    })
    fig, ax = plt.subplots(figsize=(8, 5))
    residue_counts.plot(kind="bar", ax=ax, edgecolor="black")
    ax.set_ylabel("Number of unique protein residues")
    ax.set_title(f"Unique interacting protein residues — {label}")
    ax.tick_params(axis="x", rotation=15)
    ax.grid(axis="y", linestyle="--", alpha=0.35)
    ax.set_axisbelow(True)
    for container in ax.containers:
        ax.bar_label(container, padding=3, fontsize=9)
    fig.tight_layout()
    path = output_dir / f"{label}_interacting_residues.png"
    fig.savefig(path, dpi=250, bbox_inches="tight")
    plt.close(fig)
    outputs.append(path)

    for path in outputs:
        print(f"[INFO] Plot saved: {path}")
    return outputs
