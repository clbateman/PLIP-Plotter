"""Parse a PLIP report and plot ligand-position/protein-residue interactions."""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


INTERACTION_COLORS = {
    "Hydrogen Bond": "#4A90D9",
    "Hydrophobic": "#E8A838",
    "Salt Bridge": "#6DBF67",
    "Pi-Stacking": "#B07FD4",
    "Water Bridge": "#5BC8C8",
    "Halogen Bond": "#E06C6C",
    "Metal Complex": "#A0A0A0",
}
INTERACTION_PATTERNS = {
    "Hydrogen Bond": re.compile(r"\*\*Hydrogen Bonds\*\*", re.I),
    "Hydrophobic": re.compile(r"\*\*Hydrophobic Interactions\*\*", re.I),
    "Salt Bridge": re.compile(r"\*\*Salt Bridges\*\*", re.I),
    "Pi-Stacking": re.compile(r"\*\*Pi(?:-| )Stacking\*\*", re.I),
    "Water Bridge": re.compile(r"\*\*Water Bridges\*\*", re.I),
    "Halogen Bond": re.compile(r"\*\*Halogen Bonds\*\*", re.I),
    "Metal Complex": re.compile(r"\*\*Metal Complexes\*\*", re.I),
}


# def parse_plip_file(
#     filepath: Path,
#     chain_aliases: dict[str, str] | None = None,
# ) -> dict:
#     """Return {ligand_position: {(resnr, restype, chain): set(interaction_types)}}."""
#     filepath = Path(filepath)
#     if not filepath.is_file():
#         raise FileNotFoundError(f"PLIP report not found: {filepath}")
#     aliases = {"D": "A", "E": "B"}
#     if chain_aliases:
#         aliases.update(chain_aliases)

#     text = filepath.read_text(encoding="utf-8", errors="replace")
#     interactions = defaultdict(lambda: defaultdict(set))

#     # Reports contain one block per ligand, e.g. LIG:XXX:1 (LIG).
#     sections = re.split(r"LIG:(\w+):(\d+)\s+\(LIG\)", text)
#     for i in range(1, len(sections), 3):
#         if i + 2 >= len(sections):
#             break
#         ligand_position = int(sections[i + 1])
#         body = sections[i + 2]
#         current_type = None

#         for line in body.splitlines():
#             for interaction_type, pattern in INTERACTION_PATTERNS.items():
#                 if pattern.search(line):
#                     current_type = interaction_type
#                     break

#             if current_type is None or not line.strip().startswith("|"):
#                 continue

#             parts = [p.strip() for p in line.strip().strip("|").split("|")]
#             # PLIP markdown table columns start RESNR, RESTYPE, RESCHAIN,
#             # RESNR_LIG. First cell is index 0 after stripping edge pipes.
#             if len(parts) < 4:
#                 continue
#             try:
#                 residue_number = int(parts[0])
#                 residue_type = parts[1].upper()
#                 chain = aliases.get(parts[2].upper(), parts[2].upper())
#             except (ValueError, IndexError):
#                 continue

#             interactions[ligand_position][
#                 (residue_number, residue_type, chain)
#             ].add(current_type)

#     return dict(interactions)


def parse_plip_file(
    filepath: Path,
    chain_aliases: dict[str, str] | None = None,
    peptide_chain: str = "X",) -> dict:
    """
    Parse PLIP interaction tables.

    Supports both:
      - SMALLMOLECULE reports (LIG)
      - PEPTIDE reports (--peptides)

    Returns:
        {
            peptide_position: {
                (protein_residue_number, protein_residue_type, protein_chain):
                    {interaction_types}
            }
        }
    """
    filepath = Path(filepath)

    if not filepath.is_file():
        raise FileNotFoundError(f"PLIP report not found: {filepath}")

    aliases = {"D": "A", "E": "B"}
    if chain_aliases:
        aliases.update(chain_aliases)

    text = filepath.read_text(
        encoding="utf-8",
        errors="replace",)

    interactions = defaultdict(lambda: defaultdict(set))
    current_type = None

    for line in text.splitlines():

        # Identify the current interaction type
        for interaction_type, pattern in INTERACTION_PATTERNS.items():
            if pattern.search(line):
                current_type = interaction_type
                break

        # Only process PLIP table rows
        if current_type is None or not line.strip().startswith("|"):
            continue

        parts = [p.strip() for p in line.strip().strip("|").split("|")]

        # Required columns:
        # 0 RESNR
        # 1 RESTYPE
        # 2 RESCHAIN
        # 3 RESNR_LIG
        # 4 RESTYPE_LIG
        # 5 RESCHAIN_LIG
        if len(parts) < 6:
            continue

        try:
            protein_residue_number = int(parts[0])
            protein_residue_type = parts[1].upper()
            protein_chain = parts[2].upper()

            peptide_position = int(parts[3])
            ligand_chain = parts[5].upper()

        except (ValueError, IndexError):
            # Skip table headers and separator rows
            continue

        # Keep interactions involving the selected peptide chain
        if ligand_chain != peptide_chain.upper():
            continue

        # Apply MHC chain aliases
        protein_chain = aliases.get(protein_chain, protein_chain,)

        interactions[peptide_position][(protein_residue_number, protein_residue_type, protein_chain,)].add(current_type)

    return dict(interactions)


def plot_interactions(
    interactions: dict,
    report_path: Path,
    output: Path,
    chain_a: str = "A",
    chain_b: str = "B",
) -> Path:
    """Save the two-panel interaction map."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)

    ligand_positions = sorted(interactions)
    residues_by_chain = {chain_a: set(), chain_b: set()}
    for residue_dict in interactions.values():
        for residue_key in residue_dict:
            if residue_key[2] in residues_by_chain:
                residues_by_chain[residue_key[2]].add(residue_key)

    residues_by_chain = {
        chain: sorted(values, key=lambda r: (r[0], r[1]))
        for chain, values in residues_by_chain.items()
    }

    if not ligand_positions:
        print(f"[WARN] No ligand interactions found in {report_path}; map skipped.")
        return output

    n_lig = len(ligand_positions)
    max_res = max(1, *(len(v) for v in residues_by_chain.values()))
    fig_w = max(10, 2.5 + n_lig * 0.9)
    fig_h = max(6, max_res * 0.28 * 2)
    fig, axes = plt.subplots(2, 1, figsize=(fig_w, fig_h), sharex=True)
    fig.suptitle(f"PLIP interactions — {Path(report_path).name}", fontsize=14, weight="bold")

    ligand_index = {pos: i for i, pos in enumerate(ligand_positions)}

    for ax, chain in zip(axes, (chain_a, chain_b)):
        residues = residues_by_chain[chain]
        residue_index = {res: i for i, res in enumerate(residues)}
        ax.set_facecolor("#F8F8F8")

        for x in range(n_lig + 1):
            ax.axvline(x - 0.5, color="#E0E0E0", linewidth=0.5)
        for y in range(len(residues) + 1):
            ax.axhline(y - 0.5, color="#E0E0E0", linewidth=0.5)

        for ligand_pos, residue_dict in interactions.items():
            x = ligand_index[ligand_pos]
            for residue_key, interaction_types in residue_dict.items():
                if residue_key[2] != chain or residue_key not in residue_index:
                    continue
                y = residue_index[residue_key]
                types = sorted(interaction_types)
                width = 0.8 / max(1, len(types))
                for j, interaction_type in enumerate(types):
                    offset = (j - (len(types) - 1) / 2) * width
                    rect = mpatches.FancyBboxPatch(
                        (x + offset - width / 2, y - 0.35),
                        width, 0.70,
                        boxstyle="round,pad=0.04",
                        facecolor=INTERACTION_COLORS.get(interaction_type, "#888888"),
                        edgecolor="white", linewidth=0.8,
                    )
                    ax.add_patch(rect)

        ax.set_xlim(-0.5, n_lig - 0.5)
        ax.set_ylim(-0.5, max(len(residues), 1) - 0.5)
        ax.set_yticks(range(len(residues)))
        ax.set_yticklabels([f"{r[1]}{r[0]}" for r in residues], fontsize=8)
        ax.set_title(f"Chain {chain}", loc="left", fontsize=11, weight="bold")
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", bottom=False)

    axes[-1].set_xticks(range(n_lig))
    axes[-1].set_xticklabels([f"Lig {p}" for p in ligand_positions], rotation=45, ha="right")
    axes[-1].set_xlabel("Ligand position")

    legend = [
        mpatches.Patch(facecolor=color, label=kind)
        for kind, color in INTERACTION_COLORS.items()
    ]
    fig.legend(handles=legend, title="Interaction type", loc="lower center",
               ncol=min(len(legend), 7), bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    fig.savefig(output, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[INFO] Interaction map saved: {output}")
    return output
