# PLIP-Plotter

A Python tool for parsing PLIP reports and visualising protein–ligand and protein–peptide interactions.

PLIP-Plotter runs [PLIP (Protein–Ligand Interaction Profiler)](https://github.com/pharmai/plip) on PDB structures and generates plots from the resulting interaction reports. It supports both small-molecule ligands and peptide ligands.

The aim is to make PLIP interaction results easier to explore and interpret through visual representations.

## PLIP

PLIP-Plotter uses [Protein–Ligand Interaction Profiler (PLIP)](https://github.com/pharmai/plip) for interaction detection.

PLIP is also available through the [LIP Web Server](https://plip-tool.biotec.tu-dresden.de/).

More information on PLIP:

- Salentin, S. et al. (2015). *PLIP: fully automated protein–ligand interaction profiler*. Nucleic Acids Research, 43(W1), W443–W447. [https://doi.org/10.1093/nar/gkv315](https://doi.org/10.1093/nar/gkv315)
- Schake, P., Bolz, S. N. et al. (2025). *PLIP 2025: introducing protein–protein interactions to the protein–ligand interaction profiler*. Nucleic Acids Research. [https://doi.org/10.1093/nar/gkaf361](https://doi.org/10.1093/nar/gkaf361)

## Input

The input consists of one or more PDB files containing a receptor and a ligand.

Supported ligand types:
- Small molecules
- Peptides (using PLIP peptide mode)

The tool can process:
- A single PDB file
- A folder containing multiple PDB files

## Output

For each PDB structure, PLIP-Plotter creates a separate output folder under `results/`, named after the input file.

```text
results/
└── file_name/
    ├── plip/
    │   ├── file_name_report.txt
    │   └── other PLIP output files
    └── plots/
        ├── file_name_interaction_map.png
        ├── file_name_interaction_types.png
        └── file_name_interaction_counts.png
```

The PLIP output includes the interaction report and the processed (deprotonated) PDB structure.

The following plots are generated:

- **Interaction map:** Visualises interactions between ligand positions and receptor residues, with separate panels for the two receptor chains.
- **Interaction types:** Shows the different types of interactions detected by PLIP.
- **Interaction counts:** Shows the number of interactions, including unique interacting residues and total interactions per chain.

## Usage

### Process a single PDB file

```bash
python main.py -i path/to/structure.pdb
```

Example:

```bash
python main.py -i data/R22E_B9-23.pdb
```

### Process a folder containing multiple PDB files

```bash
python main.py -i path/to/folder/
```

Example:

```bash
python main.py -i data/
```

### Optional arguments

| Argument | Description |
| -------- | ----------- |
| `-o`, `--output` | Output directory. Default: `results/` |
| `--chain-a` | Name of the first receptor chain. Default: `A` |
| `--chain-b` | Name of the second receptor chain. Default: `B` |
| `--chain-aliases` | Defines chain aliases, useful when the same receptor chain has different identifiers across structures. Default: `A:D,B:E` |
| `--skip-plip` | Skips PLIP execution and uses existing PLIP reports in the output folders. |

### 1. Clone the repository

```bash
git clone https://github.com/clbateman/PLIP-Plotter.git
cd PLIP-Plotter
```

### 2. Create and activate a virtual environment

Using Conda:

```bash
conda create -n plip-plotter python=3.10
conda activate plip-plotter
```

### 3. Install Python dependencies

```bash
pip install -r requirements.txt
```

### 4. Install PLIP

PLIP and openbabel must also be installed and accessible from your command line. Follow the [official PLIP installation instructions](https://github.com/pharmai/plip).

`pip install plip`

For conda users: 
`conda install -c conda-forge openbabel`

## Requirements

- Python 3.10
- PLIP
- Openbabel
- Python dependencies listed in `requirements.txt`