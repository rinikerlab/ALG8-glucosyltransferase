# ALG8 Membrane-Embedded Mannosyltransferase Simulation and Analysis

(This project is an successive project from https://github.com/rinikerlab/ALG_mannosyltransferase)
## Overview
This repository contains molecular dynamics (MD) simulations and trajectory analyses of the **membrane-embedded glucosyltransferase ALG8**.  
The goal of this project is to understand the donor selectivity and key interactions of ALG8 during glycosylation.

### Scientific Questions
This study aims to answer three main questions:

1. **What makes the mutant D36N keep the product** after the catalysis?   
2. **What determines ALG8’s selectivity** toward the glucose donor compared to mannose? 

To address these:
- (b) Simulations with **mannose** and **gluvose** donors and were conducted to assess donor selectivity in the pre-transfer state.
- (c) Simulations with **D36 (wildtype)** and **D36N (mutant)** and were conducted to assess product stability in the post-transfer state.
- (c) All possible protonation state of H40 ($\delta$, $\epsilon$, and double) are studied to investigate the role of this amino acid.

---

## MD Simulations

Example MD simulation input files are provided in the folder [`MDsimulations/`](MDsimulations/).  
For quick testing, the demonstration protocol uses **50× fewer steps** to allow short test runs.

You can run the short test directly by executing:

```bash
bash MD.sh
```

This script performs **energy minimization, equilibration, and production** in sequence.

### Running Full-Length Simulations
To perform full-length production runs, replace the demonstration input files with the ones in  
[`MDsimulations/longer_input_files/`](MDsimulations/longer_input_files/).

### System Composition
| Component | Description |
|------------|-------------|
| Protein | ALG8 glocosyltransferase |
| Substrates | Acceptor and donor |
| Membrane | POPC bilayer (constructed via CHARMM-GUI) |
| Solvent | OPC water |
| Ions | 0.15 M NaCl |
| Force fields | `ff19SB` (protein), `Lipid21` (lipid), `GAFF2` (substrates), `opc` (water) |

### Simulation Environment
| Item | Description |
|------|-------------|
| MD Engine | AMBER 24 |
| GPU | RTX 3090 Ti (approx. 2 h → 10 ns) |
| OS | Linux (tested), should work on other platforms supporting AMBER and Python |
| Reproducibility | Test inputs are already minimized; random seeds are not fixed |

A full list of Python dependencies is provided in [`environment.yml`](environment.yml).

---

##  Trajectory Analysis

All trajectory analysis is performed in **Python** (Jupyter Notebook) using:
- `numpy` and `pytraj` for data processing  
- `matplotlib` and `seaborn` for visualization  

The analysis scripts automatically compute:

- RMSD  
- RMSF  
- Interatomic distances  
- Interatomic angles
- Hydrogen bonds

### Usage
By specifying the folder containing topology and trajectory files, all analyses can be executed with:

```python
bash multi_analysis.sh
```

Upon completion:
- **Figures** are saved in the [`Figures/`](Figures/) directory.  
- **Analyzed objects** (in `.pkl` format) are saved in the [`analysis/objects/`](objects/) directory for faster future access.

The `.pkl` files contain precomputed analysis results, allowing users to regenerate figures without reloading trajectories.

All plots included in the related publication are automatically produced during analysis.

---

##  Repository Structure
```
ALG9_project/
│
├── MDsimulations/
│   ├── input_files/                # Example minimized test inputs
│   ├── longer_input_files/         # Full-length production inputs
│   └── MD.sh                       # Run script for test simulations
│
├── analysis/
│   ├── ALG9_initial.ipynb          # Jupyter notebook for the analysis on the initial runs
│                                     with a guessed donor binding pose
│   ├── ALG9.ipynb                  # Jupyter notebook for the analysis on the later runs 
|                                     from the binding pose revealed from the initial runs
│   ├── objects/                    # Stored analyzed data (.pkl files)
│   └── Figures/                    # Auto-generated figures
│
├── environment.yml                 # Python environment specification
└── README.md                       # This file
```

---

## Notes for Users

- To analyze your own trajectories, modify the input paths in the notebook or analysis script to point to your topology and trajectory files.  
- The workflow has been tested on Linux and should be portable to macOS or Windows if AMBER and Python dependencies are available.  
- GPU acceleration (CUDA) is highly recommended for production-length simulations.
---

## Trajectories of the Published Work
Due to the large filesize (> 20GB without solvents), MD trajectories are provided upon request at sriniker@ethz.ch.

---

---


> **Citation:**  
> *To be added once the manuscript DOI becomes available.*

---