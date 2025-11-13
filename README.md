#  postgwas-tools

**postgwas-tools** is a Python package that provides a collection of tools for **post-GWAS (Genome-Wide Association Study) analysis**, annotation, and visualization.  
It simplifies steps such as plotting, SNP annotation, and result parsing for large-scale GWAS projects.

---

## Installation

### Option 1 — From source
```bash
git clone https://github.com/neurospin/postgwas-tools.git
cd postgwas-tools
python3 -m venv venv
. venv/bin/activate
pip install -e .
```

### Option 2 — From PyPI (planned)
```bash
pip install postgwas-tools
```

### Option 3 - remote install from Neurospin to irene (preferred)
To install the phecovbabel pacakge as a pcocc-rs image on irene!
 - install the installing tools from https://github.com/neurospin-brainomics/InstallerTGCC
 - run the following commands

```bash
git clone https://github.com/neurospin/postgwas-tools.git
cd postgwas-tools
# push-docker-tgcc.sh is a command from InstallTGCC
push-docker-tgcc.sh -s `pwd` -d postgwas -t 1.0 -r n4h00001rs
# connect TGCC
pcocc-rs run n4h00001rs:postgwas-1.0 -h
```

---

## Package Structure

```
src/postgwas_tools/
├── annot/
│   ├── fuma/                 # FUMA-related post-GWAS visualization
│   ├── ldsc/                 # LD Score Regression utilities
│   ├── magma/                # MAGMA annotation helpers
│   ├── replication/          # Replication and cross-cohort checks
│   └── utils.py              # Common helper functions
└── mostest/                  # MOSTest-specific result processors
```

---

## Documentation & Workflows

For detailed usage examples, including **QQ plots**, **locus-based analysis**, and **replication workflows**, see:
 [WORKFLOW.md](WORKFLOW.md)

---

## Dependencies

- Python ≥ 3.9  
- pandas, numpy, matplotlib, seaborn, scipy, qmplot, h5py

All dependencies are automatically installed when you run `pip install -e .`.

---

## Citation

If you use `postgwas-tools` in your work, please cite:

> Dufournet A. *postgwas-tools: a Python suite for streamlined post-GWAS analysis*. CEA NeuroSpin, 2025.

---

## Links

- [LICENSE](LICENSE)  
- [Homepage](https://github.com/Dufouranto0/postgwas-tools)  
- [Issue Tracker](https://github.com/Dufouranto0/postgwas-tools/issues)