<div align="center">

<img src="docs/hoodini_logo_github.svg" alt="Hoodini Logo" width="400"/>

### Large-scale gene neighborhood analyses that feel like magic

[![Python](https://img.shields.io/badge/Python-3.10+-blue?logo=python&logoColor=white)](https://www.python.org/)
[![Conda](https://img.shields.io/badge/Conda-coming_soon-lightgrey?logo=anaconda&logoColor=white)](https://bioconda.github.io/)
[![PyPI](https://img.shields.io/badge/PyPI-coming_soon-lightgrey?logo=pypi&logoColor=white)](https://pypi.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

[**🎮 Live Demo**](https://pentamorfico.github.io/hoodini-docs/demo) · [**📖 Documentation**](https://pentamorfico.github.io/hoodini-docs/docs/hoodini) · [**🖼️ Gallery**](https://pentamorfico.github.io/hoodini-docs/gallery) · [**🧪 Colab**](https://colab.research.google.com/github/pentamorfico/hoodini-colab/blob/main/hoodini_colab.ipynb)

</div>

---

## 🧬 What is Hoodini?

**Hoodini** is a gene-centric comparative genomics toolkit that fetches public assemblies, extracts gene neighborhoods, runs pairwise protein and nucleotide comparisons, annotates neighborhoods with defense systems and mobile elements, and builds phylogenetic trees — all with GPU-accelerated interactive visualization.

<br/>

<div align="center">

| 🚀 **Scales** | ⚡ **Fast** | 🔬 **Annotations** | 🎨 **Visualization** |
|:---:|:---:|:---:|:---:|
| 1000s of genomes | Minutes not hours | PADLOC, DefenseFinder, CCTyper | Publication-ready SVG |

</div>

<br/>

<div align="center">
  <img src="docs/hoodini-viz-export - 2026-01-14T051133.869.svg" alt="Example visualization" width="90%"/>
</div>

<br/>

---

## ✨ Key Features

- 📥 **Automated data retrieval** — Fetches assemblies from NCBI using protein or nucleotide accessions
- 🧬 **Neighborhood extraction** — Configurable genomic windows around target genes
- 🔗 **Protein clustering** — Groups homologous proteins for synteny comparison
- 📊 **Pairwise comparisons** — AAI (amino acid) and ANI (nucleotide) similarities
- 🌳 **Tree construction** — Phylogenetic trees from sequence identity
- 🛡️ **Defense annotations** — PADLOC, DefenseFinder, CCTyper, geNomad
- 🎨 **Interactive visualization** — Self-contained HTML with 50+ color palettes

---

## 🚀 Quick Start

```bash
# Single protein query
hoodini run --input WP_012345678.1 --output results

# With protein comparisons and phylogenetic tree
hoodini run --input proteins.txt --output results --prot-links --tree-mode aai_tree

# Full analysis with annotations
hoodini run --input proteins.txt --output results \
  --prot-links --tree-mode aai_tree \
  --padloc --deffinder --cctyper --genomad \
  --num-threads 16
```

📖 See the [Tutorial](https://pentamorfico.github.io/hoodini-docs/docs/hoodini/tutorial) for a complete walkthrough.

---

## 📦 Installation

Hoodini requires Python packages and bioinformatics tools. The recommended methods handle all dependencies.

> ⚠️ **Note:** Bioconda and PyPI packages are coming soon. Use the development installation below.

<table>
<tr>
<td width="50%">

**Using pixi (recommended)**

```bash
git clone https://github.com/pentamorfico/hoodini.git
cd hoodini
pixi install
pixi run hoodini download databases
```

</td>
<td width="50%">

**Using mamba/conda**

```bash
git clone https://github.com/pentamorfico/hoodini.git
cd hoodini
mamba env create -f environment.yml
mamba activate hoodini
pip install -e .
hoodini download databases
```

</td>
</tr>
</table>

<details>
<summary><strong>Python-only installation (uv/pip)</strong></summary>

<br/>

> ⚠️ This only installs Python packages. Bioinformatics tools must be in your PATH.

**Using uv:**

```bash
git clone https://github.com/pentamorfico/hoodini.git
cd hoodini
uv sync
uv run hoodini download databases
```

**Using pip:**

```bash
git clone https://github.com/pentamorfico/hoodini.git
cd hoodini
pip install -e .
hoodini download databases
```

</details>

<details>
<summary><strong>Docker</strong></summary>

<br/>

> ⚠️ Docker image available but not fully tested. Please report any issues.

```bash
docker volume create hoodini-data
docker run --rm -v hoodini-data:/app/src/hoodini/data \
  pentamorfico/hoodini:latest hoodini download databases
docker run --rm -v hoodini-data:/app/src/hoodini/data -v $(pwd):/work \
  pentamorfico/hoodini:latest hoodini run --input /work/proteins.txt --output /work/results
```

</details>

📖 See [Installation Guide](https://pentamorfico.github.io/hoodini-docs/docs/hoodini/installation) for detailed instructions.

---

## 🛠️ Usage

```
hoodini run      Run the main pipeline
hoodini download Download required databases
```

<details>
<summary><strong>Input Options</strong></summary>

| Option | Description |
|--------|-------------|
| `--input ID\|FILE` | Single accession or file with one per line |
| `--inputsheet FILE` | TSV with accessions and custom metadata |

</details>

<details>
<summary><strong>Neighborhood Extraction</strong></summary>

| Option | Description |
|--------|-------------|
| `--win-mode` | `win_genes` (gene count) or `win_nts` (nucleotides) |
| `--win INT` | Window size (default: 10 genes or 10000 nt) |
| `--sorfs` | Re-annotate small ORFs |

</details>

<details>
<summary><strong>Comparisons & Trees</strong></summary>

| Option | Description |
|--------|-------------|
| `--prot-links` | All-vs-all protein similarities |
| `--nt-links` | Pairwise nucleotide alignments |
| `--tree-mode` | `aai_tree` or `ani_tree` |

</details>

<details>
<summary><strong>Annotations</strong></summary>

| Option | Description |
|--------|-------------|
| `--padloc` | Defense systems (PADLOC) |
| `--deffinder` | Defense systems (DefenseFinder) |
| `--cctyper` | CRISPR-Cas typing |
| `--genomad` | Mobile genetic elements |
| `--domains LIST` | Domain databases |

</details>

📖 Full reference: [CLI Documentation](https://pentamorfico.github.io/hoodini-docs/docs/hoodini/cli-reference)

---

## 📁 Output

Hoodini generates a `hoodini-viz/` folder with:
- Self-contained HTML viewer
- Newick tree
- TSV and Parquet data files

📖 See [Outputs Guide](https://pentamorfico.github.io/hoodini-docs/docs/hoodini/outputs) for details.

---

## 📚 Learn More

| Resource | Description |
|----------|-------------|
| [📖 Documentation](https://pentamorfico.github.io/hoodini-docs/docs/hoodini) | Full documentation |
| [🎮 Live Demo](https://pentamorfico.github.io/hoodini-docs/demo) | Interactive examples |
| [🖼️ Gallery](https://pentamorfico.github.io/hoodini-docs/gallery) | Real-world examples from publications |
| [🧪 Colab](https://colab.research.google.com/github/pentamorfico/hoodini-colab/blob/main/hoodini_colab.ipynb) | Run in Google Colab |
| [📦 hoodini-viz](https://github.com/pentamorfico/hoodini-viz) | Visualization library (npm) |

---

## 🙏 Acknowledgments

Hoodini is inspired by excellent tools in the field:

- [GCsnap](https://github.com/JoanaMPereira/GCsnap) — Gene context visualization
- [FlaGs](https://github.com/GCA-VH-lab/FlaGs) — Flanking genes analysis
- [Taxonium](https://github.com/theosanderson/taxonium) — Large trees visualization
- [clinker](https://github.com/gamcil/clinker) — Gene cluster comparison
- [gggenes](https://github.com/wilkox/gggenes) — Gene arrow maps in R
- [gggenomes](https://github.com/thackl/gggenomes) — Comparative genomics visualization

---

## 📄 Citation

> [Citation pending publication]

## 📜 License

MIT License. See [LICENSE](LICENSE) file.
