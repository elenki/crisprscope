# CRISPRScope

[Work in progress `DRAFT`]

**CRISPRScope** is a bioinformatics toolkit designed to streamline the analysis of single-cell multiplex CRISPR screens. It can process raw targeted amplicon sequencing data to directly quantify genomic edits [TBD: `, resolving mutational co-occurrence and zygosity`] at single-cell resolution. The primary output is a standardized AnnData (`.h5ad`) object, integrating raw allele data with quality metrics and experimental metadata for robust downstream analysis and visualization.

---

## Background

Traditional high-content CRISPR screens often infer genomic edits indirectly from transcriptomic (RNA) data, which can be influenced by expression levels and other cellular states.

**CRISPRScope** is built to use targeted single-cell DNA sequencing to provide direct, unambiguous evidence of genomic edits. [TBD: `This allows for the precise quantification of:

- **Multiplexed Edits**: Which combinations of genes are edited in the same cell.
- **Zygosity**: Whether one (heterozygous) or both (homozygous) copies of a gene are edited.
- **Clonal Evolution**: How the frequency of these combinations changes under selective pressure.`]

---

## Features

- **End-to-End Pipeline**: Handles processing of raw sequencing data to produce a ready-to-use AnnData object.
- **Standardized Output**: Creates a multi-layered `AnnData` object that includes:
  - **Cell Metadata**: Information about each cell, such as cell type, treatment, and experimental conditions.
  - **Amplicon Metadata**: Details about the amplicons used in the experiment, including gene targets and primer sequences.
  - [TBD - ` **Zygosity Layer**: Indicates whether each allele is edited (1) or not (0), allowing for easy identification of heterozygous and homozygous edits.`]
  - [WIP - `**Raw Allele Sequences**: Provides the actual sequences of the alleles, enabling detailed analysis of the edits.`]
- **Command-Line Interface**: Simple, easy-to-use CLI for running the pipeline.

---

## Installation

### Clone the repository:

```bash
git clone https://github.com/elenki/crisprscope.git
cd crisprscope
````

### Create and activate a conda environment (recommended):

```bash
conda create -n crisprscope_env python=3.9
conda activate crisprscope_env
```

### Install CRISPRScope in editable mode:

This command will install the necessary dependencies and make the `crisprscope` command available in your terminal.

```bash
pip install -e .
```

---

## Usage

The primary function of CRISPRScope is the `build` command, which constructs the AnnData object.

```bash
crisprscope build --input <path_to_data_dir> --output <path_to_output_file.h5ad>
```

---

## Example

To run the pipeline on the example data provided in this project, use the following command from the project’s root directory:

```bash
crisprscope build \
  --input <path_to_data_dir> \
  --output crisprscope_output.h5ad \
  --verbose
```

This will process the input directory and create a file named `crisprscope_output.h5ad` in your current directory. This file can then be loaded into Python or R for downstream analysis.

---

## Loading the Output in Python

```python
import anndata as ad

# Load the object
adata = ad.read_h5ad("crisprscope_output.h5ad")

# Explore the data
print(adata)
print("\nCell Metadata (.obs):")
print(adata.obs.head())
print("\nAmplicon Metadata (.var):")
print(adata.var.head())
print("\nZygosity Layer (.layers['zygosity']):")
print(adata.layers['zygosity'][:5, :5])
```

---

## Project Status

The **Integrator module** (`build` command) is functionally complete.
Future work will include the implementation of the **Orchestrator module** (`run` command) to automate the execution of CRISPResso2 from raw FASTQ files.