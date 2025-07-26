# CRISPRScope

**CRISPRScope** is a bioinformatics toolkit designed to streamline the analysis of single-cell multiplex CRISPR screens. It processes raw targeted amplicon sequencing data to directly quantify genomic edits, resolving mutational co-occurrence and zygosity at single-cell resolution. The primary output is a standardized AnnData (`.h5ad`) object, integrating raw allele data with quality metrics and experimental metadata for robust downstream analysis and visualization.

---

## The Core Innovation

Traditional high-content CRISPR screens often infer genomic edits indirectly from transcriptomic (RNA) data, which can be influenced by expression levels and other cellular states.

**CRISPRScope** is built on the methodological advance of using targeted single-cell DNA sequencing to provide direct, unambiguous evidence of genomic edits. This allows for the precise quantification of:

- **Multiplexed Edits**: Which combinations of genes are edited in the same cell.
- **Zygosity**: Whether one (heterozygous) or both (homozygous) copies of a gene are edited.
- **Clonal Evolution**: How the frequency of these combinations changes under selective pressure.

---

## Features

- **End-to-End Pipeline**: Handles data processing from CRISPResso2 outputs to an analysis-ready object.  
- **Standardized Output**: Creates a multi-layered `AnnData` object, the standard for single-cell analysis.  
- **Rich Data Layers**: Includes layers for modification percentage, read counts, zygosity status, and raw allele sequences.  
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