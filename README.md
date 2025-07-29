# CRISPRScope

**CRISPRScope** is a bioinformatics toolkit for the analysis of single-cell multiplex CRISPR screens. It processes targeted single-cell dna sequencing data to directly quantify genomic edits, resolving mutational co-occurrence and zygosity at single-cell resolution. The primary output is a standardized `anndata` (`.h5ad`) object, integrating raw allele data with quality metrics and experimental metadata for robust, vectorized downstream analysis.

-----

## Background

High-content crispr screens are a powerful tool for functional genomics. However, traditional methods often infer genomic edits indirectly from transcriptomic (rna) data, which can be confounded by expression levels and other cellular states.

**crisprscope** is built on the principle of **direct genomic evidence**. By using targeted single-cell dna sequencing data, it provides an unambiguous readout of the genetic perturbations in each cell. This allows for the precise, quantitative analysis of complex genetic experiments, including:

  * **Multiplexed Edits**: Systematically identifying which combinations of genes are edited in the same cell.
  * **Zygosity**: Classifying edits as heterozygous (`wt/mut`), homozygous (`mut/mut`), or compound heterozygous (`mut/mut2`).
  * **Clonal Evolution**: Tracking the frequency of specific alleles and combinatorial genotypes under selective pressure.

-----

## The `AnnData` Output Structure

The core output of `crisprscope` is a multi-layered `anndata` object designed for computational analysis. This structure organizes the complex data from a single-cell crispr experiment into a single, portable `.h5ad` file.

  * `adata.X` [`float32`]: The primary data matrix containing the **modification percentage** for each cell (obs) and amplicon (var).
  * `adata.obs` (`pd.DataFrame`): Cell metadata, including quality scores (`amplicon score`, `read count`) and other per-cell metrics.
  * `adata.var` (`pd.DataFrame`): Amplicon metadata, including the reference sequence, guide information, and run-wide qc statistics (`editing_efficiency`, `total_coverage`).
  * `adata.layers['counts']` [`int32`]: A matrix of the **total read counts** for each cell and amplicon.
  * `adata.layers['zygosity']` [`int8`]: An integer-coded matrix classifying the edit state. The encoding is stored in `adata.uns['zygosity_encoding']`.
  * `adata.layers['allele_seq_1']` / `['allele_seq_2']` [`bytes`]: Matrices containing the raw dna sequence of the first and second most abundant alleles.
  * `adata.layers['allele_freq_1']` / `['allele_freq_2']` [`float32`]: Matrices containing the corresponding frequencies of the top two alleles.
  * `adata.uns` (`dict`): Unstructured metadata, including the run configuration (`crisprscope_config`) and the zygosity encoding map.

-----

## Installation

First, clone the repository and navigate into the directory:

```bash
git clone https://github.com/elenki/crisprscope.git
cd crisprscope
```

Next, create and activate a conda environment (recommended):

```bash
conda create -n crisprscope_env python=3.9
conda activate crisprscope_env
```

Finally, install `crisprscope` and its dependencies:

```bash
# For general use
pip install .

# For development (editable install)
pip install -e .
```

-----

## Usage

`crisprscope` provides a command-line interface for its two main modules.

### Integrator (`build` command)

The primary function is `build`, which constructs the `anndata` object from a completed crispr-esso2 analysis directory.

```bash
crisprscope build --input <path_to_data_dir> \
                  --output <path_to_output.h5ad> \
                  --config <path_to_config.yaml> \
                  --verbose
```

**Arguments:**

  * `--input`: (Required) Path to the input directory containing `settings.txt` and the `settings.txt.crispresso.filtered` subdirectory.
  * `--output`: (Required) Path for the output `.h5ad` file.
  * `--config`: (Optional) Path to a custom `config.yaml` file. If not provided, the package default is used.
  * `--verbose`: Enable detailed logging output.

### Orchestrator (`run` command)

*(Coming soon)* The `run` command will orchestrate the entire pipeline, starting from raw `.fastq` files to automate demultiplexing and execution of crispresso2.

-----

## Example Analysis in Python

The structured `anndata` object enables powerful and fast queries.

```python
import anndata as ad
import pandas as pd
import numpy as np

# Load the object
adata = ad.read_h5ad("crisprscope_output.h5ad")

# --- Example 1: Find all cells with a heterozygous edit in 'trp53' ---
trp53_idx = adata.var_names.get_loc('trp53')
het_code = 1 # From adata.uns['zygosity_encoding']

# A fast, vectorized boolean mask
is_het_trp53 = adata.layers['zygosity'][:, trp53_idx] == het_code
het_cells = adata.obs_names[is_het_trp53]

print(f"Found {len(het_cells)} cells with heterozygous trp53 edits.")

# --- Example 2: Find the most common edited allele for 'chd2' ---
chd2_idx = adata.var_names.get_loc('chd2')
ref_seq = adata.var.loc['chd2']['sequence']

# Get all top alleles for chd2, excluding the reference sequence
top_alleles = adata.layers['allele_seq_1'][:, chd2_idx]
edited_alleles = [a.decode('utf-8') for a in top_alleles if a.decode('utf-8') != ref_seq]

if edited_alleles:
    most_common_edit = pd.Series(edited_alleles).value_counts().idxmax()
    print(f"The most common edited allele for chd2 is: {most_common_edit}")

```

-----

## Configuration

Methodological parameters (e.g., zygosity thresholds, read structure) are controlled by a `config.yaml` file. A default is provided with the package, but you can supply your own using the `--config` flag.

-----

## Project Status

  * ✅ **Integrator Module (`build`)**: The core logic is complete, tested, and ready for use.
  * ⏳ **Performance & Scalability**: Phase 2 is underway to benchmark and optimize the integrator for very large datasets.
  * ⚙️ **Orchestrator Module (`run`)**: In development.