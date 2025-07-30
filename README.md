# CRISPRSCope

**CRISPRSCope** is a high-performance bioinformatics toolkit for analyzing single-cell multiplex CRISPR screens. It provides an end-to-end workflow, from raw sequencing data to an analysis-ready `AnnData` object. By using targeted single-cell **DNA sequencing**, CRISPRSCope directly quantifies genomic edits, precisely resolving mutational co-occurrence and zygosity at single-cell resolution.

-----

## Background

High-content CRISPR screens are a powerful tool for functional genomics. However, traditional methods often infer genomic edits indirectly from transcriptomic (RNA) data, which can be confounded by expression levels and other cellular states <sup>[1](https://pubmed.ncbi.nlm.nih.gov/37214176/)</sup>.

**CRISPRSCope** is built on the principle of **direct genomic evidence**. It provides an unambiguous readout of the genetic perturbations in each cell, enabling precise, quantitative analysis of complex genetic experiments, including:

  * **Multiplexed Edits**: Systematically identifying which combinations of genes are edited in the same cell <sup>[2](https://pubmed.ncbi.nlm.nih.gov/33081820/)</sup>.
  * **Zygosity**: Classifying edits as heterozygous (`WT/Mut`), homozygous (`Mut/Mut`), or compound heterozygous (`Mut/Mut2`).
  * **Clonal Evolution**: Tracking the frequency of specific alleles and combinatorial genotypes under selective pressure.

### References:
1. Bock C, Datlinger P, Chardon F, Coelho MA, Dong MB, Lawson KA, Lu T, Maroc L, Norman TM, Song B, Stanley G, Chen S, Garnett M, Li W, Moffat J, Qi LS, Shapiro RS, Shendure J, Weissman JS, Zhuang X. **High-content CRISPR screening**. Nat Rev Methods Primers. 2022;2(1):9. doi: 10.1038/s43586-022-00098-7. Epub 2022 Feb 10. PMID: 37214176; PMCID: PMC10200264.

<!-- https://pubmed.ncbi.nlm.nih.gov/33081820/ -->
2. Ten Hacken E, Clement K, Li S, Hernández-Sánchez M, Redd R, Wang S, Ruff D, Gruber M, Baranowski K, Jacob J, Flynn J, Jones KW, Neuberg D, Livak KJ, Pinello L, Wu CJ. **High throughput single-cell detection of multiplex CRISPR-edited gene modifications**. Genome Biol. 2020 Oct 20;21(1):266. doi: 10.1186/s13059-020-02174-1. PMID: 33081820; PMCID: PMC7574538.

<!-- https://pubmed.ncbi.nlm.nih.gov/36468984/ -->
<!-- 3. Ten Hacken E, Sewastianik T, Yin S, Hoffmann GB, Gruber M, Clement K, Penter L, Redd RA, Ruthen N, Hergalant S, Sholokhova A, Fell G, Parry EM, Broséus J, Guieze R, Lucas F, Hernández-Sánchez M, Baranowski K, Southard J, Joyal H, Billington L, Regis FFD, Witten E, Uduman M, Knisbacher BA, Li S, Lyu H, Vaisitti T, Deaglio S, Inghirami G, Feugier P, Stilgenbauer S, Tausch E, Davids MS, Getz G, Livak KJ, Bozic I, Neuberg DS, Carrasco RD, Wu CJ. **In Vivo Modeling of CLL Transformation to Richter Syndrome Reveals Convergent Evolutionary Paths and Therapeutic Vulnerabilities**. Blood Cancer Discov. 2023 Mar 1;4(2):150-169. doi: 10.1158/2643-3230.BCD-22-0082. PMID: 36468984; PMCID: PMC9975769. -->

-----

## Workflow

CRISPRSCope consists of two main modules that can be run sequentially or independently.

```
[Raw FASTQs (R1/R2)] ──> │ Orchestrator │ ──> [Per-Cell FASTQs] 
      .. ──> (CRISPResso2) ──> [CRISPResso2 Outputs] 
        .. ──> │ Integrator │ ──> [Analysis-Ready AnnData (.h5ad)]
```

1.  **Orchestrator (`run`)**: Processes raw, paired-end FASTQ files. Its high-performance demultiplexer sorts reads into cell-specific FASTQ files, preparing them for analysis. *(CRISPResso2 execution is under development)*.
2.  **Integrator (`build`)**: Parses the outputs from a completed CRISPResso2 run, aggregates data across all cells and amplicons, and builds the final, structured `AnnData` object.

-----

## The `AnnData` Output Structure

The core output of the `build` command is a multi-layered `AnnData` object designed for computational analysis.

  * `adata.X` [`float32`]: Primary data matrix containing the **modification percentage** for each cell (obs) and amplicon (var).
  * `adata.obs` (`pd.DataFrame`): Cell metadata, including quality scores (`Amplicon Score`, `Read Count`).
  * `adata.var` (`pd.DataFrame`): Amplicon metadata, including reference sequence and run-wide QC statistics (`editing_efficiency`, `total_coverage`).
  * `adata.layers['counts']` [`int32`]: Matrix of the **total read counts** for each cell and amplicon.
  * `adata.layers['zygosity']` [`int8`]: Integer-coded matrix classifying the edit state. The encoding is stored in `adata.uns['zygosity_encoding']`.
  * `adata.layers['allele_seq_1']` / `['allele_seq_2']` [`bytes`]: Matrices containing the raw DNA sequence of the top two most abundant alleles.
  * `adata.layers['allele_freq_1']` / `['allele_freq_2']` [`float32`]: Matrices containing the corresponding frequencies of the top two alleles.
  * `adata.uns` (`dict`): Unstructured metadata, including the run configuration (`crisprscope_config`).

-----

## Installation

Clone the repository and navigate into the directory:

```bash
git clone https://github.com/elenki/crisprscope.git
cd crisprscope
```

Create and activate a conda environment (recommended):

```bash
conda create -n crisprscope_env python=3.9
conda activate crisprscope_env
```

Install `crisprscope` and its dependencies (`anndata`, `pandas`, `pyyaml`, `dnaio`):

```bash
# For general use
pip install .

# For development (includes pytest, jupyter, matplotlib)
pip install -e .[dev]
```

-----

## Usage

CRISPRSCope provides a command-line interface for its two main modules.

### Integrator (`build` command)

Constructs the `AnnData` object from a completed CRISPResso2 analysis directory.

```bash
crisprscope build --input <path_to_data_dir> \
                  --output <path_to_output.h5ad> \
                  --config <path_to_config.yaml> \
                  --verbose
```

  * `--input`: (Required) Path to the directory containing `settings.txt`.
  * `--output`: (Required) Path for the output `.h5ad` file.
  * `--config`: (Optional) Path to a custom `config.yaml`. The package default is used otherwise.

### Orchestrator (`run` command)

Processes raw sequencing data. The high-performance demultiplexer is complete.

```bash
crisprscope run --settings <path_to_settings.txt> \
                --config <path_to_config.yaml> \
                --threads <int> \
                --verbose
```

  * `--settings`: (Required) Path to the run's `settings.txt` file, which specifies input FASTQs and barcode files.
  * `--threads`: (Optional) Number of parallel processes to use for future steps.

-----

## Configuration

Methodological parameters (e.g., zygosity thresholds, read structure, performance tuning) are controlled by a `config.yaml` file. A default is provided with the package, and you can supply your own using the `--config` flag.

```yaml
# Default config.yaml structure
analysis_parameters:
  zygosity:
    wt_max_mod_pct: 20
    # ...
read_structure:
  mission_bio_v2:
    barcode_1_start: 5
    # ...
performance_parameters:
  demultiplexer_chunk_size: 1000000
```

-----

## Project Status

  * ✅ **Integrator Module (`build`)**: The core logic is complete, tested, and performance-validated to scale linearly with cell count.
  * ✅ **Orchestrator Module (`run`)**: The high-performance, scalable demultiplexing step is complete and tested.
  * ⚙️ **Next Steps**: Implementation of the remaining orchestrator steps: parallelized CRISPResso2 execution and results summarization.