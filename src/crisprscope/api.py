"""
CRISPRScope: High-level API for building and saving AnnData objects.
"""

import logging
from pathlib import Path
from typing import Dict, Any
import yaml
from importlib import resources
import shutil

from anndata import AnnData

from .integrator.loaders import (
    load_settings,
    load_amplicons,
    load_editing_summary,
    load_quality_scores,
    load_crispresso_alleles,
)
from .integrator.builder import CRISPRScopeAnnDataBuilder

logger = logging.getLogger(__name__)

def _load_config(config_path: str = None) -> Dict[str, Any]:
    """Loads a YAML config, using the package default if no path is provided."""
    if config_path:
        config_file = Path(config_path)
        logger.info(f"Loading user-provided config from: {config_file}")
        if not config_file.is_file():
            raise FileNotFoundError(f"Custom config file not found: {config_path}")
    else:
        logger.info("No custom config provided. Loading package default.")
        config_file = resources.files('crisprscope').joinpath('config.yaml')

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def build_anndata(base_path: str, output_path: str = None, config_path: str = None) -> AnnData:
    """
    Builds and optionally saves a CRISPRScope AnnData object from a
    CRISPResso2 output directory.
    """
    logger.info(f"--- Starting CRISPRScope Integrator ---")
    
    config = _load_config(config_path)
    logger.info(f"Processing data from: {base_path}")

    # --- 1. Load all data sources ---
    logger.info("Step 1/3: Loading primary data files and parsing alleles...")
    data_path = Path(base_path)
    allele_paths = [] # ensure this is defined
    try:
        settings = load_settings(data_path / "settings.txt")
        amplicon_file = Path(settings.get('amplicons'))
        amplicons_df = load_amplicons(amplicon_file)
        summary_df = load_editing_summary(data_path / "settings.txt.filteredEditingSummary.txt")
        scores_df = load_quality_scores(data_path / "settings.txt.amplicon_score.txt")
        # this now returns a list of paths to temporary parquet files
        allele_paths = load_crispresso_alleles(data_path / "settings.txt.crispresso.filtered")
        logger.info("✅ All data loaded successfully.")
    except Exception as e:
        logger.error(f"❌ Failed during data loading: {e}", exc_info=True)
        raise

    adata = None
    try:
        # --- 2. Build the AnnData object ---
        logger.info("Step 2/3: Building AnnData object...")
        builder = CRISPRScopeAnnDataBuilder(
            config=config,
            settings=settings,
            amplicons=amplicons_df,
            editing_summary=summary_df,
            quality_scores=scores_df,
            allele_parquet_paths=allele_paths # use the new keyword argument
        )
        adata = builder.build()
        logger.info("✅ AnnData object built successfully.")

        # --- 3. Save the object if an output path is provided ---
        if output_path:
            logger.info(f"Step 3/3: Saving AnnData object to: {output_path}")
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(output_file, compression="gzip")
            logger.info(f"✅ AnnData object saved successfully.")
    
    finally:
        # --- 4. cleanup ---
        # this `finally` block ensures we ALWAYS clean up the temp files,
        # even if the builder or save steps fail.
        if allele_paths:
            temp_dir_to_remove = allele_paths[0].parent
            logger.info(f"Cleaning up temporary directory: {temp_dir_to_remove}")
            try:
                shutil.rmtree(temp_dir_to_remove)
                logger.info("✅ Cleanup complete.")
            except Exception as e:
                logger.warning(f"⚠️ Could not clean up temporary directory {temp_dir_to_remove}: {e}")
            
    logger.info(f"--- CRISPRScope Integrator Finished ---")
    return adata