"""
CRISPRScope: High-level API for building and saving AnnData objects.
"""

import logging
from pathlib import Path
from typing import Dict, Any
import yaml
from importlib import resources

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
    """
    Loads a YAML config, using the package default if no path is provided.
    
    This helper function centralizes the config loading logic.
    """
    if config_path:
        config_file = Path(config_path)
        logger.info(f"Loading user-provided config from: {config_file}")
        if not config_file.is_file():
            raise FileNotFoundError(f"Custom config file not found: {config_path}")
    else:
        # This is the modern, robust way to find a file inside your installed package.
        # It works correctly for wheels, eggs, and editable installs.
        logger.info("No custom config provided. Loading package default.")
        config_file = resources.files('crisprscope').joinpath('config.yaml')

    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def build_anndata(base_path: str, output_path: str = None, config_path: str = None) -> AnnData:
    """
    Builds and optionally saves a CRISPRScope AnnData object from a
    CRISPResso2 output directory.

    This is the main user-facing function for the Integrator module.

    Args:
        base_path: The path to the directory containing the 'settings.txt'
                   and 'settings.txt.crispresso.filtered' subdirectories.
        output_path: Optional. If provided, the path to save the final
                     .h5ad object.
        config_path: Optional. Path to a custom YAML config file. If not
                     provided, the package default is used.

    Returns:
        The fully constructed AnnData object.
    """
    logger.info(f"--- Starting CRISPRScope Integrator ---")
    
    # --- 0. Load Configuration ---
    try:
        config = _load_config(config_path)
    except Exception as e:
        logger.error(f"❌ Failed during configuration loading: {e}", exc_info=True)
        raise
    
    logger.info(f"Processing data from: {base_path}")

    # --- 1. Load all data sources ---
    logger.info("Step 1/2: Loading all data files...")
    data_path = Path(base_path)
    try:
        settings = load_settings(data_path / "settings.txt")
        amplicon_file = Path(settings.get('amplicons'))
        amplicons_df = load_amplicons(amplicon_file)
        summary_df = load_editing_summary(data_path / "settings.txt.filteredEditingSummary.txt")
        scores_df = load_quality_scores(data_path / "settings.txt.amplicon_score.txt")
        alleles_data = load_crispresso_alleles(data_path / "settings.txt.crispresso.filtered")
        logger.info("✅ All data loaded successfully.")
    except Exception as e:
        logger.error(f"❌ Failed during data loading: {e}", exc_info=True)
        raise

    # --- 2. Build the AnnData object ---
    logger.info("Step 2/2: Building AnnData object...")
    try:
        builder = CRISPRScopeAnnDataBuilder(
            config=config,
            settings=settings,
            amplicons=amplicons_df,
            editing_summary=summary_df,
            quality_scores=scores_df,
            alleles_data=alleles_data
        )
        adata = builder.build()
        logger.info("✅ AnnData object built successfully.")
    except Exception as e:
        logger.error(f"❌ Failed during AnnData construction: {e}", exc_info=True)
        raise

    # --- 3. Save the object if an output path is provided ---
    if output_path:
        logger.info(f"Saving AnnData object to: {output_path}")
        try:
            output_file = Path(output_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            adata.write_h5ad(output_file, compression="gzip")
            logger.info(f"✅ AnnData object saved successfully.")
        except Exception as e:
            logger.error(f"❌ Failed to save AnnData object: {e}", exc_info=True)
            raise
            
    logger.info(f"--- CRISPRScope Integrator Finished ---")
    return adata