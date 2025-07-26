"""
CRISPRScope: AnnData Builder

This module contains the CRISPRScopeAnnDataBuilder class, which is responsible
for taking parsed data and assembling it into a structured AnnData object.
"""

import pandas as pd
import numpy as np
from anndata import AnnData
from typing import Dict, Any
import logging

from .utils import extract_amplicon_names

logger = logging.getLogger(__name__)

class CRISPRScopeAnnDataBuilder:
    """
    Assembles a CRISPRScope AnnData object from loaded data sources.
    """
    def __init__(self, settings: Dict[str, Any], amplicons: pd.DataFrame, 
                 editing_summary: pd.DataFrame, quality_scores: pd.DataFrame):
        """
        Initializes the builder with all necessary data.

        Args:
            settings: Parsed data from settings.txt.
            amplicons: Parsed data from amplicons.txt.
            editing_summary: Parsed data from filteredEditingSummary.txt.
            quality_scores: Parsed data from amplicon_score.txt.
        """
        logger.info("Initializing AnnData Builder...")
        self.settings = settings
        self.amplicons = amplicons
        self.editing_summary = editing_summary
        self.quality_scores = quality_scores

        # Establish the core dimensions of the AnnData object
        self.cell_barcodes = self.editing_summary.index.tolist()
        self.amplicon_names = extract_amplicon_names(self.editing_summary)
        
        self.n_obs = len(self.cell_barcodes)
        self.n_vars = len(self.amplicon_names)
        logger.info(f"Builder configured for {self.n_obs} cells and {self.n_vars} amplicons.")

    def _build_obs(self) -> pd.DataFrame:
        """Constructs the .obs (cell metadata) DataFrame."""
        logger.info("Building .obs DataFrame (cell metadata)...")
        
        # Start with the index we already have
        obs_df = pd.DataFrame(index=self.cell_barcodes)
        
        # Join the quality scores data. We use .reindex to ensure the order is correct.
        obs_df = obs_df.join(self.quality_scores)
        
        # Fill any missing values for cells that might not be in the scores file
        obs_df.fillna({
            'Amplicon Score': 0.0,
            'Read Count': 0,
            'Barcode Rank': -1,
            'Color': 'Unknown'
        }, inplace=True)

        logger.info("Finished building .obs DataFrame.")
        return obs_df

    def _build_var(self) -> pd.DataFrame:
        """Constructs the .var (amplicon metadata) DataFrame."""
        logger.info("Building .var DataFrame (amplicon metadata)...")
        
        # Use the amplicons DataFrame as the base, reindexed to our official list
        var_df = self.amplicons.reindex(self.amplicon_names)
        
        # We can pre-populate columns for statistics we'll calculate later
        var_df['total_coverage'] = 0
        var_df['average_coverage'] = 0.0
        var_df['editing_efficiency'] = 0.0
        
        # Calculate summary statistics for each amplicon
        for amp in self.amplicon_names:
            count_col = f'totCount.{amp}'
            mod_col = f'modPct.{amp}'
            
            if count_col in self.editing_summary.columns:
                counts = self.editing_summary[count_col].dropna()
                var_df.loc[amp, 'total_coverage'] = int(counts.sum())
                var_df.loc[amp, 'average_coverage'] = counts.mean()

            if mod_col in self.editing_summary.columns:
                 # Calculate efficiency only on cells with coverage
                mod_pcts = self.editing_summary[mod_col].dropna()
                var_df.loc[amp, 'editing_efficiency'] = mod_pcts.mean()

        logger.info("Finished building .var DataFrame.")
        return var_df

    def build(self) -> AnnData:
        """
        Builds the final AnnData object.

        Returns:
            The fully constructed AnnData object.
        """
        logger.info("Assembling final AnnData object...")
        
        obs_df = self._build_obs()
        var_df = self._build_var()

        # For now, the primary data matrix X will be a placeholder.
        # We will populate it and the layers in the next step.
        placeholder_x = np.zeros((self.n_obs, self.n_vars), dtype=np.float32)
        
        adata = AnnData(
            X=placeholder_x,
            obs=obs_df,
            var=var_df
        )

        # Store the run settings in the .uns slot
        adata.uns['crisprscope_settings'] = self.settings
        
        logger.info("AnnData object skeleton created successfully.")
        return adata
