"""
CRISPRScope: AnnData Builder

This module contains the CRISPRScopeAnnDataBuilder class, which is responsible
for taking parsed data and assembling it into a structured AnnData object.
"""

import pandas as pd
import numpy as np
from anndata import AnnData
from typing import Dict, Any, Tuple
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
        
        obs_df = pd.DataFrame(index=self.cell_barcodes)
        obs_df = obs_df.join(self.quality_scores)
        
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
        
        var_df = self.amplicons.reindex(self.amplicon_names)
        
        var_df['total_coverage'] = 0
        var_df['average_coverage'] = 0.0
        var_df['editing_efficiency'] = 0.0
        
        for amp in self.amplicon_names:
            count_col = f'totCount.{amp}'
            mod_col = f'modPct.{amp}'
            
            if count_col in self.editing_summary.columns:
                counts = self.editing_summary[count_col].dropna()
                var_df.loc[amp, 'total_coverage'] = int(counts.sum())
                var_df.loc[amp, 'average_coverage'] = counts.mean()

            if mod_col in self.editing_summary.columns:
                mod_pcts = self.editing_summary[mod_col].dropna()
                var_df.loc[amp, 'editing_efficiency'] = mod_pcts.mean()

        logger.info("Finished building .var DataFrame.")
        return var_df

    def _build_data_matrices(self) -> Tuple[np.ndarray, np.ndarray]:
        """
        Constructs the core data matrices (X and layers) from the editing summary.

        Returns:
            A tuple containing:
            - mod_pct_matrix (np.ndarray): Matrix for adata.X (n_obs x n_vars).
            - counts_matrix (np.ndarray): Matrix for adata.layers['counts'] (n_obs x n_vars).
        """
        logger.info("Building core data matrices (X, layers)...")
        
        # Create ordered lists of the column names we need to extract
        mod_pct_cols = [f'modPct.{amp}' for amp in self.amplicon_names]
        counts_cols = [f'totCount.{amp}' for amp in self.amplicon_names]

        # Select the data in the correct order.
        # .reindex ensures that if a column is missing, it's created with NaNs.
        mod_pct_df = self.editing_summary.reindex(columns=mod_pct_cols)
        counts_df = self.editing_summary.reindex(columns=counts_cols)

        # Convert to NumPy arrays and fill missing values with 0
        mod_pct_matrix = mod_pct_df.fillna(0).to_numpy(dtype=np.float32)
        counts_matrix = counts_df.fillna(0).to_numpy(dtype=np.int32)

        logger.info("Finished building core data matrices.")
        return mod_pct_matrix, counts_matrix

    def build(self) -> AnnData:
        """
        Builds the final AnnData object.

        Returns:
            The fully constructed AnnData object.
        """
        logger.info("Assembling final AnnData object...")
        
        obs_df = self._build_obs()
        var_df = self._build_var()
        mod_pct_matrix, counts_matrix = self._build_data_matrices()
        
        adata = AnnData(
            X=mod_pct_matrix,
            obs=obs_df,
            var=var_df
        )

        # Add the counts matrix to the layers
        adata.layers['counts'] = counts_matrix
        
        # Store the run settings in the .uns slot
        adata.uns['crisprscope_settings'] = self.settings
        
        logger.info("AnnData object created successfully with populated data layers.")
        return adata