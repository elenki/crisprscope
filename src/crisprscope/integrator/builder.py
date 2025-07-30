"""
CRISPRScope: AnnData Builder

This module contains the CRISPRScopeAnnDataBuilder class, which is responsible
for taking parsed data and assembling it into a structured AnnData object
using efficient, vectorized operations and iterative loading for scalability.
"""

import pandas as pd
import numpy as np
from anndata import AnnData
from typing import Dict, Any, List
from pathlib import Path
import logging

from .utils import extract_amplicon_names

logger = logging.getLogger(__name__)

class CRISPRScopeAnnDataBuilder:
    """
    Assembles a CRISPRScope AnnData object from loaded data sources.
    """
    def __init__(self, config: Dict[str, Any], settings: Dict[str, Any],
                 amplicons: pd.DataFrame, editing_summary: pd.DataFrame,
                 quality_scores: pd.DataFrame, allele_parquet_paths: List[Path]):
        """
        Initializes the builder with all necessary data and configuration.
        """
        logger.info("Initializing AnnData Builder...")
        self.config = config
        self.settings = settings
        self.amplicons = amplicons
        self.editing_summary = editing_summary
        self.quality_scores = quality_scores
        self.allele_parquet_paths = allele_parquet_paths

        self.cell_barcodes = self.editing_summary.index.tolist()
        self.amplicon_names = extract_amplicon_names(self.editing_summary)
        
        self.n_obs = len(self.cell_barcodes)
        self.n_vars = len(self.amplicon_names)

        # create index maps for efficient numpy array population
        self._cell_map = {name: i for i, name in enumerate(self.cell_barcodes)}
        self._amp_map = {name: i for i, name in enumerate(self.amplicon_names)}
        
        logger.info(f"Builder configured for {self.n_obs} cells and {self.n_vars} amplicons.")

    def _build_obs(self) -> pd.DataFrame:
        """Constructs the .obs (cell metadata) DataFrame."""
        logger.info("Building .obs DataFrame (cell metadata)...")
        obs_df = pd.DataFrame(index=self.cell_barcodes)
        obs_df.index = obs_df.index.astype(str)
        obs_df = obs_df.join(self.quality_scores)
        obs_df.fillna({
            'Amplicon Score': 0.0,
            'Read Count': 0,
            'Barcode Rank': -1,
            'Color': 'Unknown'
        }, inplace=True)
        return obs_df

    def _build_var(self) -> pd.DataFrame:
        """Constructs the .var (amplicon metadata) DataFrame."""
        logger.info("Building .var DataFrame (amplicon metadata)...")
        var_df = self.amplicons.reindex(self.amplicon_names)
        
        count_cols = [f'totCount.{amp}' for amp in self.amplicon_names]
        mod_cols = [f'modPct.{amp}' for amp in self.amplicon_names]
        
        valid_count_cols = [col for col in count_cols if col in self.editing_summary.columns]
        valid_mod_cols = [col for col in mod_cols if col in self.editing_summary.columns]

        if valid_count_cols:
            counts_matrix = self.editing_summary[valid_count_cols]
            counts_matrix.columns = [c.split('.', 1)[1] for c in counts_matrix.columns]
            var_df['total_coverage'] = counts_matrix.sum(axis=0).astype(int)
            var_df['average_coverage'] = counts_matrix.mean(axis=0)

        if valid_mod_cols:
            mod_matrix = self.editing_summary[valid_mod_cols]
            mod_matrix.columns = [c.split('.', 1)[1] for c in mod_matrix.columns]
            var_df['editing_efficiency'] = mod_matrix.mean(axis=0)
            
        var_df.fillna(0, inplace=True)
        return var_df

    def _build_main_layers(self) -> Dict[str, np.ndarray]:
        """Constructs the core data matrices (X and layers) from the editing summary."""
        mod_pct_cols = [f'modPct.{amp}' for amp in self.amplicon_names]
        counts_cols = [f'totCount.{amp}' for amp in self.amplicon_names]
        
        mod_pct_df = self.editing_summary.reindex(columns=mod_pct_cols)
        counts_df = self.editing_summary.reindex(columns=counts_cols)
        
        mod_pct_matrix = mod_pct_df.fillna(0).to_numpy(dtype=np.float32)
        counts_matrix = counts_df.fillna(0).to_numpy(dtype=np.int32)
        
        return {'X': mod_pct_matrix, 'counts': counts_matrix}

    def _build_allele_layers_iteratively(self) -> Dict[str, np.ndarray]:
        """
        Iteratively processes Parquet files to build allele layers, keeping
        memory usage low. This is the new, scalable approach.
        """
        logger.info("Iteratively building allele layers from intermediate Parquet files...")
        
        shape = (self.n_obs, self.n_vars)
        allele_layers = {
            'allele_seq_1': np.full(shape, '', dtype=object),
            'allele_freq_1': np.zeros(shape, dtype=np.float32),
            'allele_seq_2': np.full(shape, '', dtype=object),
            'allele_freq_2': np.zeros(shape, dtype=np.float32),
            'second_allele_freq': np.zeros(shape, dtype=np.float32)
        }

        max_len_seen = 1
        for path in self.allele_parquet_paths:
            df = pd.read_parquet(path)
            if df.empty:
                continue
            
            df['total_reads'] = df.groupby(['cell_barcode', 'amplicon_name'], observed=True)['count'].transform('sum')
            df['frequency'] = (df['count'] / df['total_reads']) * 100
            df.sort_values('frequency', ascending=False, inplace=True)
            
            for (cell_bc, amp_name), group in df.groupby(['cell_barcode', 'amplicon_name'], observed=True):
                row_idx = self._cell_map.get(str(cell_bc))
                col_idx = self._amp_map.get(str(amp_name))

                if row_idx is None or col_idx is None:
                    continue

                top1 = group.iloc[0]
                allele_layers['allele_seq_1'][row_idx, col_idx] = top1['allele_sequence']
                allele_layers['allele_freq_1'][row_idx, col_idx] = top1['frequency']
                max_len_seen = max(max_len_seen, len(top1['allele_sequence']))

                if len(group) > 1:
                    top2 = group.iloc[1]
                    allele_layers['allele_seq_2'][row_idx, col_idx] = top2['allele_sequence']
                    allele_layers['allele_freq_2'][row_idx, col_idx] = top2['frequency']
                    allele_layers['second_allele_freq'][row_idx, col_idx] = top2['frequency']
                    max_len_seen = max(max_len_seen, len(top2['allele_sequence']))
        
        for key in ['allele_seq_1', 'allele_seq_2']:
            allele_layers[key] = np.char.encode(allele_layers[key].astype(str), 'utf-8').astype(f'S{max_len_seen}')

        logger.info("Allele layers built successfully.")
        return allele_layers

    def build(self) -> AnnData:
        """Builds the final AnnData object."""
        logger.info("Assembling final AnnData object...")
        
        obs_df = self._build_obs()
        var_df = self._build_var()
        main_layers = self._build_main_layers()
        allele_layers = self._build_allele_layers_iteratively()

        logger.info("Calculating zygosity using vectorized, configured parameters...")
        mod_pct_matrix = main_layers['X']
        second_allele_freq = allele_layers['second_allele_freq']
        
        zyg_params = self.config['analysis_parameters']['zygosity']
        cond_wt = mod_pct_matrix <= zyg_params['wt_max_mod_pct']
        cond_hom = (mod_pct_matrix >= zyg_params['hom_min_mod_pct']) & (second_allele_freq < zyg_params['compound_het_min_allele2_pct'])
        cond_comp_het = (mod_pct_matrix >= zyg_params['hom_min_mod_pct']) & (second_allele_freq >= zyg_params['compound_het_min_allele2_pct'])
        
        zygosity_matrix = np.select(
            [cond_wt, cond_hom, cond_comp_het], 
            [0, 2, 3], 
            default=1
        ).astype(np.int8)

        adata = AnnData(X=main_layers['X'], obs=obs_df, var=var_df)
        adata.layers['counts'] = main_layers['counts']
        adata.layers['zygosity'] = zygosity_matrix
        
        for key, value in allele_layers.items():
            if key != 'second_allele_freq':
                adata.layers[key] = value

        adata.uns['crisprscope_settings'] = self.settings
        adata.uns['crisprscope_config'] = self.config
        
        encoding_map = {
            0: 'WT/WT', 1: 'WT/Mut', 2: 'Mut/Mut (Homozygous)', 
            3: 'Mut/Mut2 (Compound Het)', -1: 'NoData'
        }
        adata.uns['zygosity_encoding'] = {str(k): v for k, v in encoding_map.items()}
        
        logger.info("AnnData object created successfully with all data layers.")
        return adata