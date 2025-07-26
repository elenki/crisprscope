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
                 editing_summary: pd.DataFrame, quality_scores: pd.DataFrame,
                 alleles_data: Dict[tuple, Dict[str, int]]):
        """
        Initializes the builder with all necessary data.

        Args:
            settings: Parsed data from settings.txt.
            amplicons: Parsed data from amplicons.txt.
            editing_summary: Parsed data from filteredEditingSummary.txt.
            quality_scores: Parsed data from amplicon_score.txt.
            alleles_data: Parsed allele counts from CRISPResso output.
        """
        logger.info("Initializing AnnData Builder...")
        self.settings = settings
        self.amplicons = amplicons
        self.editing_summary = editing_summary
        self.quality_scores = quality_scores
        self.alleles_data = alleles_data

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
        obs_df.fillna({'Amplicon Score': 0.0, 'Read Count': 0, 'Barcode Rank': -1, 'Color': 'Unknown'}, inplace=True)
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
            count_col, mod_col = f'totCount.{amp}', f'modPct.{amp}'
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
        """Constructs the core data matrices (X and layers) from the editing summary."""
        logger.info("Building core data matrices (X, layers)...")
        mod_pct_cols = [f'modPct.{amp}' for amp in self.amplicon_names]
        counts_cols = [f'totCount.{amp}' for amp in self.amplicon_names]
        mod_pct_df = self.editing_summary.reindex(columns=mod_pct_cols)
        counts_df = self.editing_summary.reindex(columns=counts_cols)
        mod_pct_matrix = mod_pct_df.fillna(0).to_numpy(dtype=np.float32)
        counts_matrix = counts_df.fillna(0).to_numpy(dtype=np.int32)
        logger.info("Finished building core data matrices.")
        return mod_pct_matrix, counts_matrix

    def _build_advanced_layers(self) -> Tuple[np.ndarray, np.ndarray]:
        """Constructs the advanced zygosity and allele layers."""
        logger.info("Building advanced layers (zygosity, alleles)...")
        
        # Zygosity encoding: 0=WT/WT, 1=WT/Mut, 2=Mut/Mut, 3=Mut/Mut2, -1=NoData
        zygosity_matrix = np.full((self.n_obs, self.n_vars), -1, dtype=np.int8)
        alleles_matrix = np.empty((self.n_obs, self.n_vars), dtype=object)

        # Get the reference sequence for each amplicon for comparison
        ref_seqs = self.amplicons['sequence'].to_dict()

        for i, cell_bc in enumerate(self.cell_barcodes):
            for j, amp_name in enumerate(self.amplicon_names):
                allele_summary = self.alleles_data.get((cell_bc, amp_name))
                
                if not allele_summary:
                    alleles_matrix[i, j] = "NoData"
                    continue

                total_reads = sum(allele_summary.values())
                ref_seq = ref_seqs.get(amp_name, "")
                ref_reads = allele_summary.get(ref_seq, 0)
                mod_reads = total_reads - ref_reads
                mod_pct = (mod_reads / total_reads) * 100 if total_reads > 0 else 0

                # --- Zygosity Classification ---
                if mod_pct < 20:
                    zygosity_matrix[i, j] = 0 # WT/WT
                elif 20 <= mod_pct < 80:
                    zygosity_matrix[i, j] = 1 # WT/Mut
                else: # mod_pct >= 80
                    sorted_alleles = sorted(allele_summary.items(), key=lambda item: item[1], reverse=True)
                    if len(sorted_alleles) > 1:
                        second_allele_freq = (sorted_alleles[1][1] / total_reads) * 100
                        if second_allele_freq < 20:
                            zygosity_matrix[i, j] = 2 # Mut/Mut (Homozygous)
                        else:
                            zygosity_matrix[i, j] = 3 # Mut/Mut2 (Compound Heterozygous)
                    else:
                        zygosity_matrix[i, j] = 2 # Mut/Mut (Only one allele type found)
                
                # --- Allele String Formatting ---
                sorted_alleles = sorted(allele_summary.items(), key=lambda item: item[1], reverse=True)[:2]
                allele_str_parts = []
                for allele, count in sorted_alleles:
                    allele_type = "WT" if allele == ref_seq else "ALT"
                    freq = (count / total_reads) * 100
                    allele_str_parts.append(f"{allele_type}:{allele}:{freq:.1f}%")
                alleles_matrix[i, j] = "|".join(allele_str_parts)

        logger.info("Finished building advanced layers.")
        return zygosity_matrix, alleles_matrix

    def build(self) -> AnnData:
        """Builds the final AnnData object."""
        logger.info("Assembling final AnnData object...")
        
        obs_df = self._build_obs()
        var_df = self._build_var()
        mod_pct_matrix, counts_matrix = self._build_data_matrices()
        zygosity_matrix, alleles_matrix = self._build_advanced_layers()
        
        adata = AnnData(X=mod_pct_matrix, obs=obs_df, var=var_df)
        adata.layers['counts'] = counts_matrix
        adata.layers['zygosity'] = zygosity_matrix
        adata.layers['alleles'] = alleles_matrix
        
        adata.uns['crisprscope_settings'] = self.settings
        adata.uns['zygosity_encoding'] = {
            0: 'WT/WT', 1: 'WT/Mut', 2: 'Mut/Mut (Homozygous)', 
            3: 'Mut/Mut2 (Compound Het)', -1: 'NoData'
        }
        
        logger.info("AnnData object created successfully with all data layers.")
        return adata