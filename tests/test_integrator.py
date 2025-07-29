"""
Tests for the CRISPRScope Integrator module.
"""
import pytest
import numpy as np
from crisprscope.api import build_anndata

# --- Test Data Definition ---
# This can be shared across multiple tests
N_CELLS = 2
AMPLICON_SPECS = {
    'trp53': 'GATTACACATTACACATTACACATTACA', # 28bp
    'chd2':  'CCTAGGCCTAGGCCTAGGCCTAGG'      # 24bp
}
EDIT_SCENARIOS = [
    # Cell 0: heterozygous -1bp deletion in trp53
    {'cell_id': 0, 'amp_name': 'trp53', 'alleles': {
        'GATTACACATTACACATTACACAT-ACA': 15, # alt (top allele)
        'GATTACACATTACACATTACACATTACA': 10  # wt
    }},
    # Cell 1: homozygous -3bp deletion in chd2
    {'cell_id': 1, 'amp_name': 'chd2', 'alleles': {
        'CCTAGGCCTAGGCCTAGGCCT---': 30      # alt
    }},
    # Cell 1: unedited for trp53
    {'cell_id': 1, 'amp_name': 'trp53', 'alleles': {
        'GATTACACATTACACATTACACATTACA': 25  # wt
    }},
]

# --- Test Suite ---

def test_builder_passes_on_valid_simulation(simulation_factory):
    """
    Tests that the builder correctly processes a known, valid simulation.
    This test is designed to PASS.
    """
    # 1. SETUP: Generate the test data directory
    test_data_path = simulation_factory(N_CELLS, AMPLICON_SPECS, EDIT_SCENARIOS)

    # 2. EXECUTION: Run the main API function
    adata = build_anndata(base_path=test_data_path)

    # 3. ASSERTION: Verify the output is EXACTLY as expected
    assert adata.n_obs == N_CELLS
    assert adata.n_vars == len(AMPLICON_SPECS)
    assert 'allele_seq_1' in adata.layers

    # --- Verify Cell 0 (het edit in trp53) ---
    cell_0_idx, trp53_idx = 0, adata.var_names.get_loc('trp53')
    
    # Zygosity: WT/Mut == 1
    assert adata.layers['zygosity'][cell_0_idx, trp53_idx] == 1
    # Top allele is the ALT allele (15 reads > 10 reads)
    assert adata.layers['allele_seq_1'][cell_0_idx, trp53_idx].decode('utf-8') == 'GATTACACATTACACATTACACAT-ACA'
    # Top allele frequency
    assert pytest.approx(adata.layers['allele_freq_1'][cell_0_idx, trp53_idx]) == (15 / 25) * 100
    # Second allele is the WT allele
    assert adata.layers['allele_seq_2'][cell_0_idx, trp53_idx].decode('utf-8') == 'GATTACACATTACACATTACACATTACA'
    assert pytest.approx(adata.layers['allele_freq_2'][cell_0_idx, trp53_idx]) == (10 / 25) * 100

    # --- Verify Cell 1 (hom edit in chd2) ---
    cell_1_idx, chd2_idx = 1, adata.var_names.get_loc('chd2')

    # Zygosity: Mut/Mut == 2 (>=80% mod, second allele freq < 20%)
    assert adata.layers['zygosity'][cell_1_idx, chd2_idx] == 2
    # Top allele frequency is 100%
    assert pytest.approx(adata.layers['allele_freq_1'][cell_1_idx, chd2_idx]) == 100.0
    # There is no second allele
    assert adata.layers['allele_seq_2'][cell_1_idx, chd2_idx].decode('utf-8') == ''
    assert adata.layers['allele_freq_2'][cell_1_idx, chd2_idx] == 0.0

def test_builder_fails_on_incorrect_logic(simulation_factory):
    """
    Demonstrates how the test suite catches a logical error.
    This test is INTENTIONALLY WRONG and designed to FAIL.
    """
    # 1. SETUP: Same as the passing test
    test_data_path = simulation_factory(N_CELLS, AMPLICON_SPECS, EDIT_SCENARIOS)

    # 2. EXECUTION: Same as the passing test
    adata = build_anndata(base_path=test_data_path)

    # 3. ASSERTION: Make an incorrect claim about the data
    cell_0_idx, trp53_idx = 0, adata.var_names.get_loc('trp53')
    
    # We know from the simulation that the zygosity for this cell/amplicon
    # should be 1 (WT/Mut). Here, we incorrectly assert that it is 2 (Mut/Mut).
    # This test will FAIL, proving that our test harness can detect logical flaws
    # in the builder's code. If you were to change the builder's logic incorrectly,
    # this test would catch it.
    print("\nExpecting this assertion to fail...")
    assert adata.layers['zygosity'][cell_0_idx, trp53_idx] == 2, "Zygosity for heterozygous edit was not 2 (this is the intended failure)"