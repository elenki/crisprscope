"""
CRISPRScope: Utility functions for the Integrator module.
"""
import pandas as pd
from typing import List
import logging

logger = logging.getLogger(__name__)

def extract_amplicon_names(editing_df: pd.DataFrame) -> List[str]:
    """
    Extracts unique amplicon names from the editing summary DataFrame columns.

    It identifies columns matching the 'totCount.<amplicon_name>' or
    'modPct.<amplicon_name>' pattern and returns a sorted, unique list
    of the amplicon names.

    Args:
        editing_df: The DataFrame loaded from the editing summary file.

    Returns:
        A sorted list of unique amplicon names.
    """
    amplicons = set()
    for col in editing_df.columns:
        if col.startswith(('totCount.', 'modPct.')):
            # Splits 'totCount.amplicon_name' into ['totCount', 'amplicon_name']
            # and takes the second part.
            amp_name = col.split('.', 1)[1]
            amplicons.add(amp_name)
    
    amplicon_list = sorted(list(amplicons))
    logger.info(f"Extracted {len(amplicon_list)} unique amplicon names.")
    return amplicon_list