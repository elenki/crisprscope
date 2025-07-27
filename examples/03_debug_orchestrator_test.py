"""
CRISPRScope Orchestrator: Exploratory Demultiplexing Debugger

This script searches for known barcodes within raw FASTQ reads to discover
the correct barcode positions and orientation (forward vs. reverse complement).
"""
import gzip
from pathlib import Path
import sys

# Add the src directory to the Python path
project_root = Path(__file__).resolve().parent.parent
sys.path.append(str(project_root / "src"))

# --- Configuration ---
N_READS_TO_SCAN = 1000  # Scan a larger number of reads to increase chances of finding a match
N_BARCODES_TO_TEST = 100 # Use a subset of known barcodes for the search
BASE_DATA_PATH = Path("/uufs/chpc.utah.edu/common/home/u6046470/clement/projects/20230818_scCRISPR/data/20200804_BaF3_revision")
R1_PATH = BASE_DATA_PATH / "data/BaF3-NSG_S1_L001_R1_001.fastq.gz"
BARCODES_PATH = BASE_DATA_PATH / "v2-barcodes.txt"

def reverse_complement(seq: str) -> str:
    """Computes the reverse complement of a DNA sequence."""
    complement_map = str.maketrans("ATCGN", "TAGCN")
    return seq.translate(complement_map)[::-1]

def main():
    print("\n-------------------------------------------------")
    print("--- CRISPRScope Exploratory Debugger ---")
    print("-------------------------------------------------\n")

    # 1. Load the set of known barcodes
    try:
        with open(BARCODES_PATH, 'r') as f:
            all_known_barcodes = {line.strip() for line in f}
        
        # Create a test set and a reverse-complemented set
        barcodes_to_test = list(all_known_barcodes)[:N_BARCODES_TO_TEST]
        rc_barcodes_to_test = [reverse_complement(bc) for bc in barcodes_to_test]

        print(f"✅ Loaded {len(all_known_barcodes)} known barcodes.")
        print(f"   Searching for the first {N_BARCODES_TO_TEST} barcodes (and their RCs)...")
        print(f"   Sample barcode: {barcodes_to_test[0]} -> RC: {rc_barcodes_to_test[0]}")

    except Exception as e:
        print(f"❌ Failed to load barcodes file: {e}")
        return

    # 2. Scan the first N reads from the R1 FASTQ file
    print(f"\n--- Scanning first {N_READS_TO_SCAN} reads from {R1_PATH.name} ---\n")
    found_match = False
    try:
        with gzip.open(R1_PATH, 'rt') as f:
            for i in range(N_READS_TO_SCAN):
                header = f.readline()
                if not header:
                    break
                sequence = f.readline().strip()
                f.readline()  # Skip '+'
                f.readline()  # Skip quality

                # Search for forward barcodes
                for bc in barcodes_to_test:
                    position = sequence.find(bc)
                    if position != -1:
                        print(f"✅ SUCCESS: Found FORWARD barcode '{bc}' in read {i+1} at position {position}.")
                        found_match = True
                        break # Stop searching this read
                
                if found_match:
                    break # Stop scanning more reads

                # Search for reverse-complemented barcodes
                for rc_bc in rc_barcodes_to_test:
                    position = sequence.find(rc_bc)
                    if position != -1:
                        print(f"✅ SUCCESS: Found REVERSE COMPLEMENT barcode '{rc_bc}' in read {i+1} at position {position}.")
                        found_match = True
                        break # Stop searching this read
                
                if found_match:
                    break # Stop scanning more reads

        if not found_match:
            print(f"❌ FAILED: No matches found in the first {N_READS_TO_SCAN} reads.")
            print("   This suggests a fundamental issue with either the barcode list or the read structure.")

    except Exception as e:
        print(f"❌ An error occurred while reading the FASTQ file: {e}")

if __name__ == "__main__":
    main()
