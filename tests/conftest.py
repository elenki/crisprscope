"""
Pytest fixtures and simulation engine for CRISPRScope testing.
"""

import pytest
import pandas as pd
import numpy as np
import gzip
import yaml
from pathlib import Path
import random
import string

class CRISPRScopeSimulator:
    """
    Generates a complete, temporary directory structure mimicking the output
    of a CRISPRScope orchestrator run, based on a defined ground truth.
    """
    def __init__(self, output_dir: Path):
        self.output_dir = output_dir
        self.crispresso_dir = self.output_dir / "settings.txt.crispresso.filtered"
        self.crispresso_dir.mkdir(parents=True, exist_ok=True)
        self.cell_barcodes = []
        self.config = {}

    def _generate_barcodes(self, n_cells: int):
        """Generates a list of unique, random 18bp cell barcodes."""
        barcodes = set()
        while len(barcodes) < n_cells:
            bc = ''.join(random.choices("ATCG", k=18))
            barcodes.add(bc)
        self.cell_barcodes = sorted(list(barcodes))
        # Write the barcodes to a file for the orchestrator to use
        with open(self.output_dir / "barcodes.txt", "w") as f:
            for bc in self.cell_barcodes:
                f.write(f"{bc}\n")

    def _write_config_file(self):
        """Writes a default config.yaml file for the simulation."""
        self.config = {
            'analysis_parameters': {
                'zygosity': {
                    'wt_max_mod_pct': 20,
                    'het_max_mod_pct': 80,
                    'hom_min_mod_pct': 80,
                    'compound_het_min_allele2_pct': 20
                }
            },
            'read_structure': {
                'mission_bio_v2': {
                    'barcode_1_start': 5,
                    'barcode_1_end': 14,
                    'barcode_2_start': 20,
                    'barcode_2_end': 29
                }
            },
            'performance_parameters': {
                'demultiplexer_chunk_size': 1_000_000
            }
        }
        with open(self.output_dir / "config.yaml", "w") as f:
            yaml.dump(self.config, f)
        
    def _write_settings_file(self, amplicons_path: Path):
        r1_path = self.output_dir / "R1.fastq.gz"
        r2_path = self.output_dir / "R2.fastq.gz"
        barcodes_path = self.output_dir / "barcodes.txt"
        settings_content = (
            f"r1\t{r1_path}\n"
            f"r2\t{r2_path}\n"
            f"barcodes\t{barcodes_path}\n"
            f"amplicons\t{amplicons_path}\n"
        )
        with open(self.output_dir / "settings.txt", "w") as f:
            f.write(settings_content)

    def _write_amplicons_file(self, amplicon_specs: dict) -> Path:
        path = self.output_dir / "amplicons.txt"
        with open(path, "w") as f:
            for name, seq in amplicon_specs.items():
                f.write(f"{name}\t{seq}\tNA\n")
        return path

    def _write_raw_fastqs(self, edit_scenarios: list):
        """Generates raw R1 and R2 fastq files from edit scenarios."""
        r1_path = self.output_dir / "R1.fastq.gz"
        r2_path = self.output_dir / "R2.fastq.gz"
        
        bc_params = self.config['read_structure']['mission_bio_v2']
        read_len = 50 # Arbitrary length for the rest of the read

        with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
            read_counter = 0
            for scenario in edit_scenarios:
                full_barcode = self.cell_barcodes[scenario['cell_id']]
                bc1 = full_barcode[:9]
                bc2 = full_barcode[9:]

                for allele_seq, count in scenario['alleles'].items():
                    for _ in range(count):
                        # Construct R1 read with embedded barcodes
                        r1_seq_list = ['N'] * read_len
                        r1_seq_list[bc_params['barcode_1_start']:bc_params['barcode_1_end']] = list(bc1)
                        r1_seq_list[bc_params['barcode_2_start']:bc_params['barcode_2_end']] = list(bc2)
                        r1_seq = "".join(r1_seq_list)
                        
                        # R2 read is just the allele
                        r2_seq = allele_seq

                        # Write FASTQ records
                        read_id = f"@READ_{read_counter}"
                        qual = "F" * len(r1_seq)
                        f1.write(f"{read_id}\n{r1_seq}\n+\n{qual}\n")
                        
                        qual = "F" * len(r2_seq)
                        f2.write(f"{read_id}\n{r2_seq}\n+\n{qual}\n")
                        
                        read_counter += 1

    def _write_crispresso_fastqs(self, edit_scenarios: list, amplicon_specs: dict):
        """Generates and writes the gzipped FASTQ files for each amplicon."""
        # First, organize scenarios by amplicon for efficient file writing
        scenarios_by_amp = {name: [] for name in amplicon_specs}
        for scenario in edit_scenarios:
            scenarios_by_amp[scenario['amp_name']].append(scenario)

        for amp_name, scenarios in scenarios_by_amp.items():
            amp_dir = self.crispresso_dir / f"CRISPResso_on_{amp_name}"
            amp_dir.mkdir(exist_ok=True)
            fastq_path = amp_dir / "CRISPResso_output.fastq.gz"
            
            with gzip.open(fastq_path, "wt") as f:
                for scenario in scenarios:
                    cell_bc = self.cell_barcodes[scenario['cell_id']]
                    for allele_seq, count in scenario['alleles'].items():
                        for i in range(count):
                            read_id = f"@READ_{i}:{cell_bc}:{amp_name}"
                            qual_line = "F" * len(allele_seq)
                            f.write(f"{read_id}\n{allele_seq}\n+\n{qual_line}\n")

    def _write_summary_files(self, edit_scenarios: list, amplicon_specs: dict):
        """Generates and writes the summary and score files."""
        # Create a ground truth DataFrame
        records = []
        for scenario in edit_scenarios:
            cell_bc = self.cell_barcodes[scenario['cell_id']]
            total_reads = sum(scenario['alleles'].values())
            ref_seq = amplicon_specs[scenario['amp_name']]
            wt_reads = scenario['alleles'].get(ref_seq, 0)
            mod_reads = total_reads - wt_reads
            mod_pct = (mod_reads / total_reads * 100) if total_reads > 0 else 0
            records.append({
                'cell_barcode': cell_bc,
                'amplicon': scenario['amp_name'],
                'totCount': total_reads,
                'modPct': mod_pct
            })
        
        if not records:
            # Handle case with no edits
            df = pd.DataFrame(columns=['cell_barcode']).set_index('cell_barcode')
        else:
            df = pd.DataFrame(records)
            # Pivot to the wide format expected by the loader
            df = df.pivot_table(
                index='cell_barcode', 
                columns='amplicon', 
                values=['totCount', 'modPct']
            ).fillna(0)
            df.columns = [f'{val}.{col}' for val, col in df.columns]

        df.to_csv(self.output_dir / "settings.txt.filteredEditingSummary.txt", sep="\t")

        # Create a dummy quality score file
        score_df = pd.DataFrame({
            'Amplicon Score': [random.uniform(50, 90) for _ in self.cell_barcodes],
            'Read Count': [random.randint(500, 1000) for _ in self.cell_barcodes],
            'Barcode Rank': range(1, len(self.cell_barcodes) + 1),
            'Color': 'High_score_High_reads'
        }, index=self.cell_barcodes)
        score_df.index.name = "Unnamed: 0"
        score_df.to_csv(self.output_dir / "settings.txt.amplicon_score.txt", sep="\t")

    def generate(self, n_cells: int, amplicon_specs: dict, edit_scenarios: list, generate_raw_fastq: bool = False):
        """Main method to orchestrate the generation of all test files."""
        self._generate_barcodes(n_cells)
        self._write_config_file()
        amplicons_path = self._write_amplicons_file(amplicon_specs)
        self._write_settings_file(amplicons_path)
        
        if generate_raw_fastq:
            self._write_raw_fastqs(edit_scenarios)
        
        # Always generate the post-crispresso files for integrator tests
        self._write_crispresso_fastqs(edit_scenarios, amplicon_specs)
        self._write_summary_files(edit_scenarios, amplicon_specs)


@pytest.fixture
def simulation_factory(tmp_path):
    """A pytest fixture that returns a factory for creating simulation directories."""
    def _create_simulation(n_cells: int, amplicon_specs: dict, edit_scenarios: list, generate_raw_fastq: bool = False) -> Path:
        sim_path = tmp_path / "simulation"
        sim_path.mkdir()
        simulator = CRISPRScopeSimulator(output_dir=sim_path)
        simulator.generate(n_cells, amplicon_specs, edit_scenarios, generate_raw_fastq)
        return sim_path
        
    return _create_simulation