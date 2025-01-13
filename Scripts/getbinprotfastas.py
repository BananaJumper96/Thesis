import gzip
import argparse
from collections import defaultdict
from pathlib import Path
import logging
import sys
from typing import Dict, List, Optional

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Process contig to bin, contig to ORF, and protein FASTA files.'
    )
    parser.add_argument('--cont2bin', required=True, type=Path, help='Path to the contig to bin mapping file')
    parser.add_argument('--cont2orf', required=True, type=Path, help='Path to the contig to ORF mapping file')
    parser.add_argument('--protein_fasta', required=True, type=Path, help='Path to the protein FASTA file')
    parser.add_argument('--output_dir', required=True, type=Path, help='Directory to save the bin protein FASTA files')
    parser.add_argument('--bin_stats', type=Path, help='Path to the bin stats file')
    parser.add_argument('--max_cont', type=float, default=5.0, help='Maximum contamination allowed')
    parser.add_argument('--min_comp', type=float, default=90.0, help='Minimum completion required')
    return parser.parse_args()

def setup_logging() -> logging.Logger:
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    if not logger.handlers:
        logger.addHandler(handler)
    return logger

def read_contig_to_bin(cont2bin_path: Path, logger: logging.Logger) -> Dict[str, List[str]]:
    bin2cont_dict = defaultdict(list)
    logger.info(f"Reading contig to bin mapping from '{cont2bin_path}'")
    try:
        with gzip.open(cont2bin_path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    logger.warning(f"Skipping malformed line {line_num} in cont2bin: {parts}")
                    continue
                contig, bin_id = parts[:2]
                bin2cont_dict[bin_id].append(contig)
        logger.info(f"Loaded {len(bin2cont_dict)} bins from contig to bin mapping.")
    except Exception as e:
        logger.error(f"Failed to read contig to bin mapping: {e}")
        sys.exit(1)
    return bin2cont_dict

def read_contig_to_orf(cont2orf_path: Path, logger: logging.Logger) -> Dict[str, List[str]]:
    cont2orf_dict = defaultdict(list)
    logger.info(f"Reading contig to ORF mapping from '{cont2orf_path}'")
    try:
        with gzip.open(cont2orf_path, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    logger.warning(f"Skipping malformed line {line_num} in cont2orf: {parts}")
                    continue
                contig, orf_id = parts[0], '_'.join(parts[1:])
                cont2orf_dict[contig].append(orf_id)
        logger.info(f"Loaded ORFs for {len(cont2orf_dict)} contigs from contig to ORF mapping.")
    except Exception as e:
        logger.error(f"Failed to read contig to ORF mapping: {e}")
        sys.exit(1)
    return cont2orf_dict

def read_protein_fasta(protein_fasta_path: Path, logger: logging.Logger) -> Dict[str, str]:
    fasta_dict = {}
    logger.info(f"Reading protein FASTA from '{protein_fasta_path}'")
    try:
        with gzip.open(protein_fasta_path, 'rt') as f:
            current_orf = None
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line.startswith('>'):
                    current_orf = line[1:].split(' ')[0]
                    logger.debug(f"Processing ORF ID: {current_orf}")
                else:
                    if current_orf:
                        fasta_dict[current_orf] = line
        logger.info(f"Loaded protein sequences for {len(fasta_dict)} ORFs from protein FASTA.")
    except Exception as e:
        logger.error(f"Failed to read protein FASTA file: {e}")
        sys.exit(1)
    return fasta_dict

def read_bin_stats(bin_stats_path: Path, logger: logging.Logger) -> Dict[str, Dict[str, float]]:
    stats_dict = defaultdict(dict)
    logger.info(f"Reading bin stats from '{bin_stats_path}'")
    try:
        with gzip.open(bin_stats_path, 'rt') as f:
            header = next(f, None)  # Skip header
            if header is None:
                logger.warning("Bin stats file is empty.")
                return stats_dict
            for line_num, line in enumerate(f, 2):
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    logger.warning(f"Skipping malformed line {line_num} in bin_stats: {parts}")
                    continue
                bin_id, comp_str, cont_str = parts[:3]
                try:
                    completion = float(comp_str)
                    contamination = float(cont_str)
                    stats_dict[bin_id]['completion'] = completion
                    stats_dict[bin_id]['contamination'] = contamination
                except ValueError:
                    logger.warning(
                        f"Invalid numerical values in bin_stats for bin '{bin_id}' on line {line_num}: {parts[1:3]}"
                    )
        logger.info(f"Loaded stats for {len(stats_dict)} bins from bin stats file.")
    except Exception as e:
        logger.error(f"Failed to read bin stats file: {e}")
        sys.exit(1)
    return stats_dict

def write_bin_fasta(
    bin_id: str,
    contigs: List[str],
    cont2orf_dict: Dict[str, List[str]],
    fasta_dict: Dict[str, str],
    output_dir: Path,
    logger: logging.Logger
) -> bool:
    bin_fasta_path = output_dir / f"{bin_id}.faa"
    try:
        with bin_fasta_path.open('w') as f:
            for contig in contigs:
                orfs = cont2orf_dict.get(contig, [])
                for orf_id in orfs:
                    sequence = fasta_dict.get(orf_id)
                    if sequence:
                        f.write(f">{orf_id}\n{sequence}\n")
                    else:
                        logger.warning(
                            f"No sequence found for ORF '{orf_id}' in bin '{bin_id}'. "
                            "This is likely not a problem, since the contig to 'ORF' mapping contains other structural element gene IDs."
                        )
        logger.info(f"Written FASTA for bin '{bin_id}' with {len(contigs)} contigs.")
        return True
    except Exception as e:
        logger.error(f"Failed to write FASTA for bin '{bin_id}': {e}")
        return False

def main():
    args = parse_args()
    logger = setup_logging()
    logger.info(f"Script started with arguments: {args}")

    bin2cont_dict = read_contig_to_bin(args.cont2bin, logger)
    cont2orf_dict = read_contig_to_orf(args.cont2orf, logger)
    fasta_dict = read_protein_fasta(args.protein_fasta, logger)

    stats_dict: Dict[str, Dict[str, float]] = {}
    if args.bin_stats:
        if args.bin_stats.is_file():
            stats_dict = read_bin_stats(args.bin_stats, logger)
        else:
            logger.error(f"Bin stats file '{args.bin_stats}' does not exist.")
            sys.exit(1)
    else:
        logger.info("No bin stats file provided. Proceeding without bin filtering.")

    # Ensure output directory exists
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Output directory '{args.output_dir}' is ready.")
    except Exception as e:
        logger.error(f"Failed to create output directory '{args.output_dir}': {e}")
        sys.exit(1)

    # Write protein FASTA files per bin
    logger.info(f"Writing protein FASTA files per bin to '{args.output_dir}'")
    bins_written = 0
    for bin_id, contigs in bin2cont_dict.items():
        if args.bin_stats:
            bin_stats = stats_dict.get(bin_id, {})
            completion = bin_stats.get('completion', 0.0)
            contamination = bin_stats.get('contamination', float('inf'))
            if completion < args.min_comp:
                logger.debug(f"Skipping bin '{bin_id}' due to low completion: {completion:.2f}%")
                continue
            if contamination > args.max_cont:
                logger.debug(f"Skipping bin '{bin_id}' due to high contamination: {contamination:.2f}%")
                continue

        success = write_bin_fasta(
            bin_id, contigs, cont2orf_dict, fasta_dict, args.output_dir, logger
        )
        if success:
            bins_written += 1

    logger.info(f"Finished writing FASTA files for {bins_written} bins.")
    logger.info("Script completed successfully.")

if __name__ == '__main__':
    main()
