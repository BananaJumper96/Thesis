import gzip
import argparse
from collections import defaultdict
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Process contig to bin, contig to ORF, and protein fasta files.')
    parser.add_argument('--cont2bin', required=True, help='Path to the contig to bin mapping file')
    parser.add_argument('--cont2orf', required=True, help='Path to the contig to ORF mapping file')
    parser.add_argument('--protein_fasta', required=True, help='Path to the protein fasta file')
    parser.add_argument('--output_dir', required=True, help='Directory to save the bin protein FASTA files')
    return parser.parse_args()

def main():
    args = parse_args()

    bin2cont_dict = defaultdict(list)
    cont2orf_dict = defaultdict(list)
    fasta_dict = defaultdict(str)

    # Read contig to bin mapping
    with gzip.open(args.cont2bin, 'rb') as f:
        for line in f:
            line = line.decode('utf-8').strip().split('\t')
            bin2cont_dict[line[1]].append(line[0])

    # Read contig to ORF mapping
    with gzip.open(args.cont2orf, 'rb') as f:
        for line in f:
            line = line.decode('utf-8').strip().split('\t')
            contig = line[0]
            orf_id = '_'.join(line[1:])
            cont2orf_dict[contig].append(orf_id)

    # Read protein fasta file
    with gzip.open(args.protein_fasta, 'rb') as f:
        for line in f:
            line = line.decode('utf-8').strip()
            if line.startswith('>'):
                orf_id = line[1:].split(' ')[0]
            else:
                fasta_dict[orf_id] = line

    # Ensure output directory exists
    os.makedirs(args.output_dir, exist_ok=True)

    # Write protein fastas per bin
    for bin_id, contigs in bin2cont_dict.items():
        bin_fasta_file = os.path.join(args.output_dir, f'{bin_id}.faa')
        with open(bin_fasta_file, 'w') as f:
            for contig in contigs:
                for orf_id in cont2orf_dict[contig]:
                    if fasta_dict[orf_id]:
                        f.write(f'>{orf_id}\n{fasta_dict[orf_id]}\n')

if __name__ == '__main__':
    main()