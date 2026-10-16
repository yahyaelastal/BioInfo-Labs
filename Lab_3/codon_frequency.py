"""
Codon Frequency Analyzer
Reads a FASTA file (virus or bacteria genome), counts codon frequencies
using a sliding window of 3, and reports the top 10 codons with their
associated amino acids.
"""

import os
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# Standard genetic code (DNA codon table)
CODON_TABLE = {
    'TTT': 'Phe (F)', 'TTC': 'Phe (F)',
    'TTA': 'Leu (L)', 'TTG': 'Leu (L)',
    'CTT': 'Leu (L)', 'CTC': 'Leu (L)', 'CTA': 'Leu (L)', 'CTG': 'Leu (L)',
    'ATT': 'Ile (I)', 'ATC': 'Ile (I)', 'ATA': 'Ile (I)',
    'ATG': 'Met (M)',
    'GTT': 'Val (V)', 'GTC': 'Val (V)', 'GTA': 'Val (V)', 'GTG': 'Val (V)',
    'TCT': 'Ser (S)', 'TCC': 'Ser (S)', 'TCA': 'Ser (S)', 'TCG': 'Ser (S)',
    'CCT': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
    'ACT': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
    'GCT': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
    'TAT': 'Tyr (Y)', 'TAC': 'Tyr (Y)',
    'TAA': 'Stop (*)', 'TAG': 'Stop (*)',
    'CAT': 'His (H)', 'CAC': 'His (H)',
    'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
    'AAT': 'Asn (N)', 'AAC': 'Asn (N)',
    'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
    'GAT': 'Asp (D)', 'GAC': 'Asp (D)',
    'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
    'TGT': 'Cys (C)', 'TGC': 'Cys (C)',
    'TGA': 'Stop (*)',
    'TGG': 'Trp (W)',
    'CGT': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)',
    'AGT': 'Ser (S)', 'AGC': 'Ser (S)',
    'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
    'GGT': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)',
}


def read_fasta(filepath):
    """Parse a FASTA file and return a dict of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(filepath, 'r') as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header is not None:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]  # strip the '>'
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_header is not None:
        sequences[current_header] = ''.join(current_seq)

    return sequences


def sliding_window_codons(sequence):
    """
    Count codons using a sliding window of size 3.
    Window slides one nucleotide at a time: positions 0-2, 1-3, 2-4, ...
    Only complete 3-nt windows are counted; ambiguous bases (N, etc.) are skipped.
    """
    counts = Counter()
    valid_bases = set('ACGT')

    for i in range(len(sequence) - 2):
        codon = sequence[i:i + 3]
        if all(b in valid_bases for b in codon):
            counts[codon] += 1

    return counts


def get_amino_acid(codon):
    return CODON_TABLE.get(codon, 'Unknown')


def display_results(header, counts, top_n=10):
    total = sum(counts.values())
    top = counts.most_common(top_n)

    print()
    print('=' * 65)
    print(f' Genome: {header[:60]}')
    print(f' Total sliding-window codons counted: {total:,}')
    print('=' * 65)
    print(f'  {"Rank":<5} {"Codon":<8} {"Count":>10} {"Frequency":>12}  {"Amino Acid"}')
    print('-' * 65)

    for rank, (codon, count) in enumerate(top, 1):
        freq = count / total * 100
        aa = get_amino_acid(codon)
        print(f'  {rank:<5} {codon:<8} {count:>10,} {freq:>11.4f}%  {aa}')

    print('=' * 65)
    print()


def plot_results(header, counts, top_n=10):
    total = sum(counts.values())
    top = counts.most_common(top_n)

    codons = [c for c, _ in top]
    freqs  = [cnt / total * 100 for _, cnt in top]
    labels = [f"{c}\n{get_amino_acid(c)}" for c in codons]

    # Color bars by amino acid so synonymous codons share a color
    aa_list = [get_amino_acid(c) for c in codons]
    unique_aa = list(dict.fromkeys(aa_list))
    cmap = plt.colormaps['tab10'].resampled(len(unique_aa))
    color_map = {aa: cmap(i) for i, aa in enumerate(unique_aa)}
    colors = [color_map[aa] for aa in aa_list]

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(range(top_n), freqs, color=colors, edgecolor='white', linewidth=0.8)

    # Value labels on top of bars
    for bar, freq in zip(bars, freqs):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.03,
                f'{freq:.2f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')

    ax.set_xticks(range(top_n))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel('Frequency (%)', fontsize=12)
    ax.set_xlabel('Codon  (Amino Acid)', fontsize=12)
    ax.set_title(f'Top {top_n} Codon Frequencies\n{header[:80]}', fontsize=13, fontweight='bold')
    ax.set_ylim(0, max(freqs) * 1.18)
    ax.yaxis.grid(True, linestyle='--', alpha=0.5)
    ax.set_axisbelow(True)

    # Legend for amino acids
    patches = [mpatches.Patch(color=color_map[aa], label=aa) for aa in unique_aa]
    ax.legend(handles=patches, title='Amino Acid', bbox_to_anchor=(1.01, 1),
              loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.show()


def main():
    filepath = 'GCF_009858895.2_ASM985889v3_genomic.fna'

    if not os.path.isfile(filepath):
        print(f'Error: file not found: {filepath}')
        return

    print(f'\nReading FASTA file: {filepath}')
    sequences = read_fasta(filepath)

    if not sequences:
        print('Error: no sequences found in the FASTA file.')
        return

    print(f'Found {len(sequences)} sequence(s).\n')

    for header, seq in sequences.items():
        print(f'Processing sequence ({len(seq):,} bp) ...')
        counts = sliding_window_codons(seq)
        display_results(header, counts, top_n=10)
        plot_results(header, counts, top_n=10)


if __name__ == '__main__':
    main()
