# Lab 4 – aDNA Random Sampling & Naive Sequence Assembly

## Overview
This lab simulates shotgun sequencing by taking 2,000 random reads from the SARS-CoV-2 Wuhan-Hu-1 genome (NC_045512.2, 29,903 bp) and attempting to reassemble the original sequence using a naive greedy overlap algorithm.

## Files
- `GCF_009858895.2_ASM985889v3_genomic.fna` — reference genome (NCBI)
- `lab4_aDNA_assembly.py` — main script (sampling, assembly, plots)
- `figures/lab4_results.png` — output charts

## How to run
```bash
pip install matplotlib
python lab4_aDNA_assembly.py
```

## Results
| Metric | Value |
|---|---|
| Genome length | 29,903 bp |
| Reads sampled | 2,000 (100–120 bp each) |
| Bases covered | 99.9% |
| Average depth | 7.38× |
| Contigs after assembly | 104 (from 200 reads) |
| Longest contig | 548 bp |

## Main Problem with the Algorithm

The naive greedy overlap assembler has **five fundamental problems**:

### 1. O(n³) Time Complexity — the biggest problem
The algorithm compares every read against every other read to find the best overlap, then repeats after each merge. This makes it **O(n²) per step and O(n³) overall**.

With 2,000 reads that means roughly **8 billion operations**. On 200 reads it already takes ~22 seconds. Scaling to the full 2,000 would take hours, and real datasets have millions of reads — making this approach completely impractical.

### 2. Incomplete Coverage
Even with 2,000 random reads, 38 bp of the genome were never sampled. Those regions are permanently unrecoverable — no algorithm can reconstruct sequence it never saw. Random sampling is inherently uneven.

### 3. Repeat Regions Cause Wrong Merges
SARS-CoV-2 contains repeated sequence motifs (e.g. transcription regulatory sequences). Two reads from *different* locations can look identical, causing the assembler to merge them incorrectly and produce **chimeric contigs** — sequences that never existed in the real genome.

### 4. Reads Have No Known Position
In real sequencing the positional information is lost. The assembler must compare every pair blindly with no way to verify whether an overlap is genuine or a coincidence.

### 5. No Error Correction
Ancient DNA carries deamination damage (C→T and G→A mutations). A naive character-match overlap finder treats damaged bases as mismatches, breaking valid overlaps and further fragmenting the assembly.

## How Real Assemblers Solve This
- **de Bruijn graphs** — break reads into k-mers and build a graph; assembly becomes finding an Eulerian path (polynomial time)
- **Reference-guided mapping** — align reads to a known reference genome (BWA, Bowtie2)
- **aDNA-specific tools** — mapDamage, EAGER handle deamination before assembly
