# Lab 3 – Codon Frequency Analyzer

Counts codon frequencies across all 6 reading frames (3 forward + 3 reverse complement) of the SARS-CoV-2 genome, identifies the top 10 most frequent codons, and uses a local AI model to suggest foods low in the associated amino acids.

## Files

| File | Description |
|---|---|
| `codon_frequency.py` | Main script |
| `GCF_009858895.2_ASM985889v3_genomic.fna` | SARS-CoV-2 Wuhan-Hu-1 genome (NC_045512.2, NCBI) |

## Dependencies

```bash
pip install matplotlib transformers torch
```

## How to run

```bash
python codon_frequency.py
```

## What it does

1. **Parses** the FASTA file and extracts the genome sequence
2. **Counts codons** across all 6 reading frames
3. **Ranks** codons from most to least frequent and displays the top 10 as a table and bar chart
4. **Queries TinyLlama** (local AI model, ~1.1 GB, downloaded once) with the top 10 codons and their amino acids, and returns common foods lowest in those amino acids

## AI model notes

- Model: `TinyLlama/TinyLlama-1.1B-Chat-v1.0` (Hugging Face)
- Downloaded automatically on first run to `C:\Users\<you>\.cache\huggingface\hub\`
- Runs entirely offline — no API key needed
- First run takes 3–15 minutes depending on internet speed; subsequent runs load in ~20 seconds
