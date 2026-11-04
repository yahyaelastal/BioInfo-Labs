# Lab 6 — Transposable Element Simulator & Detector

## Overview

This lab has two parts:

- **Part 1 (`lab6.py`):** Generates an artificial DNA sequence, simulates the insertion of transposable elements covering all three structural cases, runs a detection algorithm, and visualizes the results.
- **Part 2 (`lab6_part2.py`):** Reads a real genome from a `.fna` file and scans it for transposon candidates using inverted repeat signatures, with deduplication and genome-wide visualization.

---

## What is a Transposable Element?

A transposable element (TE) is a DNA sequence that can move (transpose) from one location in the genome to another. Its structure is:

```
5' ── [Inverted Repeat] ── [Body] ── [RC of Inverted Repeat] ── 3'
```

- **Inverted Repeat (IR):** a short sequence at the left boundary.
- **Body:** the internal sequence of the transposon.
- **Reverse Complement IR (RC-IR):** the reverse complement of the IR, placed at the right boundary.

The IR and RC-IR pair act as the "signature" that the detection algorithm looks for.

---

## Program Structure

### 1. DNA Generation
A random DNA sequence of ~300 bases is generated using the four nucleotides (A, T, G, C).

### 2. Transposon Insertion
Six transposons are inserted to cover all three structural cases:

| Transposon | Case | Description |
|------------|------|-------------|
| TE1 | 1 | Independent |
| TE2 | 1 | Independent |
| TE3 | 2 | Outer — nested case |
| TE4 | 2 | Inner — nested inside TE3 |
| TE5 | 3 | Overlapping with TE6 |
| TE6 | 3 | Overlapping with TE5 |

### 3. Detection Algorithm
The detector searches the DNA for each known IR sequence and locates the matching RC-IR downstream. Each (IR → RC-IR) pair defines the start and end of a transposon.

### 4. Case Classification
After detection, each found transposon is classified:
- **Independent (Case 1):** no spatial relationship with others.
- **Nested (Case 2):** one transposon's range falls completely inside another's.
- **Overlapping (Case 3):** two transposons share a partial boundary region.

---

## The Three Cases Explained

### Case 1 — Independent Transposons
TE1 and TE2 are inserted into separate, non-overlapping regions of the DNA.
Detection is straightforward: find IR, find RC-IR, done.

```
DNA: ───[TE1]─────────────[TE2]────────────────
```

### Case 2 — Nested Transposons
TE4 is inserted **inside the body** of TE3, after TE3 is already in the DNA.

```
DNA: ───[TE3_IR ── [TE4] ── TE3_RC-IR]──────────
```

The detector correctly identifies both:
- TE3 by its outer IR pair (wider span)
- TE4 by its inner IR pair (narrower span, inside TE3)

Both are classified as **nested** because one range is fully contained within the other.

### Case 3 — Overlapping Transposons (and what actually happens)
TE5 is inserted first. TE6 is then inserted at a position **8 bases before the end of TE5**, intending to create an overlap.

**What was intended:**
```
DNA: ───[TE5_start ───────────[TE6_start─── TE5_end]──── TE6_end]───
```

**What actually happens in the sequence:**
When TE6 is inserted at position 406 (inside TE5's right boundary), it physically pushes TE5's closing RC-IR further downstream in the sequence:

```
Before TE6 insertion:
  TE5: positions 380–413  (RC-IR at 408–413)

After TE6 insertion at 406:
  TE6 body (32 bp) is inserted at 406
  TE5's RC-IR is pushed from 408 → 440
  TE5 appears to span: 380–445
  TE6 spans: 406–437  ← entirely inside 380–445
```

**Result:** the detector finds TE5 spanning 380–445 (66 bp instead of the original 34 bp), and TE6 (406–437) falls completely inside that range. The case classifier therefore labels both as **nested**, not overlapping.

---

## Key Finding — Overlapping Transposons Appear Nested

This is the most important observation of this lab:

> When one transposon inserts into the boundary region of another, it displaces the closing inverted repeat of the first transposon further downstream. In the resulting sequence, the displaced boundaries make the two transposons appear nested rather than overlapping.

This is exactly the kind of edge case the professor highlighted: *"intersect 2 or 3 transposons to see how the algorithm behaves in those situations."*

The algorithm is not wrong — it correctly finds where the IR and RC-IR are in the actual sequence. The behavior reveals a real challenge in transposon detection: **overlapping insertions distort each other's boundaries**, making automated classification unreliable without additional context.

---

## Output

### Terminal
- Ground truth table: position, length, IR sequences, and case for each inserted TE.
- Detected table: same attributes as found by the algorithm, plus case classification.

### Graph (`transposons.png`)
Two panels:
- **Top:** Ground truth — colored by case (green = independent, blue = nested, orange = overlapping).
- **Bottom:** Detected — what the algorithm found, with the same color coding.

Dark blocks at each end of a bar represent the inverted repeat regions.

---

## How to Run

**Part 1 — Synthetic DNA simulator:**
```bash
python lab6.py
```

**Part 2 — Real genome scanner:**
```bash
python lab6_part2.py
```

**Requirements:** Python 3, matplotlib

```bash
pip install matplotlib
```

---

## Part 2 — Real Genome Scanner

### File
`lab6_part2.py` — reads any `.fna` (FASTA nucleotide) file and detects transposon candidates.

### Genome Used
`GCF_009858895.2_ASM985889v3_genomic.fna` — SARS-CoV-2 (Wuhan-Hu-1), complete genome, 29,903 bp.

### Detection Parameters

| Parameter | Value | Reason |
|-----------|-------|--------|
| IR length | 8 bp | Longer than 6 bp to reduce random matches in a 30kb genome |
| Min body length | 10 bp | As specified — minimum meaningful transposon body |
| Max TE length | 500 bp | Caps spurious long-range coincidental matches |
| Max IR frequency | 3 occurrences | IRs appearing many times are likely random repeats, not TE signatures |

### Algorithm Steps
1. Build an index of all 8-mers and their positions in the genome.
2. For each 8-mer that appears ≤3 times, check if its reverse complement also appears ≤3 times downstream.
3. If the gap (body) is ≥10 bp and total length ≤500 bp, record it as a candidate.
4. **Deduplicate:** cluster candidates whose start and end positions differ by ≤8 bp (sliding-window artifacts from the same real site) and keep only the shortest-body representative per cluster.
5. Classify each remaining candidate as independent, nested, or overlapping based on spatial relationships.

### Output
- **Terminal:** full list of candidates with position, length, IR sequences, and case.
- **`transposons_genome.png`:** two-panel figure:
  - *Top:* genome map with candidates as colored bars in stacked lanes.
  - *Bottom:* density histogram (candidates per 500 bp window) showing hotspots.

### Results on SARS-CoV-2
- **216 candidates** detected after deduplication.
- Density hotspots observed around positions ~6,500–7,500, ~12,500–13,500, and ~22,000–24,000.

### Key Biological Finding

> **SARS-CoV-2 is an RNA virus and does not have transposable elements.**

All 216 candidates are coincidental inverted repeat matches in the nucleotide sequence — not real transposons. This result demonstrates a critical principle in bioinformatics:

**Algorithmic output must always be interpreted in biological context.** The detector works correctly — it finds every IR pair meeting the structural criteria — but the genome chosen has no actual transposons to find. Running the same algorithm on a eukaryotic genome (plant, insect, or mammal) would yield biologically meaningful candidates, since those genomes contain large families of known TEs (LINE, SINE, DNA transposons, etc.).

The density plot is still informative: regions with many IR-like patterns in SARS-CoV-2 may correspond to RNA secondary structure elements or replication signals, which is a separate but interesting biological observation.

---

## Author
Yahya El-Astal — Bioinformatics Lab 6
