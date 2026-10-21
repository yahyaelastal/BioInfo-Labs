"""
Lab 4 – aDNA Random Sampling & Naive Sequence Assembly
Genome : SARS-CoV-2 Wuhan-Hu-1  (NC_045512.2)
Requires: pip install matplotlib
"""

import random
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# ── 1. PARSE FASTA ────────────────────────────────────────────────────────────
def parse_fasta(path):
    seq = []
    with open(path) as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip().upper())
    return "".join(seq)

genome = parse_fasta("GCF_009858895.2_ASM985889v3_genomic.fna")
print(f"Genome loaded: {len(genome):,} bp")

# ── 2. RANDOM SAMPLING ────────────────────────────────────────────────────────
random.seed(42)
reads = []
for _ in range(2000):
    length = random.randint(100, 120)
    start  = random.randint(0, len(genome) - length)
    reads.append((start, length))

print(f"Reads sampled: {len(reads):,}")

# ── 3. COVERAGE ───────────────────────────────────────────────────────────────
coverage = np.zeros(len(genome), dtype=int)
for start, length in reads:
    coverage[start:start + length] += 1

bases_covered   = int((coverage > 0).sum())
bases_uncovered = len(genome) - bases_covered
print(f"Covered: {bases_covered:,} bp ({bases_covered/len(genome)*100:.1f}%)")
print(f"Uncovered: {bases_uncovered} bp — permanently unrecoverable")
print(f"Avg depth: {coverage.mean():.2f}×")

# ── 4. NAIVE GREEDY ASSEMBLY (first 200 reads) ────────────────────────────────
def find_overlap(a, b, min_ov=30):
    for i in range(max(0, len(a) - len(b)), len(a)):
        if b.startswith(a[i:]) and len(a[i:]) >= min_ov:
            return len(a[i:])
    return 0

pool = [genome[s:s+l] for s, l in reads[:200]]
progress = []
it = 0
print("Running assembler on 200 reads (O(n³) — capped for demo)…")
while True:
    it += 1
    best, bi, bj, bov = 29, -1, -1, 0
    for i in range(len(pool)):
        for j in range(len(pool)):
            if i == j: continue
            ov = find_overlap(pool[i], pool[j])
            if ov > best:
                best, bi, bj, bov = ov, i, j, ov
    if bi == -1:
        break
    merged = pool[bi] + pool[bj][bov:]
    pool.append(merged)
    for idx in sorted([bi, bj], reverse=True):
        pool.pop(idx)
    progress.append((it, len(pool), max(len(s) for s in pool)))

iters   = [p[0] for p in progress]
n_cont  = [p[1] for p in progress]
longest = [p[2] for p in progress]
print(f"Assembly done: {len(pool)} contigs, longest = {max(len(s) for s in pool)} bp")

# ── 5. COMPLEXITY SCALING ─────────────────────────────────────────────────────
measured_time = 22.57   # seconds from actual run on 200 reads
k  = measured_time / 200**3
ns    = [100, 200, 500, 1000, 2000, 5000, 10000]
times = [k * n**3 for n in ns]

def human_time(s):
    if s < 60:   return f"{s:.1f}s"
    if s < 3600: return f"{s/60:.0f}min"
    if s < 86400: return f"{s/3600:.1f}h"
    return f"{s/86400:.0f}d"

# ── 6. PLOTS ──────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 8))
fig.suptitle("Lab 4 – aDNA Random Sampling & Assembly Results", fontsize=13, fontweight="bold")

# Fig 1 – per-base coverage
ax = axes[0, 0]
smooth = np.convolve(coverage, np.ones(300)/300, mode="same")
ax.fill_between(range(len(genome)), smooth, alpha=0.3, color="steelblue")
ax.plot(smooth, color="steelblue", linewidth=0.7, label="Coverage (300-bp smooth)")
ax.axhline(coverage.mean(), color="orange", linewidth=1.2,
           linestyle="--", label=f"Mean {coverage.mean():.1f}×")
ax.fill_between(range(len(genome)), 0, smooth.max(),
                where=(coverage == 0), color="red", alpha=0.7, label="Uncovered")
ax.set_title("Per-base Coverage Across Genome")
ax.set_xlabel("Genomic position (bp)")
ax.set_ylabel("Coverage depth (×)")
ax.legend(fontsize=7)

# Fig 2 – depth histogram
ax = axes[0, 1]
ax.hist(coverage, bins=range(0, int(coverage.max()) + 2),
        color="steelblue", edgecolor="white", linewidth=0.4)
ax.axvline(coverage.mean(), color="orange", linestyle="--",
           linewidth=1.5, label=f"Mean {coverage.mean():.1f}×")
ax.set_title("Coverage Depth Distribution")
ax.set_xlabel("Depth (×)")
ax.set_ylabel("Number of bases")
ax.legend(fontsize=8)

# Fig 3 – assembly progress
ax = axes[1, 0]
ax.plot(iters, n_cont, color="steelblue", linewidth=1.5, label="Contigs remaining")
ax2 = ax.twinx()
ax2.plot(iters, longest, color="darkorange", linewidth=1.5, label="Longest contig (bp)")
ax.set_title("Greedy Assembler Progress (200 reads)")
ax.set_xlabel("Merge iteration")
ax.set_ylabel("Contigs remaining", color="steelblue")
ax2.set_ylabel("Longest contig (bp)", color="darkorange")
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1 + lines2, labels1 + labels2, fontsize=7)

# Fig 4 – O(n³) scaling
ax = axes[1, 1]
bar_colors = ["orange" if n == 200 else "steelblue" for n in ns]
bars = ax.bar([str(n) for n in ns], times, color=bar_colors, edgecolor="white")
for bar, t in zip(bars, times):
    ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() * 1.1,
            human_time(t), ha="center", va="bottom", fontsize=7.5)
ax.set_yscale("log")
ax.set_title("Estimated O(n³) Assembly Time")
ax.set_xlabel("Number of reads")
ax.set_ylabel("Time (log scale)")
ax.text(1.5, times[1] * 3, "← measured", fontsize=7, color="orange")

plt.tight_layout()
Path("figures").mkdir(exist_ok=True)
plt.savefig("figures/lab4_results.png", dpi=150)
print("Saved → figures/lab4_results.png")
plt.show()
