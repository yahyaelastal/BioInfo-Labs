import math
import matplotlib.pyplot as plt

# ── Parameters ─────────────────────────────────────────────────────────────────
SEQUENCE = "TGACGTATA"
NA_CONC = 0.050  # 50 mM in Mol/L (M)


def tm_simple(seq):
    """Tm = 4*(G+C) + 2*(A+T)"""
    gc = seq.count("G") + seq.count("C")
    at = seq.count("A") + seq.count("T")
    return 4 * gc + 2 * at


def tm_salt(seq, na):
    """Tm = 81.5 + 16.6*log10([Na+]) + 0.41*(%GC) - 600/length"""
    n = len(seq)
    gc = seq.count("G") + seq.count("C")
    pct_gc = (gc / n) * 100
    return 81.5 + 16.6 * math.log10(na) + 0.41 * pct_gc - (600 / n)


# ── Calculations ───────────────────────────────────────────────────────────────
seq = SEQUENCE.upper().strip()
n = len(seq)
gc_count = seq.count("G") + seq.count("C")
pct_gc = (gc_count / n) * 100

t1 = tm_simple(seq)
t2 = tm_salt(seq, NA_CONC)

# ── Plotting ───────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(13, 5))
gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.4])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[:, 1])

# 1. Base Frequency
bases = ["A", "T", "G", "C"]
counts = [seq.count(b) for b in bases]
colors = ["#2ecc71", "#e74c3c", "#e67e22", "#3498db"]
bars1 = ax1.bar(bases, counts, color=colors, edgecolor="black", linewidth=0.5)
ax1.set_title(f"Base Frequency ({n} bp)", fontsize=10, fontweight="bold")
ax1.set_ylabel("Count")
ax1.set_ylim(0, max(counts) + 1)
for bar, c in zip(bars1, counts):
    ax1.text(bar.get_x() + bar.get_width() / 2, c + 0.05, str(c), ha="center", va="bottom", fontsize=9)

# 2. Tm Formula Comparison
labels = ["Simple", "Salt-adj."]
values = [t1, t2]
bars2 = ax2.bar(labels, values, color=["#9b59b6", "#1abc9c"], edgecolor="black", linewidth=0.5)
ax2.set_title("Melting Temperature Comparison", fontsize=10, fontweight="bold")
ax2.set_ylabel("Tm (°C)")
margin = max(6, abs(t1 - t2) * 2 + 2)
ax2.set_ylim(min(values) - margin, max(values) + margin)
for bar, v in zip(bars2, values):
    ax2.text(bar.get_x() + bar.get_width() / 2, v + margin * 0.05, f"{v:.1f} °C", ha="center", fontsize=9, fontweight="bold")

# 3. Tm vs %GC Content Curve
gc_pcts = list(range(0, 101))
curve_simple = [4 * (n * p / 100) + 2 * (n - (n * p / 100)) for p in gc_pcts]
curve_salt = [81.5 + 16.6 * math.log10(NA_CONC) + 0.41 * p - (600 / n) for p in gc_pcts]

ax3.plot(gc_pcts, curve_simple, color="#9b59b6", linewidth=2, label="Simple: 4(G+C) + 2(A+T)")
ax3.plot(gc_pcts, curve_salt, color="#1abc9c", linewidth=2, label=f"Salt-adj.: [Na⁺]={NA_CONC*1000:.0f} mM")

# Current Sequence Markers
ax3.scatter([pct_gc], [t1], color="#9b59b6", s=70, zorder=5)
ax3.scatter([pct_gc], [t2], color="#1abc9c", s=70, zorder=5)
ax3.axvline(pct_gc, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax3.annotate(
    f"Current Seq\n({pct_gc:.1f}% GC)",
    xy=(pct_gc, (t1 + t2) / 2),
    xytext=(pct_gc + 6, (t1 + t2) / 2),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="gray", lw=0.8),
)

ax3.set_title(f"Tm vs. %GC (Length = {n} bp)", fontsize=11, fontweight="bold")
ax3.set_xlabel("% GC Content", fontsize=10)
ax3.set_ylabel("Tm (°C)", fontsize=10)
ax3.grid(True, linestyle="--", alpha=0.5)
ax3.legend(fontsize=9)

plt.tight_layout()
plt.show()