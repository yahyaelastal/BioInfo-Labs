import matplotlib
matplotlib.use("TkAgg")          # interactive window; change to "Agg" to suppress
import matplotlib.pyplot as plt

# ── Sequence ──────────────────────────────────────────────────────────────────
S = "CGGACTGATCTATCT" + "A" * 27 + "CGTAGCATCTATCGATCTATCTAGCGATCTA"
SW_LEN = 30

print(f"Sequence length : {len(S)}")
print(f"Sliding window  : {SW_LEN}")
print(f"Number of windows: {len(S) - SW_LEN + 1}\n")

# ── CG content (%) ────────────────────────────────────────────────────────────
# Returns percentage of C+G bases in the window
def cg_content(window):
    return (window.count("C") + window.count("G")) / len(window) * 100


# ── Information content σ(S) ──────────────────────────────────────────────────
# σ(S) = 100 - 100 × (1/(|S|-1)) × Σ_{v=1}^{|S|-1} [ matches(v) / (|S|-v) ]
# where matches(v) = number of positions i where S[i] == S[i+v]
# High σ → high complexity; Low σ → repetitive/low information
def information_content(s):
    n = len(s)
    total = 0.0
    for v in range(1, n):
        matches = sum(1 for i in range(n - v) if s[i] == s[i + v])
        total += matches / (n - v)
    return 100.0 - (total / (n - 1)) * 100.0


# ── Center of weight ──────────────────────────────────────────────────────────
# Weighted average of window positions by their signal value
def center_of_weight(values):
    total = sum(values)
    return sum(i * v for i, v in enumerate(values)) / total if total else 0.0


# ── Sliding window scan ───────────────────────────────────────────────────────
positions = list(range(len(S) - SW_LEN + 1))
cg_vals   = [cg_content(S[i : i + SW_LEN])         for i in positions]
ic_vals   = [information_content(S[i : i + SW_LEN]) for i in positions]

avg_cg = sum(cg_vals) / len(cg_vals)
avg_ic = sum(ic_vals) / len(ic_vals)
print(f"Average CG content : {avg_cg:.2f}  (board reference ≈ 29.27)")
print(f"Average IC         : {avg_ic:.2f}  (board reference ≈ 27.53)")

cow_cg = center_of_weight(cg_vals)
cow_ic = center_of_weight(ic_vals)
print(f"Center of weight CG: {cow_cg:.2f}")
print(f"Center of weight IC: {cow_ic:.2f}")

# ── Chart 1 : CG% and IC per sliding window ───────────────────────────────────
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

ax1.plot(positions, cg_vals, color="steelblue", linewidth=1.5)
ax1.axhline(avg_cg, color="steelblue", linestyle="--", linewidth=1,
            label=f"avg = {avg_cg:.2f}%")
ax1.set_ylabel("CG content (%)")
ax1.set_title(f"DNA Pattern Analysis — CG Content (SW={SW_LEN})")
ax1.legend()
ax1.grid(alpha=0.3)

ax2.plot(positions, ic_vals, color="tomato", linewidth=1.5)
ax2.axhline(avg_ic, color="tomato", linestyle="--", linewidth=1,
            label=f"avg = {avg_ic:.2f}")
ax2.set_ylabel("Information Content σ")
ax2.set_xlabel("Window start position")
ax2.set_title(f"DNA Pattern Analysis — Information Content (SW={SW_LEN})")
ax2.legend()
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("pattern_chart.png", dpi=120)
print("\nChart 1 saved → pattern_chart.png")
plt.show()

# ── Chart 2 : Center of weight for each pattern ───────────────────────────────
fig2, ax = plt.subplots(figsize=(6, 4))
labels  = ["CG Center", "IC Center"]
centers = [cow_cg, cow_ic]
colors  = ["steelblue", "tomato"]
bars = ax.bar(labels, centers, color=colors, width=0.4)
for bar, val in zip(bars, centers):
    ax.text(bar.get_x() + bar.get_width() / 2, val + 0.3,
            f"{val:.2f}", ha="center", fontsize=11)
ax.set_ylabel("Window position")
ax.set_title("Center of Weight of Each Pattern")
ax.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig("centers_chart.png", dpi=120)
print("Chart 2 saved → centers_chart.png")
plt.show()
