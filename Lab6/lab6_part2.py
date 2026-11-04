import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict

# ── Parameters ─────────────────────────────────────────────────────────────────
FNA_FILE   = "GCF_009858895.2_ASM985889v3_genomic.fna"
IR_LEN     = 8      # longer IR = fewer false positives on a real genome
MIN_BODY   = 10     # minimum body length between IR and RC-IR
MAX_TE_LEN = 500    # ignore candidates longer than this (likely spurious)
MAX_IR_OCC = 3      # only use IRs that appear at most this many times in genome

# ── Read FASTA ─────────────────────────────────────────────────────────────────
def read_fasta(path):
    header, seq = "", []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                header = line[1:]
            else:
                seq.append(line.upper())
    return header, "".join(seq)

# ── Helpers ────────────────────────────────────────────────────────────────────
def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return "".join(comp.get(b, 'N') for b in reversed(seq))

def is_palindrome(seq):
    return seq == reverse_complement(seq)

# ── Index all k-mer positions ──────────────────────────────────────────────────
def build_kmer_index(dna, k):
    index = defaultdict(list)
    for i in range(len(dna) - k + 1):
        index[dna[i:i+k]].append(i)
    return index

# ── Scan for transposon candidates ────────────────────────────────────────────
def find_transposons(dna, ir_len, min_body, max_len, max_ir_occ):
    print(f"Building {ir_len}-mer index over {len(dna)} bp genome...")
    index = build_kmer_index(dna, ir_len)
    print(f"Index built. Scanning for IR pairs...")

    candidates = []
    seen = set()

    for ir, positions in index.items():
        # skip palindromes (IR == RC-IR, would match itself)
        if is_palindrome(ir):
            continue
        # skip IRs that appear too many times (random noise)
        if len(positions) > max_ir_occ:
            continue

        rc_ir = reverse_complement(ir)
        rc_positions = index.get(rc_ir, [])
        if not rc_positions or len(rc_positions) > max_ir_occ:
            continue

        for start in positions:
            body_start = start + ir_len
            for end_of_rcir in rc_positions:
                # RC-IR must come AFTER the IR
                if end_of_rcir <= body_start + min_body:
                    continue
                te_end = end_of_rcir + ir_len - 1
                te_len = te_end - start + 1
                body_len = end_of_rcir - body_start
                if body_len < min_body or te_len > max_len:
                    continue

                key = (start, te_end)
                if key in seen:
                    continue
                seen.add(key)

                candidates.append({
                    "start":    start,
                    "end":      te_end,
                    "length":   te_len,
                    "body_len": body_len,
                    "ir":       ir,
                    "rc_ir":    rc_ir,
                })

    candidates.sort(key=lambda x: x["start"])
    return candidates

# ── Classify spatial relationships ────────────────────────────────────────────
def classify(candidates):
    for d in candidates:
        s, e = d["start"], d["end"]
        case = 1
        for other in candidates:
            if other is d:
                continue
            os, oe = other["start"], other["end"]
            if os < s and oe > e:
                case = 2   # d is nested inside other
                break
            if s < os and e > oe:
                case = 2   # d contains other
                break
            if s < os < e < oe:
                case = 3
                break
            if os < s < oe < e:
                case = 3
                break
        d["case"] = case
    return candidates

# ── Deduplicate: merge near-identical candidates (sliding-window artifacts) ────
def deduplicate(candidates, ir_len):
    """
    Cluster candidates whose start AND end positions differ by <= ir_len.
    From each cluster keep only the one with the shortest body (most specific match).
    """
    if not candidates:
        return []
    clusters = []
    used = [False] * len(candidates)
    for i, a in enumerate(candidates):
        if used[i]:
            continue
        cluster = [a]
        used[i] = True
        for j, b in enumerate(candidates):
            if used[j]:
                continue
            if abs(a["start"] - b["start"]) <= ir_len and abs(a["end"] - b["end"]) <= ir_len:
                cluster.append(b)
                used[j] = True
        # keep the candidate with the shortest body in this cluster
        best = min(cluster, key=lambda x: x["body_len"])
        clusters.append(best)
    clusters.sort(key=lambda x: x["start"])
    return clusters

# ── Visualization ──────────────────────────────────────────────────────────────
CASE_COLORS = {1: "#4CAF50", 2: "#2196F3", 3: "#FF9800"}
CASE_LABELS = {1: "Independent", 2: "Nested", 3: "Overlapping"}

def draw_genome_map(candidates, dna_len, title, out_file):
    fig, (ax_map, ax_density) = plt.subplots(
        2, 1, figsize=(18, 7),
        gridspec_kw={"height_ratios": [3, 1]}
    )
    fig.suptitle(title, fontsize=12, fontweight="bold")

    # ── Top panel: genome map with lane stacking ───────────────────────────────
    ax_map.barh(0, dna_len, left=0, height=0.25, color="#DDDDDD", zorder=1)
    ax_map.text(dna_len / 2, 0, f"Genome ({dna_len:,} bp)",
                ha="center", va="center", fontsize=7, color="#555", zorder=2)

    lanes = []
    def get_lane(s, e):
        for i, lane in enumerate(lanes):
            if all(e < ls or s > le for ls, le in lane):
                lane.append((s, e))
                return i + 0.5
        lanes.append([(s, e)])
        return len(lanes) - 1 + 0.5

    for i, te in enumerate(candidates):
        s, e = te["start"], te["end"]
        y = get_lane(s, e)
        color = CASE_COLORS[te["case"]]
        ax_map.barh(y, te["length"], left=s, height=0.35,
                    color=color, alpha=0.85, edgecolor="black", linewidth=0.5, zorder=3)
        ax_map.barh(y, IR_LEN, left=s, height=0.35,
                    color="black", alpha=0.3, zorder=4)
        ax_map.barh(y, IR_LEN, left=e - IR_LEN + 1, height=0.35,
                    color="black", alpha=0.3, zorder=4)
        # label if wide enough
        if te["length"] > dna_len * 0.012:
            ax_map.text(s + te["length"] / 2, y,
                        f"#{i+1}\n{te['length']}bp",
                        ha="center", va="center", fontsize=5,
                        fontweight="bold", color="white", zorder=5)

    max_y = max(len(lanes), 1) + 0.7
    ax_map.set_xlim(-500, dna_len + 500)
    ax_map.set_ylim(-0.4, max_y)
    ax_map.set_yticks([])
    ax_map.set_xlabel("Genome Position (bp)")
    ax_map.set_title("Detected TE Candidates", fontsize=10)
    ax_map.spines[["top", "right", "left"]].set_visible(False)

    legend_handles = [mpatches.Patch(color=c, label=CASE_LABELS[k])
                      for k, c in CASE_COLORS.items()]
    legend_handles.append(mpatches.Patch(color="black", alpha=0.3, label="Inverted Repeat"))
    ax_map.legend(handles=legend_handles, loc="upper right", fontsize=8, frameon=False)

    # ── Bottom panel: coverage density ────────────────────────────────────────
    bin_size = 500
    bins = range(0, dna_len + bin_size, bin_size)
    counts = [0] * (len(list(bins)) - 1)
    for te in candidates:
        for b in range(len(counts)):
            if te["start"] < (b + 1) * bin_size and te["end"] >= b * bin_size:
                counts[b] += 1

    colors_density = []
    for b, count in enumerate(counts):
        pos = b * bin_size
        case_at_pos = 1
        for te in candidates:
            if te["start"] <= pos <= te["end"]:
                case_at_pos = te["case"]
                break
        colors_density.append(CASE_COLORS.get(case_at_pos, "#4CAF50"))

    ax_density.bar(
        [b * bin_size for b in range(len(counts))],
        counts, width=bin_size * 0.9,
        color="#5C8DD6", alpha=0.7, edgecolor="none"
    )
    ax_density.set_xlim(-500, dna_len + 500)
    ax_density.set_xlabel("Genome Position (bp)")
    ax_density.set_ylabel("TE Count\nper 500bp", fontsize=8)
    ax_density.set_title("Candidate Density Along Genome", fontsize=10)
    ax_density.spines[["top", "right"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_file, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Graph saved to {out_file}")

# ── Main ───────────────────────────────────────────────────────────────────────
print("=" * 65)
print("  TRANSPOSABLE ELEMENT DETECTOR — Real Genome (.fna)")
print("=" * 65)

header, dna = read_fasta(FNA_FILE)
print(f"\nGenome : {header}")
print(f"Length : {len(dna):,} bp")
print(f"\nParameters:")
print(f"  IR length        : {IR_LEN} bp")
print(f"  Min body length  : {MIN_BODY} bp")
print(f"  Max TE length    : {MAX_TE_LEN} bp")
print(f"  Max IR frequency : {MAX_IR_OCC} occurrences\n")

candidates = find_transposons(dna, IR_LEN, MIN_BODY, MAX_TE_LEN, MAX_IR_OCC)
candidates = deduplicate(candidates, IR_LEN)
candidates = classify(candidates)

CASE_NAMES = {1: "independent", 2: "nested", 3: "overlapping"}

print(f"\n[ DETECTED TRANSPOSON CANDIDATES ]\n")
if not candidates:
    print("  No candidates found. Try relaxing parameters (larger MAX_IR_OCC or shorter IR_LEN).")
else:
    for i, te in enumerate(candidates, 1):
        print(f"  Candidate #{i}  —  {CASE_NAMES[te['case']]}")
        print(f"    Start     : {te['start']}")
        print(f"    End       : {te['end']}")
        print(f"    Length    : {te['length']} bp  (body: {te['body_len']} bp)")
        print(f"    IR        : 5'-{te['ir']}-3'")
        print(f"    RC-IR     : 5'-{te['rc_ir']}-3'")
        print()

print(f"Total candidates : {len(candidates)}")
print("=" * 65)

draw_genome_map(candidates, len(dna),
                f"Transposon Candidates — {header[:60]}",
                "transposons_genome.png")
