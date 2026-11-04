import random
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

BASES = ['A', 'T', 'G', 'C']

def random_dna(length):
    return ''.join(random.choices(BASES, k=length))

def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp[b] for b in reversed(seq))

def build_transposon(body_len=20, ir_len=6):
    """Build a transposon: [IR][body][RC(IR)]"""
    ir = random_dna(ir_len)
    body = random_dna(body_len)
    rc_ir = reverse_complement(ir)
    return ir + body + rc_ir, ir, rc_ir

def insert(dna, pos, sequence):
    return dna[:pos] + sequence + dna[pos:]

# ── Generate base DNA ──────────────────────────────────────────────────────────
random.seed(42)
dna = random_dna(300)

known_irs = []     # list of (ir, rc_ir) for every inserted transposon
transposons = []
ir_len = 6

print("=" * 60)
print("  TRANSPOSABLE ELEMENT SIMULATOR & DETECTOR")
print("=" * 60)

# ── Case 1: two independent transposons ───────────────────────────────────────
te1, ir1, rc_ir1 = build_transposon(body_len=20, ir_len=ir_len)
pos1 = 30
dna = insert(dna, pos1, te1)
known_irs.append((ir1, rc_ir1))
transposons.append({
    "name": "TE1 (independent)", "case": 1,
    "start": pos1, "end": pos1 + len(te1) - 1,
    "length": len(te1), "ir": ir1, "rc_ir": rc_ir1,
})

te2, ir2, rc_ir2 = build_transposon(body_len=18, ir_len=ir_len)
pos2 = 120
dna = insert(dna, pos2, te2)
known_irs.append((ir2, rc_ir2))
transposons.append({
    "name": "TE2 (independent)", "case": 1,
    "start": pos2, "end": pos2 + len(te2) - 1,
    "length": len(te2), "ir": ir2, "rc_ir": rc_ir2,
})

# ── Case 2: nested — TE4 inside TE3 ───────────────────────────────────────────
te3, ir3, rc_ir3 = build_transposon(body_len=40, ir_len=ir_len)
pos3 = 200
dna = insert(dna, pos3, te3)
known_irs.append((ir3, rc_ir3))

te4, ir4, rc_ir4 = build_transposon(body_len=14, ir_len=ir_len)
inner_offset = ir_len + 10
pos4 = pos3 + inner_offset
dna = insert(dna, pos4, te4)
known_irs.append((ir4, rc_ir4))

transposons.append({
    "name": "TE3 (outer — nested)", "case": 2,
    "start": pos3, "end": pos3 + len(te3) + len(te4) - 1,
    "length": len(te3) + len(te4), "ir": ir3, "rc_ir": rc_ir3,
    "note": f"contains TE4 at offset {inner_offset}",
})
transposons.append({
    "name": "TE4 (inner — nested)", "case": 2,
    "start": pos4, "end": pos4 + len(te4) - 1,
    "length": len(te4), "ir": ir4, "rc_ir": rc_ir4,
    "note": "nested inside TE3",
})

# ── Case 3: overlapping ───────────────────────────────────────────────────────
te5, ir5, rc_ir5 = build_transposon(body_len=22, ir_len=ir_len)
pos5 = 380
dna = insert(dna, pos5, te5)
known_irs.append((ir5, rc_ir5))

overlap = 8
pos6 = pos5 + len(te5) - overlap
te6, ir6, rc_ir6 = build_transposon(body_len=20, ir_len=ir_len)
dna = insert(dna, pos6, te6)
known_irs.append((ir6, rc_ir6))

transposons.append({
    "name": "TE5 (overlapping)", "case": 3,
    "start": pos5, "end": pos5 + len(te5) - 1,
    "length": len(te5), "ir": ir5, "rc_ir": rc_ir5,
    "note": f"overlaps with TE6 by {overlap} bases",
})
transposons.append({
    "name": "TE6 (overlapping)", "case": 3,
    "start": pos6, "end": pos6 + len(te6) - 1,
    "length": len(te6), "ir": ir6, "rc_ir": rc_ir6,
    "note": f"overlaps with TE5 by {overlap} bases",
})

# ── Detection — search only for known IR sequences ────────────────────────────
def detect_transposons(dna, known_irs, ir_len=6):
    """
    For each known (IR, RC-IR) pair, find the IR in the DNA then
    look for its RC-IR downstream — that pair is a transposon boundary.
    """
    found = []
    for ir, rc_ir in known_irs:
        # find every occurrence of this IR
        start = 0
        while True:
            pos = dna.find(ir, start)
            if pos == -1:
                break
            # look for the closing RC-IR after the IR
            close = dna.find(rc_ir, pos + ir_len)
            if close != -1:
                end = close + ir_len - 1
                found.append({
                    "start": pos, "end": end,
                    "length": end - pos + 1,
                    "ir": ir, "rc_ir": rc_ir,
                })
            start = pos + 1
    found.sort(key=lambda x: x["start"])
    return found

detected = detect_transposons(dna, known_irs, ir_len)

# ── Classify detected cases ────────────────────────────────────────────────────
def classify_case(d, all_detected):
    """Assign case 1/2/3 by checking spatial relationships."""
    s, e = d["start"], d["end"]
    for other in all_detected:
        if other is d:
            continue
        os, oe = other["start"], other["end"]
        if os < s and oe > e:          # d is fully inside other → nested inner
            return 2
        if s < os and e > oe:          # d fully contains other → nested outer
            return 2
        if s < os < e < oe:            # partial overlap
            return 3
        if os < s < oe < e:
            return 3
    return 1

for d in detected:
    d["case"] = classify_case(d, detected)

# ── Print ground truth ────────────────────────────────────────────────────────
print("\n[ GROUND TRUTH — Inserted Transposons ]\n")
for t in transposons:
    print(f"  {t['name']}")
    print(f"    Case   : {t['case']}")
    print(f"    Start  : {t['start']}")
    print(f"    End    : {t['end']}")
    print(f"    Length : {t['length']} bp")
    print(f"    IR     : 5'-{t['ir']}-3'")
    print(f"    RC-IR  : 5'-{t['rc_ir']}-3'")
    if "note" in t:
        print(f"    Note   : {t['note']}")
    print()

# ── Print detections ──────────────────────────────────────────────────────────
CASE_NAMES = {1: "independent", 2: "nested", 3: "overlapping"}
print("[ DETECTED Transposons ]\n")
for i, d in enumerate(detected, 1):
    print(f"  Detection #{i}  —  {CASE_NAMES[d['case']]}")
    print(f"    Start  : {d['start']}")
    print(f"    End    : {d['end']}")
    print(f"    Length : {d['length']} bp")
    print(f"    IR     : 5'-{d['ir']}-3'")
    print(f"    RC-IR  : 5'-{d['rc_ir']}-3'")
    print()

print(f"Total DNA length : {len(dna)} bp")
print(f"Transposons found: {len(detected)}")
print("=" * 60)

# ── Visualization ──────────────────────────────────────────────────────────────
CASE_COLORS = {1: "#4CAF50", 2: "#2196F3", 3: "#FF9800"}
CASE_LABELS = {1: "Independent", 2: "Nested", 3: "Overlapping"}

fig, (ax_truth, ax_det) = plt.subplots(2, 1, figsize=(14, 6))
fig.suptitle("Transposable Elements in Artificial DNA", fontsize=14, fontweight="bold")
dna_len = len(dna)

def draw_panel(ax, elements, title, numbered=False):
    ax.barh(0, dna_len, left=0, height=0.3, color="#CCCCCC", zorder=1)
    ax.text(dna_len / 2, 0, f"DNA ({dna_len} bp)",
            ha="center", va="center", fontsize=8, color="#444", zorder=2)

    # assign vertical lanes to avoid overlap in display
    lanes = []
    def get_lane(start, end):
        for i, lane in enumerate(lanes):
            if all(end < s or start > e for s, e in lane):
                lane.append((start, end))
                return i + 0.6
        lanes.append([(start, end)])
        return len(lanes) - 1 + 0.6

    for idx, elem in enumerate(elements):
        s, e = elem["start"], elem["end"]
        length = elem["length"]
        case = elem.get("case", 1)
        y = get_lane(s, e)
        color = CASE_COLORS[case]

        ax.barh(y, length, left=s, height=0.38,
                color=color, alpha=0.88, edgecolor="black", linewidth=0.8, zorder=3)
        # IR blocks
        ax.barh(y, ir_len, left=s, height=0.38, color="black", alpha=0.35, zorder=4)
        ax.barh(y, ir_len, left=e - ir_len + 1, height=0.38, color="black", alpha=0.35, zorder=4)

        label = elem.get("name", f"Det#{idx+1}").split("(")[0].strip() if not numbered \
                else f"Det#{idx+1}"
        ax.text(s + length / 2, y,
                f"{label}\n{s}–{e}\n{length} bp",
                ha="center", va="center", fontsize=6.5,
                fontweight="bold", color="white", zorder=5)

    max_y = len(lanes) + 0.9 if lanes else 1.2
    ax.set_xlim(-15, dna_len + 15)
    ax.set_ylim(-0.4, max_y)
    ax.set_yticks([])
    ax.set_xlabel("Position (bp)")
    ax.set_title(title, fontsize=11)
    ax.spines[["top", "right", "left"]].set_visible(False)

draw_panel(ax_truth, transposons, "Ground Truth (Inserted Transposons)", numbered=False)
draw_panel(ax_det,   detected,    "Detected Transposons (IR Scanning)",  numbered=True)

legend_handles = [mpatches.Patch(color=c, label=CASE_LABELS[k]) for k, c in CASE_COLORS.items()]
legend_handles.append(mpatches.Patch(color="black", alpha=0.35, label="Inverted Repeat"))
fig.legend(handles=legend_handles, loc="lower center", ncol=4,
           fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))

plt.tight_layout(rect=[0, 0.05, 1, 1])
plt.savefig("transposons.png", dpi=150, bbox_inches="tight")
plt.show()
print("\nGraph saved to transposons.png")
