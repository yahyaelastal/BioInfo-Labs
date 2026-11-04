import urllib.request
from collections import Counter
import matplotlib.pyplot as plt

ACCESSION = "L09137"
URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text")


def get_sequence():
    with urllib.request.urlopen(URL) as r:
        lines = r.read().decode().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">")).upper()


def is_primitive(unit):
    n = len(unit)
    for d in range(1, n):
        if n % d == 0 and unit[:d] * (n // d) == unit:
            return False
    return True


def find_repeats(seq, unit_min=3, unit_max=5, total_min=6, total_max=10):
    found = []
    n = len(seq)
    for u in range(unit_min, unit_max + 1):
        i = 0
        while i + u <= n:
            unit = seq[i:i + u]
            if not is_primitive(unit):
                i += 1
                continue
            if i >= u and seq[i - u:i] == unit:
                i += 1
                continue
            count = 1
            j = i + u
            while j + u <= n and seq[j:j + u] == unit:
                count += 1
                j += u
            total = count * u
            if count >= 2 and total_min <= total <= total_max:
                found.append({"start": i + 1, "end": i + total, "unit": unit,
                              "copies": count, "length": total})
            i += 1
    found.sort(key=lambda r: r["start"])
    return found


def main():
    seq = get_sequence()
    print(f"Sequence {ACCESSION}: {len(seq)} nt")
    assert 1000 <= len(seq) <= 3000, "sequence must be 1000-3000 nt"

    repeats = find_repeats(seq)

    print(f"\nFound {len(repeats)} repeats\n")
    print(f"{'Start':>6} {'End':>6}  {'Unit':<6} {'Copies':>6} {'Length':>6}  Repeat")
    for r in repeats:
        print(f"{r['start']:>6} {r['end']:>6}  {r['unit']:<6} {r['copies']:>6} "
              f"{r['length']:>6}  {r['unit'] * r['copies']}")

    by_unit_len = Counter(len(r["unit"]) for r in repeats)
    by_motif = Counter(r["unit"] for r in repeats).most_common(10)
    colors = {3: "#1f77b4", 4: "#ff7f0e", 5: "#2ca02c"}

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])

    units = [3, 4, 5]
    ax1.bar([str(u) for u in units], [by_unit_len.get(u, 0) for u in units],
            color=[colors[u] for u in units])
    ax1.set_xlabel("Repeated unit length (bases)")
    ax1.set_ylabel("Number of repeats")
    ax1.set_title("Repeats by unit length")
    for k, u in enumerate(units):
        ax1.text(k, by_unit_len.get(u, 0), str(by_unit_len.get(u, 0)),
                 ha="center", va="bottom")

    if by_motif:
        names = [m for m, _ in by_motif]
        vals = [c for _, c in by_motif]
        ax2.barh(names[::-1], vals[::-1], color="#7f7f7f")
    ax2.set_xlabel("Occurrences")
    ax2.set_title("Most frequent repeated units")

    for u in units:
        pts = [r for r in repeats if len(r["unit"]) == u]
        ax3.scatter([r["start"] for r in pts], [u] * len(pts),
                    color=colors[u], s=40, label=f"unit = {u} bases")
    ax3.set_xlim(0, len(seq))
    ax3.set_yticks(units)
    ax3.set_ylim(2.5, 5.5)
    ax3.set_xlabel("Position in sequence (nt)")
    ax3.set_ylabel("Unit length")
    ax3.set_title("Position of repeats along the sequence")
    ax3.grid(True, alpha=0.3)
    ax3.legend(loc="upper right")

    fig.suptitle(f"Tandem repeats in {ACCESSION} ({len(seq)} nt)", fontsize=14)
    plt.tight_layout()
    plt.savefig("ex2_repeats.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()