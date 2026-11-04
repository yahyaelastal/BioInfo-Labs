import random
import math
import urllib.request
import matplotlib.pyplot as plt

ACCESSION = "L09137"
URL = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
       f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text")


def get_sequence():
    with urllib.request.urlopen(URL) as r:
        lines = r.read().decode().splitlines()
    return "".join(l.strip() for l in lines if not l.startswith(">")).upper()


def take_samples(seq, n=10, lo=100, hi=3000):
    samples = []
    for _ in range(n):
        length = random.randint(lo, min(hi, len(seq)))
        start = random.randint(0, len(seq) - length)
        samples.append(seq[start:start + length])
    return samples


def migration(size, gel_len=10.0, lo=100, hi=3000):
    frac = (math.log10(size) - math.log10(lo)) / (math.log10(hi) - math.log10(lo))
    return gel_len * (0.9 - 0.75 * frac)


def main():
    seq = get_sequence()
    print(f"Sequence {ACCESSION}: {len(seq)} nt")
    assert 1000 <= len(seq) <= 3000, "sequence must be 1000-3000 nt"

    samples = take_samples(seq)
    sizes = [len(s) for s in samples]
    ladder = [100, 250, 500, 1000, 1500, 2000, 3000]

    print("\nSample  Size(nt)  Distance(cm)")
    for i, s in enumerate(sorted(sizes, reverse=True), 1):
        print(f"{i:>6}  {s:>8}  {migration(s):>10.2f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6),
                                   gridspec_kw={"width_ratios": [1, 1.2]})

    ax1.set_facecolor("#1a1a1a")
    ax1.set_xlim(0, 4.2)
    ax1.set_ylim(10, 0)
    for x, name, data in [(1, "Ladder", ladder), (3, "Samples", sizes)]:
        ax1.add_patch(plt.Rectangle((x - 0.4, 0.2), 0.8, 0.4, color="#444"))
        ax1.text(x, -0.1, name, ha="center", va="bottom", color="black",
                 fontsize=11, fontweight="bold")
        last = -1
        for size in sorted(data, reverse=True):
            y = migration(size)
            ax1.hlines(y, x - 0.4, x + 0.4, colors="#7CFC00", linewidth=4)
            ly = max(y, last + 0.25)
            last = ly
            ax1.text(x + 0.45, ly, str(size), va="center", color="w", fontsize=8)
    ax1.set_title("Simulated agarose gel (– top, + bottom)", pad=25)
    ax1.set_xticks([])
    ax1.set_ylabel("Distance migrated (cm)")

    xs = list(range(100, 3001, 20))
    ax2.plot(xs, [migration(s) for s in xs], label="Migration curve")
    ax2.scatter(sizes, [migration(s) for s in sizes], color="red", zorder=3,
                label="Samples")
    ax2.set_xscale("log")
    ax2.set_xlabel("Fragment size (nt, log scale)")
    ax2.set_ylabel("Distance migrated (cm)")
    ax2.set_title("Size vs. migration distance", pad=25)
    ax2.grid(True, alpha=0.3)
    ax2.legend()

    plt.tight_layout()
    plt.savefig("ex1_gel.png", dpi=150)
    plt.show()


if __name__ == "__main__":
    main()