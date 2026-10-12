import matplotlib.pyplot as plt


def calculate_tm(sequence):
    sequence = sequence.upper().strip()

    if not sequence or any(base not in "ATGC" for base in sequence):
        raise ValueError("Invalid DNA sequence. Use only A, T, G, C.")

    gc = sequence.count("G") + sequence.count("C")
    at = sequence.count("A") + sequence.count("T")
    tm = 4 * gc + 2 * at
    return gc, at, tm


def plot_frequency(sequence):
    bases = ["A", "T", "G", "C"]
    counts = [sequence.count(b) for b in bases]

    plt.bar(bases, counts, color=["green", "red", "orange", "blue"])
    plt.title(f"Base Frequency in {sequence}")
    plt.xlabel("Base")
    plt.ylabel("Frequency")
    plt.yticks(range(0, max(counts) + 1))

    for i, c in enumerate(counts):
        plt.text(i, c + 0.05, str(c), ha="center")

    plt.show()


S = "TGACGTATA"

gc, at, tm = calculate_tm(S)

print("Sequence:", S)
print("G + C count:", gc)
print("A + T count:", at)
print("Tm = 4(G+C) + 2(A+T) =", tm, "°C")

plot_frequency(S)