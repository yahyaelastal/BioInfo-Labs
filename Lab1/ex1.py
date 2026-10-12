# Manual DNA sequence (20 characters)
dna = "ATGCTTAGGCTAACGTTGCA"
print("DNA Sequence:", dna)

# --- Single-letter frequencies ---
print("\nSingle Nucleotide Frequency:")
print("A:", dna.count("A"))
print("T:", dna.count("T"))
print("C:", dna.count("C"))
print("G:", dna.count("G"))


# --- Function to count combinations of length k ---
def combo_frequency(sequence, k):
    freq = {}
    for i in range(len(sequence) - k + 1):
        combo = sequence[i:i+k]
        freq[combo] = freq.get(combo, 0) + 1
    return dict(sorted(freq.items(), key=lambda x: x[1], reverse=True))


# Count 2‑letter and 3‑letter combos
freq2 = combo_frequency(dna, 2)
freq3 = combo_frequency(dna, 3)


# --- Simple text chart ---
def show_chart(title, data):
    print(f"\n{title}")
    for combo, count in data.items():
        bar = "█" * count
        print(f"{combo}: {bar} ({count})")


# Display charts
show_chart("2‑Letter Combination Frequency:", freq2)
show_chart("3‑Letter Combination Frequency:", freq3)
