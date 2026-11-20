def needleman_wunsch(s1, s2, match=1, mismatch=-1, gap=-2):
    n, m = len(s1), len(s2)

    # Initialize scoring matrix
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(n + 1):
        dp[i][0] = i * gap
    for j in range(m + 1):
        dp[0][j] = j * gap

    # Fill the matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            score = match if s1[i-1] == s2[j-1] else mismatch
            dp[i][j] = max(
                dp[i-1][j-1] + score,   # diagonal
                dp[i-1][j]   + gap,      # up (gap in s2)
                dp[i][j-1]   + gap       # left (gap in s1)
            )

    # Traceback
    aligned1, aligned2 = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0:
            score = match if s1[i-1] == s2[j-1] else mismatch
            if dp[i][j] == dp[i-1][j-1] + score:
                aligned1.append(s1[i-1])
                aligned2.append(s2[j-1])
                i -= 1; j -= 1
                continue
        if i > 0 and dp[i][j] == dp[i-1][j] + gap:
            aligned1.append(s1[i-1])
            aligned2.append('-')
            i -= 1
        else:
            aligned1.append('-')
            aligned2.append(s2[j-1])
            j -= 1

    return (''.join(reversed(aligned1)),
            ''.join(reversed(aligned2)),
            dp)


def print_matrix(dp, s1, s2):
    header = "      " + "  ".join(f"{c:3}" for c in s2)
    print(header)
    for i, row in enumerate(dp):
        label = s1[i-1] if i > 0 else " "
        print(f"{label:2} " + "  ".join(f"{v:3}" for v in row))


if __name__ == "__main__":
    s1 = "ATCCGTTAC"
    s2 = "CGAT"

    print(f"Sequence 1: {s1}")
    print(f"Sequence 2: {s2}")
    print(f"Scoring : match=+1, mismatch=-1, gap=-2\n")

    a1, a2, dp = needleman_wunsch(s1, s2)

    print("Scoring Matrix:")
    print_matrix(dp, s1, s2)

    match_line = ''.join('|' if a == b else ' ' for a, b in zip(a1, a2))

    print("\nAlignment:")
    print(f"  {a1}")
    print(f"  {match_line}")
    print(f"  {a2}")
    print(f"\nAlignment score: {dp[len(s1)][len(s2)]}")
