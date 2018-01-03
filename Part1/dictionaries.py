DNA_TO_IDX = {
    'A': 0,
    'C': 1,
    'T': 2,
    'G': 3
}

RNA_TO_IDX = {
    'A': 0,
    'C': 1,
    'D': 2,
    'E': 3,
    'F': 4,
    'G': 5,
    'H': 6,
    'I': 7,
    'K': 8,
    'L': 9,
    'M': 10,
    'N': 11,
    'P': 12,
    'Q': 13,
    'R': 14,
    'S': 15,
    'T': 16,
    'V': 17,
    'W': 18,
    'Y': 19,
    '-': 20
}


def dna_to_idx_vec(vec):
    return [DNA_TO_IDX[char] for char in vec]


def rna_to_idx_vec(vec):
    return [RNA_TO_IDX[char] for char in vec]