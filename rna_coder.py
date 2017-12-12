import itertools

CODONS_AMIN_MAP = {
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UGU": "C", "UGC": "C",
    "GAU": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "UUU": "F", "UUC": "F",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAU": "H", "CAC": "H",
    "AUU": "I", "AUC": "I", "AUA": "I",
    "AAA": "K", "AAG": "K",
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUG": "M",
    "AAU": "N", "AAC": "N",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    "ACU": "T", "ACC": "U", "ACA": "T", "ACG": "T",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UGG": "W",
    "UAU": "Y", "UAC": "Y",
    "---": "-"
}

AMIN_CODONS_MAP = {k: list(v)[0][0] for k, v in itertools.groupby(CODONS_AMIN_MAP.items(), lambda pair: pair[1])}


def rna_code(codons_array):
    assert len(codons_array) % 3 == 0, "The given codons array has to have the length divisible by 3"
    codon_triplets = ["".join(codons_array[idx:idx + 3]) for idx in range(0, len(codons_array), 3)]
    return list(map(lambda codon: CODONS_AMIN_MAP[codon], codon_triplets))


def rna_decode(amins_array):
    return [AMIN_CODONS_MAP[amin] for amin in amins_array]


if __name__ == "__main__":
    example_codons = "GUCGCAUGUCAGCUAGCUACCGAUUUCGACUGAUCGGGGCAU"
    example_aminos = "ACDEFG"

    print(rna_code(example_codons))
    print(len(example_codons) / 3, len(rna_code(example_codons)))

    print(rna_decode(example_aminos))
    print(len(example_aminos), len(rna_decode(example_aminos)))

    print("dictionairy = " , len(AMIN_CODONS_MAP))
