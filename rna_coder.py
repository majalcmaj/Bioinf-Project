import itertools

CODONS_AMIN_MAP = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TGT": "C", "TGC": "C",
    "GAT": "D", "GAC": "D",
    "GAA": "E", "GAG": "E",
    "TTT": "F", "TTC": "F",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "AAA": "K", "AAG": "K",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATG": "M",
    "AAT": "N", "AAC": "N",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
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
    example_codons = "GTCGCATGTCAGCTAGCTACCGATTTCGACTGATCGGGGCAT"
    example_aminos = "ACDEFG"

    print(rna_code(example_codons))
    print(len(example_codons) / 3, len(rna_code(example_codons)))

    print(rna_decode(example_aminos))
    print(len(example_aminos), len(rna_decode(example_aminos)))

    print("dictionairy = " , len(AMIN_CODONS_MAP))
