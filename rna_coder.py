CODONS_TO_AMINO = {
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
    # "TAG": "O",  # duplicate from "|" , add this to decode!
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # "TGA": "U",  # duplicate from "|"
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "TAA": "|", "TAG": "|", "TGA": "|"
}


def rna_code(codons_array):  # input array of codons (len%3==0), return array of amino acids
    aminos_array = []
    codons_copy = codons_array  # do I need to do this? On the Internet they say I do
    while len(codons_copy) >= 3:  # if there is 2 or less characters left just ignore them
        amino = CODONS_TO_AMINO.get(codons_copy[:3])
        # amino = (amino != None) ? amino : "!" # ????????? IN JAVA IT WORKS XD
        if amino is not None:
            aminos_array.append(amino)
        else:
            aminos_array.append("!")  # there was wrong codons triplet, write "!" instead of amino acid
        codons_copy = codons_copy.replace(codons_copy[:3], '', 1)  # delete first 3 chars
    return aminos_array


def rna_decode(aminos_array):
    aminos_copy = aminos_array  # do I need to copy input?
    codons_array = []
    while aminos_copy != '':
        next_amino = aminos_copy[:1]
        for found_codon_triplet, search_amino in CODONS_TO_AMINO.items():
            if next_amino == search_amino:  # TODO: what if key not found , and special cases: "O" and "U"
                codons_array.append(found_codon_triplet)
                break
        aminos_copy = aminos_copy.replace(aminos_copy[:1], '', 1)  # delete first character
    return codons_array


if __name__ == "__main__":
    example_codons = "GTCGCATGTCAGCTAGCTACCGATTTCGACTGATCGGGGCAT"
    example_aminos = "ACDEFG"

    print(rna_code(example_codons))
    print(len(example_codons)/3, len(rna_code(example_codons)))

    print(rna_decode(example_aminos))
    print(len(example_aminos), len(rna_decode(example_aminos)))