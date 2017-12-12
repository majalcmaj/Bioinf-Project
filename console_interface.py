import numpy as np

from bio_data_parser import get_matrix, get_vector
from editional_distance import calculate_editional_distance
from local_match import calculate_local_match
from similarity import calculate_similarity

DEFAULT_FILENAMES = {
    11: "editional_distance.data",
    12: "similarity.data",
    13: "local_match.data",
    21: "similarity_rna.data",
    22: "similarity_rna.data",
    23: "similarity_rna.data"  # todo bedziemy chcieli rozne pliki? Mikel mnie zabije za ten slownik XD
}

if __name__ == "__main__":

    while True:
        mode = 0
        while mode != 1 and mode != 2:
            mode = int(input("-------------\n"
                             "Wybierz tryb:\n"
                             "1. DNA\n"
                             "2. RNA\n"))

        task = 0
        while task != 1 and task != 2 and task != 3:
            task = int(input("Wybierz zadanie:\n"
                             "1. Odleglosc edycyjna\n"
                             "2. Podobienstwo\n"
                             "3. Najlepsze zestawienie lokalne\n"))

        file_name = (input("Podaj nazwe pliku: (zostaw puste - plik domyslny)") or DEFAULT_FILENAMES[mode * 10 + task])
        file = open("data/" + file_name, 'r')
        u_sym = get_vector(file)
        v_sym = get_vector(file)

        weights = get_matrix(file, mode)
        print("Wczytano plik: " + file_name)

        if task == 1:  # BIEDA SWITCH XDD
            calculate_editional_distance(u_sym, v_sym, weights)
        if task == 2:
            calculate_similarity(u_sym, v_sym, weights)
        if task == 3:
            calculate_local_match(u_sym, v_sym, weights)
