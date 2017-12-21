import numpy as np

from Part2.bio_data_parser import get_vector, get_matrix

DNA_TO_IDX = {
    'A': 0,
    'C': 1,
    'T': 2,
    'G': 3
}

SPACE_INDEX = 4


def dna_to_idx_vec(vec):
    return [DNA_TO_IDX[char] for char in vec]


def nw_score(col, part_v, u, weights):
    result = np.zeros(col.shape, dtype=np.int32)
    result[0] = col[0] + weights[u, SPACE_INDEX]
    for i in range(1, len(part_v)):
        costs = [
            result[i - 1] + weights[part_v[i], SPACE_INDEX],
            col[i-1] + weights[part_v[i], u],
            col[i] + weights[u, SPACE_INDEX],
        ]
        result[i] = np.max(costs)
    return result


def calculate_similarity(u_sym, v_sym, weights):
    u = dna_to_idx_vec(u_sym)
    v = dna_to_idx_vec(v_sym)

    h = len(u) + 1
    w = len(v) + 1

    col = np.zeros(h, dtype=np.int32)
    for i in range(1, h):
        col[i] = col[i - 1] + weights[v[i-1], SPACE_INDEX]

    print(nw_score(col, v, u[1], weights))



if __name__ == "__main__":
    with open("data/similarity.data", 'r') as file:
        u_sym = get_vector(file)
        v_sym = get_vector(file)

        weights = get_matrix(file)

    calculate_similarity(u_sym, v_sym, weights)
