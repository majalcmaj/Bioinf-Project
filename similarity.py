import numpy as np

from bio_data_parser import get_matrix, get_vector

ASCII_TO_IDX = {
    'A': 0,
    'C': 1,
    'T': 2,
    'G': 3
}


def vec_to_idx_vec(vec):
    return [ASCII_TO_IDX[char] for char in vec]


SPACE_INDEX = 4
if __name__ == "__main__":

    file = open("data/similarity.data", 'r')
    u_sym = get_vector(file)
    v_sym = get_vector(file)
    u = vec_to_idx_vec(u_sym)
    v = vec_to_idx_vec(v_sym)

    weights = get_matrix(file)

    h = len(u) + 1
    w = len(v) + 1

    costs = np.zeros([w, h], dtype=np.int32)
    transitions = np.zeros([w, h, 3], dtype=np.bool_)

    for x in range(1, w):
        costs[x, 0] = costs[x - 1, 0] + weights[SPACE_INDEX, v[x - 1]]
        transitions[0, x, 0] = True

    for y in range(1, h):
        costs[0, y] = costs[0, y - 1] + weights[u[y - 1], SPACE_INDEX]
        transitions[y, 0, 2] = True

    for y in range(1, h):
        for x in range(1, w):
            trans_costs = np.array([
                costs[x, y - 1] + weights[u[y - 1], SPACE_INDEX],  # left
                costs[x - 1, y - 1] + weights[u[y - 1], v[x - 1]],  # diag
                costs[x - 1, y] + weights[SPACE_INDEX, v[x - 1]]  # top
            ], dtype=np.int32)
            max_trans = np.max(trans_costs)

            idxs = []
            for idx in range(3):
                if trans_costs[idx] == max_trans:
                    idxs.append(idx)
            transitions[x, y][idxs] = True
            costs[x, y] = max_trans

    print("Similarity: {}".format(costs[-1, -1]))

    u_star = []
    v_star = []

    x = w - 1
    y = h - 1

    while x != 0 or y != 0:
        trans = transitions[y, x]
        if trans[0]:  # left
            x -= 1
            u_star.append(u_sym[x])
            v_star.append('-')
        elif trans[1]:  # diag
            x -= 1
            y -= 1
            u_star.append(u_sym[x])
            v_star.append(v_sym[y])
        elif trans[2]:  # top
            y -= 1
            u_star.append('-')
            v_star.append(v_sym[y])

    print("".join(reversed(u_star)))
    print("".join(reversed(v_star)))
