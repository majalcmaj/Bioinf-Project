import numpy as np

from dictionaries import dna_to_idx_vec, rna_to_idx_vec
from rna_coder import rna_code, rna_decode


def calculate_editional_distance(u_sym, v_sym, weights):
    matrix_dim = len(weights)

    if matrix_dim == 5:  # USING DNA
        u = dna_to_idx_vec(u_sym)
        v = dna_to_idx_vec(v_sym)
    elif matrix_dim == 21:  # USING RNA
        u_amino = rna_code(u_sym)
        v_amino = rna_code(v_sym)
        u = rna_to_idx_vec(u_amino)
        v = rna_to_idx_vec(v_amino)
    else:
        print("Matrix dimension is invalid!")
        return

    space_index = matrix_dim - 1
    h = len(u) + 1
    w = len(v) + 1

    costs = np.zeros([w, h], dtype=np.int32)
    transitions = np.zeros([w, h, 3], dtype=np.bool_)

    for x in range(1, w):
        costs[x, 0] = costs[x - 1, 0] + weights[space_index, v[x - 1]]
        transitions[0, x, 0] = True

    for y in range(1, h):
        costs[0, y] = costs[0, y - 1] + weights[u[y - 1], space_index]
        transitions[y, 0, 2] = True

    for y in range(1, h):
        for x in range(1, w):
            trans_costs = np.array([
                costs[x, y - 1] + weights[u[y - 1], space_index],  # left
                costs[x - 1, y - 1] + weights[u[y - 1], v[x - 1]],  # diag
                costs[x - 1, y] + weights[space_index, v[x - 1]]  # top
            ], dtype=np.int32)
            min_trans = np.min(trans_costs)

            idxs = []
            for idx in range(3):
                if trans_costs[idx] == min_trans:
                    idxs.append(idx)
            transitions[x, y][idxs] = True
            costs[x, y] = min_trans

    print("Editional distance: {}".format(costs[-1, -1]))

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

    print("".join(reversed(rna_decode(u_star))))
    print("".join(reversed(rna_decode(v_star))))
