import numpy as np

from dictionaries import dna_to_idx_vec, rna_to_idx_vec
from rna_coder import rna_code, rna_decode


def calculate_local_match(u_sym, v_sym, weights):
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

    for y in range(1, h):
        for x in range(1, w):
            trans_costs = np.array([
                costs[x, y - 1] + weights[u[y - 1], space_index],  # left
                costs[x - 1, y - 1] + weights[u[y - 1], v[x - 1]],  # diag
                costs[x - 1, y] + weights[space_index, v[x - 1]],  # top
                0
            ], dtype=np.int32)
            max_trans = np.max(trans_costs)

            if max_trans != 0:
                idxs = []
                for idx in range(4):
                    if trans_costs[idx] == max_trans:
                        idxs.append(idx)

                transitions[x, y][idxs] = True
            costs[x, y] = max_trans

    u_star = []
    v_star = []

    flat_max = np.argmax(costs)
    y = int(flat_max / w)
    x = flat_max % w
    print("Best local match: {}".format(costs[y, x]))

    while costs[y, x] != 0:
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
