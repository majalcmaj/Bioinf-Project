import numpy as np

DNA_TO_IDX = {
    'A': 0,
    'C': 1,
    'T': 2,
    'G': 3
}


def dna_to_idx_vec(vec):
    return [DNA_TO_IDX[char] for char in vec]


def penalty_function(length):
    return - (length ** 2)


def calculate_similarity(u_sym, v_sym, weights, p):
    u = dna_to_idx_vec(u_sym)
    v = dna_to_idx_vec(v_sym)

    h = len(u) + 1
    w = len(v) + 1

    min_int = np.iinfo(np.int32).min

    A = np.empty([w, h], dtype=np.int32).fill(min_int)
    B = np.empty([w, h], dtype=np.int32).fill(min_int)
    C = np.empty([w, h], dtype=np.int32).fill(min_int)
    S = np.empty([w, h], dtype=np.int32).fill(min_int)

    S[0, 0] = 0
    penalties_row = p(np.arange(1, w))
    S[0, 1:w] = penalties_row
    A[0, 1:w] = penalties_row

    penalties_col = p(np.arange(h, 1))
    S[1:h, 0] = penalties_col
    B[1:h, 0] = penalties_col

    for i in range(1, h):
        for j in range(1, w):
            A[i, j] = np.max(np.max(np.vstack((B[i, 0:j - 1], C[i, 0:j - 1])), axis=0) + p(np.arange(j, 0, -1)))
            B[i, j] = np.max(np.max(np.vstack((A[0:i - 1, j], C[0:i - 1, j])), axis=0) + p(np.arange(i, 0, -1)))
            C[i, j] = S[i - 1, j - 1] + weights[u[i], w[j]]
            S[i, j] = np.max([A[i, j], B[i, j], C[i, j]])

    print(S[h - 1, w - 1])

            #
            # costs = np.zeros([w, h], dtype=np.int32)
            # transitions = np.zeros([w, h, 3], dtype=np.bool_)
            #
            # for x in range(1, w):
            #     costs[x, 0] = costs[x - 1, 0] + weights[space_index, v[x - 1]]
            #     transitions[0, x, 0] = True
            #
            # for y in range(1, h):
            #     costs[0, y] = costs[0, y - 1] + weights[u[y - 1], space_index]
            #     transitions[y, 0, 2] = True
            #
            # for y in range(1, h):
            #     for x in range(1, w):
            #         if y == 3 and x == 3:
            #             a = 1
            #         trans_costs = np.array([
            #             costs[x, y - 1] + weights[u[y - 1], space_index],  # left
            #             costs[x - 1, y - 1] + weights[u[y - 1], v[x - 1]],  # diag
            #             costs[x - 1, y] + weights[space_index, v[x - 1]]  # top
            #         ], dtype=np.int32)
            #         max_trans = np.max(trans_costs)
            #
            #         idxs = []
            #         for idx in range(3):
            #             if trans_costs[idx] == max_trans:
            #                 idxs.append(idx)
            #         transitions[x, y][idxs] = True
            #         costs[x, y] = max_trans
            #
            # print("Similarity: {}".format(costs[-1, -1]))
            # print(costs)
            #
            # u_star = []
            # v_star = []
            #
            # x = w - 1
            # y = h - 1
            #
            # while x != 0 or y != 0:
            #     trans = transitions[y, x]
            #     if trans[0]:  # left
            #         x -= 1
            #         u_star.append(u_sym[x])
            #         v_star.append('-')
            #     elif trans[1]:  # diag
            #         x -= 1
            #         y -= 1
            #         u_star.append(u_sym[x])
            #         v_star.append(v_sym[y])
            #     elif trans[2]:  # top
            #         y -= 1
            #         u_star.append('-')
            #         v_star.append(v_sym[y])
            #
            # if as_codons:
            #     print("".join(rna_decode(reversed(u_star))))
            #     print("".join(rna_decode(reversed(v_star))))
            # else:
            #     print("".join(reversed(u_star)))
            #     print("".join(reversed(v_star)))

if __name__ == "__main__":
    pass