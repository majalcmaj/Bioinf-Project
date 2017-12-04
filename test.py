from collections import namedtuple
import numpy as np

from parser import get_matrix, get_vector

ASCII_TO_IDX = {
    'A': 0,
    'C': 1,
    'T': 2,
    'G': 3,
    '-': 4
}


def vec_to_idx_vec(vec):
    return [ASCII_TO_IDX[char] for char in vec]


SPACE_INDEX = 4
if __name__ == "__main__":
    Node = namedtuple("Node", "weight isLeft isDiag isUp")
    # m = MyStruct("foo", "bar", "baz")

    file = open("biola.txt", 'r')
    u = vec_to_idx_vec(get_vector(file))
    v = vec_to_idx_vec(get_vector(file))

    weights = get_matrix(file)

    # print(u, v, weights)

    w, h = len(u) + 1, len(v) + 1

    costs = np.zeros([w, h], dtype=np.int32)
    transitions = np.zeros([w, h, 3], dtype=np.bool_)

    print(transitions.shape)

    for x in range(1, w):
        costs[x, 0] = costs[x-1, 0] + weights[SPACE_INDEX, u[x - 1]]

    for y in range(1, h):
        costs[0, y] = costs[0, y - 1] + weights[v[y - 1], SPACE_INDEX]

    for y in range(1, h):
        for x in range(1, w):
            trans_costs = np.array([
                costs[x, y-1] + weights[v[y - 1], SPACE_INDEX],  # top
                costs[x - 1, y] + weights[SPACE_INDEX, u[x - 1]],  # left
                costs[x - 1, y - 1] + weights[v[y - 1], u[x - 1]]  # diag
                ], dtype=np.int32)
            min_trans = np.min(trans_costs)

            idxs = []
            for idx in range(3):
                if trans_costs[idx] == min_trans:
                    idxs.append(idx)
            transitions[x, y][idxs] = True
            costs[x, y] = min_trans

    print(costs)
    #
    # # graph = np.empty(((w + 1) * (h + 1)), dtype=np.int32)
    # # graph.fill(np.iinfo(np.int32).max)
    #
    #



    # for y in range(h):
    #     for x in range(w):
    #         #graph[x][y] = weights[ascii_to_index[u[x]]][ascii_to_index[v[y]]]
    #         upValue   = weights[ascii_to_index(v[y])][ascii_to_index['U']] +
    #         leftValue = weights[ascii_to_index(u[x])][ascii_to_index['U']]
    #         diagValue = weights[ascii_to_index(v[y])][ascii_to_index[u[x]]]
    #
    #         graph[x][y] = Node()
    #
    # print(graph)
