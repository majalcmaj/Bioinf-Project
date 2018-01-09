import numpy as np

from bio_data_parser import get_vector, get_matrix


class Hirschberg:
    def __init__(self, u, v, weights):
        if len(u) < len(v):
            u, v = v, u
        self.u = self.dna_to_idx_vec(u)
        self.v = self.dna_to_idx_vec(v)
        self.weights = weights

    DNA_TO_IDX = {
        'A': 0,
        'C': 1,
        'T': 2,
        'G': 3
    }

    IDX_TO_DNA = "ACTG-"

    SPACE_INDEX = 4

    @staticmethod
    def dna_to_idx_vec(vec):
        return [Hirschberg.DNA_TO_IDX[char] for char in vec]

    def nw_score(self, cols, u, v_vec):
        cols[0, 1] = cols[0, 0] + self.weights[u, self.SPACE_INDEX]
        for i in range(0, len(v_vec)):
            cols[i + 1, 1] = np.max([
                cols[i + 1, 0] + self.weights[u, self.SPACE_INDEX],  # horizontal
                cols[i, 0] + self.weights[u, v_vec[i]],  # diagonal
                cols[i, 1] + self.weights[v_vec[i], self.SPACE_INDEX]
            ])

    def calc_middle_weights(self, cols, u_vec, v_vec):
        cols[0, 0] = 0
        for i in range(1, len(v_vec) + 1):
            cols[i, 0] = cols[i - 1, 0] + self.weights[v_vec[i - 1], self.SPACE_INDEX]
        for i in range(0, len(u_vec)):
            self.nw_score(cols, u_vec[i], v_vec)
            cols[:, 0] = cols[:, 1]

    def similarity_1_col(self, ui, vi1, vi2):
        char = self.u[ui]
        memo = vi1

        for i in range(vi1 + 1, vi2):
            if self.v[i] == char:
                memo = i

        result = []
        for i in range(vi1, vi2):
            if i == memo:
                char = self.u[ui]
            else:
                char = self.SPACE_INDEX
            result.append((self.v[i], char))
        return result
        # cols = np.zeros((v_len + 1, 2))
        # cols[0, 1] = cols[0, 0] + self.weights[self.u[ui - 1], self.SPACE_INDEX]
        # transitions = np.zeros((v_len + 1, 2, 3), dtype=np.bool)
        # transitions[0, 1, 0] = True
        # transitions[:, 0, 2] = True
        # for i in range(1, v_len + 1):
        #     cols[i, 0] = cols[i - 1, 0] + self.weights[self.v[vi], self.SPACE_INDEX]
        #     trans_costs = [
        #         cols[i, 0] + self.weights[self.u[ui - 1], self.SPACE_INDEX],  # horizontal
        #         cols[i - 1, 0] + self.weights[self.u[ui - 1], self.v[vi + i - 1]],  # diagonal
        #         cols[i - 1, 1] + self.weights[self.v[vi + i - 1], self.SPACE_INDEX]
        #     ]
        #     cols[i, 1] = max_val = np.max(trans_costs)
        #
        #     idxs = []
        #     for idx in range(3):
        #         if trans_costs[idx] == max_val:
        #             idxs.append(idx)
        #     transitions[i, 1][idxs] = True
        #
        # x = 1
        # y = len(cols) - 1
        #
        # result = []
        #
        # while x != 0 or y != 0:
        #     if x < 0 or y < 0:
        #         raise Exception()
        #     trans = transitions[y, x]
        #     if trans[0]:  # left
        #         x -= 1
        #     elif trans[1]:  # diag
        #         x -= 1
        #         y -= 1
        #     elif trans[2]:  # top
        #         y -= 1
        #     result.append((vi + y, ui + x))
        # return list(reversed(result))

    def _hirschberg(self, ui1, ui2, vi1, vi2, L):
        print("  " * L, ui1, ui2, vi1, vi2)

        if ui2 - ui1 == 0:
            result = [(self.v[v_], self.SPACE_INDEX) for v_ in range(vi1, vi2)]
            print("  " * L, "| ",
                  [(self.IDX_TO_DNA[self.v[v_]], self.IDX_TO_DNA[self.SPACE_INDEX]) for v_ in range(vi1, vi2)])
            return result, 0
        elif vi2 - vi1 == 0:
            result = [(self.SPACE_INDEX, self.u[u_]) for u_ in range(ui1, ui2)]
            print("  " * L, "- ",
                  [(self.IDX_TO_DNA[self.SPACE_INDEX], self.IDX_TO_DNA[self.u[u_]]) for u_ in range(ui1, ui2)])
            return result, 0
        elif ui2 - ui1 == 1:
            print("  " * L, "|\| ", [(self.IDX_TO_DNA[i1], self.IDX_TO_DNA[self.u[i2]]) for i1, i2 in
                                     self.similarity_1_col(ui1, vi1, vi2)])
            # return ["..."]
            return self.similarity_1_col(ui1, vi1, vi2), 0
        else:
            if ui1 == 3 and ui2 == 6 and vi1 == 4 and vi2 == 6:
                print("")
            u_mid = (ui1 + ui2) // 2
            cols = np.zeros(shape=(vi2 - vi1 + 1, 2))
            self.calc_middle_weights(cols, self.u[ui1:u_mid], self.v[vi1:vi2])
            mid_col = np.copy(cols[:, 1])
            self.calc_middle_weights(cols, self.u[u_mid:ui2][::-1], self.v[vi1:vi2][::-1])
            mid_col = mid_col + cols[::-1, 1]
            v_mid = vi1 + np.argmax(mid_col)
            # print("LEFT ", ui1, u_mid, vi1, v_mid)
            # print("RIGHT ", u_mid, ui2, v_mid, vi2)
            # print("=" * 20)
            left_points, _ = self._hirschberg(ui1, u_mid, vi1, v_mid, L + 1)
            right_points, _ = self._hirschberg(u_mid, ui2, v_mid, vi2, L + 1)
            return left_points + right_points, np.max(mid_col)

    def run(self):
        path, score = self._hirschberg(0, len(self.u), 0, len(self.v), 0)
        strs = list(zip(*path))
        return "".join([self.IDX_TO_DNA[idx] for idx in strs[0]]), "".join(
            [self.IDX_TO_DNA[idx] for idx in strs[1]]), str(score)


if __name__ == "__main__":
    with open("data/similarity2.data", 'r') as file:
        u_sym = get_vector(file)
        v_sym = get_vector(file)

        weights = get_matrix(file, 5)

    hirschberg = Hirschberg(u_sym, v_sym, weights)
    print("\n".join(hirschberg.run()))
