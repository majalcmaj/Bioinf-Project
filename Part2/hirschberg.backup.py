import numpy as np

from bio_data_parser import get_vector, get_matrix


class Hirschberg:
    def __init__(self, u, v, weights):
        self.u = self.dna_to_idx_vec(u)
        self.v = self.dna_to_idx_vec(v)
        self.weights = weights
        print(self.u, self.v, self.weights)

    DNA_TO_IDX = {
        'A': 0,
        'C': 1,
        'T': 2,
        'G': 3
    }

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
        for i in range(1, len(v_vec)):
            cols[i, 0] = cols[i - 1, 0] + self.weights[v_vec[i - 1], self.SPACE_INDEX]
        for i in range(0, len(u_vec)):
            self.nw_score(cols, u_vec[i], v_vec)
            cols[:, 0] = cols[:, 1]

    def _hirschberg(self, ui1, ui2, vi1, vi2):
        print("FROM", ui1, ui2, vi1, vi2, "\n")
        if ui2 - ui1 == 0:
            result = [(v_, ui2) for v_ in range(vi1, vi2)]
            return result
        elif vi2 - vi1 == 0:
            result = [(vi2, u_) for u_ in range(ui1, ui2)]
            return result
        elif ui2 - ui1 == 1:
            return ["SPC"]
        else:
            u_mid = (ui1 + ui2) // 2
            cols = np.zeros(shape=(vi2 - vi1 + 1, 2))
            self.calc_middle_weights(cols, self.u[ui1:u_mid], self.v[vi1:vi2])
            mid_col = np.copy(cols[:, 1])
            self.calc_middle_weights(cols, self.u[u_mid:ui2], self.v[vi2:vi1:-1])
            mid_col = mid_col + cols[::-1, 1]
            v_mid = vi1 + np.argmax(mid_col)
            print("LEFT ", ui1, u_mid, vi1, v_mid)
            print("RIGHT ", u_mid, ui2, v_mid, vi2)
            print("MID", (v_mid, u_mid))
            print("=" * 20)
            left_points = self._hirschberg(ui1, u_mid, vi1, v_mid)
            right_points = self._hirschberg(u_mid, ui2, v_mid, vi2)
            # Tu był błąd! W zlym momencie sklejalismy
            return left_points + right_points

    def run(self):
        path = self._hirschberg(0, len(self.u), 0, len(self.v))
        return path


if __name__ == "__main__":
    with open("data/similarity_maciek.data", 'r') as file:
        u_sym = get_vector(file)
        v_sym = get_vector(file)

        weights = get_matrix(file, 5)

    print(u_sym, v_sym)
    hirschberg = Hirschberg(u_sym, v_sym, weights)
    print(hirschberg.run())
