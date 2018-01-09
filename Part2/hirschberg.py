import numpy as np

from bio_data_parser import get_vector, get_matrix


class Hirschberg:
    def __init__(self, u, v, weights):
        self.u = self.dna_to_idx_vec(u)
        self.v = self.dna_to_idx_vec(v)
        self.weights = weights

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
        for i in range(1, len(cols)):
            cols[i, 1] = np.max([
                cols[i, 0] + self.weights[u, self.SPACE_INDEX],  # horizontal
                cols[i - 1, 0] + self.weights[u, v_vec[i - 1]],  # diagonal
                cols[i - 1, 1] + self.weights[v_vec[i - 1], self.SPACE_INDEX]
            ])

    def calc_middle_weights(self, cols, u_vec, v_vec):
        cols[0, 0] = 0
        for i in range(1, len(v_vec) + 1):
            cols[i, 0] = cols[i - 1, 0] + self.weights[v_vec[i - 1], self.SPACE_INDEX]
        for i in range(0, len(u_vec)):
            self.nw_score(cols, u_vec[i], v_vec)
            cols[:, 0] = cols[:, 1]

    def _hirschberg(self, u_idx, u_len, v_idx, v_len):
        if u_len == 1:
            return [(v_, u_idx) for v_ in range(v_idx, v_len)]
        elif v_len == 1:
            return [(v_idx, u_) for u_ in range(u_idx, u_len)]
        else:
            u_mid = u_idx + (u_len // 2)
            cols = np.zeros(shape=(v_len + 1, 2))
            self.calc_middle_weights(cols, self.u[u_idx:u_mid], self.v[v_idx:v_idx + v_len])
            mid_col = cols[:, 0]
            self.calc_middle_weights(cols, self.u[u_idx + u_len:u_mid - 1:-1], self.v[v_idx:v_idx + v_len][::-1])
            mid_col = mid_col + cols[::-1, 0]
            v_mid_local = np.argmax(mid_col)
            v_mid = v_mid_local + v_idx
            left_points = self._hirschberg(u_idx, u_len // 2, v_idx, v_mid_local)
            right_points = self._hirschberg(u_mid, u_len // 2, v_mid, v_len - v_mid_local)
            return left_points + [(v_mid, u_mid)] + right_points

    def run(self):
        return self._hirschberg(0, len(self.u), 0, len(self.v))


if __name__ == "__main__":
    with open("data/similarity.data", 'r') as file:
        u_sym = get_vector(file)
        v_sym = get_vector(file)

        weights = get_matrix(file, 5)

    hirschberg = Hirschberg(u_sym, v_sym, weights)
    print(hirschberg.run())
