import re
import numpy as np


def get_vector(file):
    vec = file.readline()
    vec = re.sub("actguACTGU", '', vec)
    vec = vec.upper()
    vec = vec.split('[')[1].split(']')[0]
    vec = vec.split()
    return vec


def get_matrix(file, mode):
    if mode == 1:
        dim = 5
    elif mode == 2:
        dim = 21
    else:
        print("Invalid mode")
        return

    matrix = [[np.char for x in range(dim)] for y in range(dim)]
    file.readline()
    for i in range(0, dim):
        line = file.readline()
        matrix[i][:] = line.split()

    matrix = np.array(matrix)
    matrix = matrix.astype(int)
    return matrix
