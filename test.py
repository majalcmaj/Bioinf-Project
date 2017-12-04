from Scripts.parser import get_matrix, get_vector;
from enum import Enum

"""class Index(Enum):
    A = 0
    C = 1
    T = 2
    G = 3
    U = 4"""


ascii_to_index = {'A': 0,
           'C': 1,
           'T': 2,
           'G': 3,
           'U': 4,
           }

file = open("biola.txt",'r');
u = get_vector(file);
v = get_vector(file);
weights = get_matrix(file);


print(u,v,weights);


w, h = len(u), len(v);
transit = [[0 for x in range(w)] for y in range(h)]

for y in range(h):
    for x in range(w):
        transit[y][x] = weights[ascii_to_index[u[x]]][ascii_to_index[v[y]]]

print(transit);



