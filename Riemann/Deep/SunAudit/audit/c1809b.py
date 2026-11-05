import sys, cmath, math
sys.path.insert(0, '.')
from lib import *

def det_bareiss(M):
    M = [row[:] for row in M]
    n = len(M)
    if n == 0:
        return 1
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = None
            for i in range(k + 1, n):
                if M[i][k] != 0:
                    sw = i; break
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

# self-test against Leibniz for small random matrices
import itertools, random
random.seed(1)
def det_leibniz(M):
    n = len(M); s = 0
    for perm in itertools.permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i+1, n) if perm[i] > perm[j])
        t = 1
        for i in range(n): t *= M[i][perm[i]]
        s += (-1) ** inv * t
    return s
for _ in range(300):
    n = random.randint(1, 5)
    M = [[random.randint(-3, 3) for _ in range(n)] for _ in range(n)]
    assert det_bareiss(M) == det_leibniz(M), M

def R(k, n):
    r = k % n
    return min(r, n - r)

res14 = []; res15 = []
for n in range(3, 160, 2):
    m = (n - 1) // 2
    A = [[R(i * i * j * j, n) for j in range(1, m + 1)] for i in range(1, m + 1)]
    B = [[(i * i * j * j) // n for j in range(1, m + 1)] for i in range(1, m + 1)]
    dA = det_bareiss(A); dB = det_bareiss(B)
    predA = is_prime(n) and n % 4 == 3
    predB = n == 9 or (is_prime(n) and n > 7 and n % 4 == 3)
    if (dA != 0) != predA: res14.append((n, dA))
    if (dB != 0) != predB: res15.append((n, dB))
print('6.8 (6.14) failures n<160:', res14)
print('6.8 (6.15) failures n<160:', res15)
m=3
print('n=7 B matrix', [[(i*i*j*j)//7 for j in range(1,4)] for i in range(1,4)], det_bareiss([[(i*i*j*j)//7 for j in range(1,4)] for i in range(1,4)]))
