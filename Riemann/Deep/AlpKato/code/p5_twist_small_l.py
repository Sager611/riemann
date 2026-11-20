# P5 (alpkato): KatoTwoW's second conjunct is the dual: for 2x2 W with det W = theta, theta^2 = 1,
#   W(g) (x) theta(g) = det(W(g)) W(g)  and  J W(g) J^-1 = det(W(g)) * W(g)^(-T),  J = [[0,1],[-1,0]],
# so tw W theta = J W^vee J^-1 (conjugate to the contragredient). Exact check on random integer matrices.
import sympy as sp, random
random.seed(1)
J = sp.Matrix([[0, 1], [-1, 0]])
ok = True
for _ in range(200):
    while True:
        A = sp.Matrix(2, 2, [random.randint(-9, 9) for _ in range(4)])
        if A.det() != 0: break
    lhs = J*A*J.inv(); rhs = A.det()*A.inv().T
    ok &= (sp.simplify(lhs - rhs) == sp.zeros(2, 2))
print("J W J^-1 = det(W) W^(-T) on 200 random integer matrices:", ok)
print("Hence with det W = theta, theta^2 = 1: W (x) theta = W (x) det(W) ~ J^-1 (W^vee) J ... i.e. tw W theta is isomorphic to W^vee (theta = theta^-1).")
# Used.gen and Used.l_odd: excluded primes for small e (Kato's rational vanishing needs no such exclusion)
for e in range(1, 9):
    ex = [l for l in sp.primerange(2, 60) if l == 2 or (2*e) % (l-1) == 0]
    print("e = %d: primes l < 60 excluded by l_odd / gen ((l-1) | 2e): %s" % (e, ex))
# HasMuQ l L => (l-1) | [L:Q] (l odd) or 2 | [L:Q] (l = 2): LemmaK only concerns l <= n+1
for n in [2, 4, 6, 8, 10, 12, 16, 18, 22]:
    print("n = %2d: primes l with HasMuQ possible: %s" % (n, [l for l in sp.primerange(2, n+2) if (l == 2 and n % 2 == 0) or (l > 2 and n % (l-1) == 0)]))
