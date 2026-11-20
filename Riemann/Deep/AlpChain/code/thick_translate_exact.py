# Exact check of ThickPartTranslateMean (Formalism/SizeWindowBrackets/GrowthBoundGivesRiemannHypothesis.lean:44-50)
# and of its subgroup version used by ch_core (average over u0 in A0, T2 inside A0): the identity is linear in the
# values GPsi(psi, A/b), so we compare the coefficient of each character psi exactly (Gaussian rationals).
import random
from fractions import Fraction as Fr
def cm(x, y): return (x[0]*y[0]-x[1]*y[1], x[0]*y[1]+x[1]*y[0])
def ca(x, y): return (x[0]+y[0], x[1]+y[1])
def cj(x): return (x[0], -x[1])
ZERO, ONE = (Fr(0), Fr(0)), (Fr(1), Fr(0))
def run(nH, seed, sub=None):
    rnd = random.Random(seed)
    zeta = {1: [ONE], 2: [ONE, (Fr(-1), Fr(0))], 4: [ONE, (Fr(0), Fr(1)), (Fr(-1), Fr(0)), (Fr(0), Fr(-1))]}[nH]
    G = [((e1, e2), h) for e1 in (0, 1) for e2 in (0, 1) for h in range(nH)]
    chars = [(a, b_, h) for a in (0, 1) for b_ in (0, 1) for h in range(nH)]
    psi = lambda c, g: cm(((-1) ** ((c[0]*g[0][0] + c[1]*g[0][1]) % 2) * Fr(1), Fr(0)), zeta[(c[2]*g[1]) % nH])
    inPhi = lambda g: g[0][0] == 0
    nPhi = sum(1 for g in G if inPhi(g))
    partner = lambda t, g: (((g[0][0] + 1) % 2, (g[0][1] + 1) % 2), (g[1] + t) % nH)
    thickPl = lambda t: ((0, 1), t % nH)
    A0 = list(range(nH)) if sub is None else sub
    T2 = sorted(set(rnd.choice(A0) for _ in range(rnd.randint(1, len(A0)))))
    t = rnd.randrange(nH); b = Fr(rnd.randint(2, 9), rnd.randint(1, 3))
    coef = {}
    for c in chars:
        s = ZERO
        for u0 in A0:
            thick = {thickPl(tp + u0) for tp in T2}
            for g in G:
                if not inPhi(g) and partner(t, g) in thick:
                    s = ca(s, cm(cj(psi(c, g)), (b, Fr(0))))
        coef[c] = (s[0] / (len(A0) * nPhi), s[1] / (len(A0) * nPhi))
    d2 = Fr(len(T2), nH)
    if sub is None:   # full average: RHS = (delta2 b / 2) * (G0 + G_psC - G_chC - G_chpsC)
        sign = {(0, 0, 0): 1, (0, 1, 0): 1, (1, 0, 0): -1, (1, 1, 0): -1}
        ok = all(coef[c] == (d2 * b / 2 * sign.get(c, 0), Fr(0)) for c in chars)
        return ok, len(T2), t
    # subgroup average: coefficient = (|T2| b /(2n)) * conj psi((1,0),0) * psi((0,0),t) * [psi trivial on A0]
    ok = True
    for c in chars:
        triv = all(psi(c, ((0, 0), s)) == ONE for s in A0)
        pred = ZERO if not triv else cm(cm((d2 * b / 2, Fr(0)), cj(psi(c, ((1, 0), 0)))), psi(c, ((0, 0), t)))
        ok &= coef[c] == pred
    return ok, len(T2), t
tot = 0
for nH in (1, 2, 4):
    for seed in range(5):
        ok, k, t = run(nH, seed); tot += ok
        print(f"full average H=Z/{nH} seed={seed} |T2|={k} t={t}: identity exact: {ok}")
for seed in range(5):
    ok, k, t = run(4, seed, sub=[0, 2]); tot += ok
    print(f"subgroup A0={{0,2}} in Z/4 seed={seed} |T2|={k} t={t}: coefficient = (|T2|b/2n) conj psi(sgc0) psi(0,t) [psi|A0=1]: {ok}")
print(f"PASS {tot}/20" if tot == 20 else f"FAIL {tot}/20")
