# float: the paper's own choices (K=350, c=1.82, mu0=1.40625 for every block, last block full) through the per-block
# formula of certify.py, at t = 3.8e7 -> must be <= the paper's region-1 bound 195.52 (uniform relaxations only add).
import math
t = 3.8e7; K = math.ceil(t ** 0.33559); n1 = math.floor(math.sqrt(t / (2 * math.pi))); R = n1 // K; c = 1.82
mu0 = 0.5 * 1.25 * (2.25 + 1 / (3.8e7 ** 0.33559 - 1))
def blk(r, L, mu, c):
    W = math.pi * (r + 1) ** 3 * K ** 3 / t; W3 = W ** (1 / 3)
    al = 1 / c + 32 * mu / (15 * math.sqrt(math.pi)) * math.sqrt(c + 1 / W3) + mu * c / (3 * W3) + mu / (3 * W3 ** 2)
    be = 32 / (3 * math.sqrt(math.pi * c)) + (6 / math.pi - 1) / (c * W3 ** 2) + (1 - 6 / math.pi) / W3
    Q = (L + c * W3) * (L / (c * W3) + (K / W3) * (al - 1 / c) + W3 * be)
    return 2 * math.sqrt(Q / (K * r))
T1 = 2 * sum(1 / math.sqrt(n) for n in range(1, 4 * K))
rem = 1.463 * t ** -0.25 + 0.127 * t ** -0.75
paper_style = T1 + sum(blk(r, K, mu0, c) for r in range(4, R + 1)) + rem
mur = lambda r: 0.5 * (1 + 1 / r) * (2 + 1 / r + 1 / (K - 1))
refined = T1 + sum(blk(r, K, mur(r), c) for r in range(4, R)) + min(2 * (n1 - K * R + 1) / math.sqrt(K * R), blk(R, n1 - K * R + 1, mur(R), c)) + rem
phi = 0.33559; L = math.log(t)
paper = 0.622447 * t ** (1/6) * L - 8.21438 * t ** (1/6) + 8.01143 * t ** (phi / 2) - 2.87533
print("t=3.8e7 K=%d n1=%d R=%d last block %d terms" % (K, n1, R, n1 - K * R + 1))
print("paper bound (CorputSum1, region 1): %.3f  ratio %.5f" % (paper, paper / (t ** (1/6) * L)))
print("per-block formula, paper's mu0, c, full last block: %.3f  ratio %.5f" % (paper_style, paper_style / (t ** (1/6) * L)))
print("per-block mu_r + short last block (same K, c): %.3f  ratio %.5f" % (refined, refined / (t ** (1/6) * L)))
t3only = T1 + sum(blk(r, K, mu0, c) for r in range(4, R)) + min(2 * (n1 - K * R + 1) / math.sqrt(K * R), blk(R, n1 - K * R + 1, mu0, c)) + rem
print("short last block only (paper's mu0, c): %.3f  ratio %.5f" % (t3only, t3only / (t ** (1/6) * L)))
