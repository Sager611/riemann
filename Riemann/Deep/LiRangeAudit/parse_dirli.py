
import re
from fractions import Fraction as F
def t0(q):
    return max(F(10**8, q), F(75*10**6, q) + 200) if q % 2 == 0 else max(F(10**8, q), F(375*10**5, q) + 200)
rows = {}
for L in open('Riemann/Deep/DirLi/results_q1000.txt'):
    m = re.match(r'q (\d+) t0 (\S+) N_half (\d+) r_half (\S+) N_38 (\d+) r_38 (\S+)', L)
    if m: rows[int(m[1])] = (F(m[2]), int(m[3]), int(m[5]))
need = [q for q in range(3, 1001) if q % 4 != 2]
print("q1000: rows", len(rows), "needed", len(need), "missing", [q for q in need if q not in rows][:5], "t0 mismatches", sum(1 for q in rows if rows[q][0] != t0(q)))
rh = {q: rows[q][1] / t0(q) ** 2 for q in rows}; r3 = {q: rows[q][2] / t0(q) ** 2 for q in rows}
qs = sorted(rows)
print("  min r_half %.5f at q=%d, min r_38 %.5f at q=%d" % (float(min(rh.values())), min(rh, key=rh.get), float(min(r3.values())), min(r3, key=r3.get)))
nonmono_N = [(qs[i], qs[i+1]) for i in range(len(qs)-1) if rows[qs[i+1]][1] > rows[qs[i]][1]]
nonmono_r = [(qs[i], qs[i+1]) for i in range(len(qs)-1) if rh[qs[i+1]] > rh[qs[i]]]
print("  N_half increases between consecutive q: %d cases %s; r_half increases: %d cases %s" % (len(nonmono_N), nonmono_N[:4], len(nonmono_r), nonmono_r[:4]))
print("  every q<=1000: N_half >= 41.30e16/q^2:", all(rows[q][1] >= F(4130, 100) * 10**16 / q**2 for q in qs), "; N_38 >= 55.94e16/q^2:", all(rows[q][2] >= F(5594, 100) * 10**16 / q**2 for q in qs), "; q=3:", rows[3][1], rows[3][2])
bl = []
for L in open('Riemann/Deep/DirLi/results_blocks.txt'):
    m = re.match(r'block (\w+) (\d+) (\d+) t0min (\S+) N_half (\d+) N_38 (\d+)', L)
    if m: bl.append((m[1], int(m[2]), int(m[3]), F(m[4]), int(m[5]), int(m[6])))
cov = set()
for k, q1, q2, t, nh, n3 in bl:
    assert t == t0(q2), (k, q1, q2)
    assert q2 <= 1.02 * q1 + 1e-9
    for q in range(q1, q2 + 1):
        if (k == 'odd' and q % 2 == 1) or (k == 'even' and q % 4 == 0): cov.add(q)
needb = [q for q in range(1001, 400001) if q % 4 != 2]
print("blocks:", len(bl), "covered", len(cov), "needed", len(needb), "missing", len(set(needb) - cov), "extra", len(cov - set(needb)))
ch = min(bl, key=lambda b: F(b[4]) / t0(b[1]) ** 2); c3 = min(bl, key=lambda b: F(b[5]) / t0(b[1]) ** 2)
print("  min N_half/t0(q1)^2 = %.5f block %s; min N_38/t0(q1)^2 = %.5f block %s" % (float(F(ch[4]) / t0(ch[1]) ** 2), ch[:3], float(F(c3[5]) / t0(c3[1]) ** 2), c3[:3]))
# uniform claim: for each q in block, N_block >= 26 t0(q)^2 and >= 35 t0(q)^2 (t0(q) <= t0(q1))
okh = all(F(b[4]) >= F(26) * t0(q) ** 2 for b in bl for q in (b[1], b[2])); ok3 = all(F(b[5]) >= F(35) * t0(q) ** 2 for b in bl for q in (b[1], b[2]))
print("  N_block >= 26 t0(q)^2 and >= 35 t0(q)^2 at both ends of every block:", okh, ok3, "; t0 nonincreasing in q within parity (q<=4e5):", all(t0(q + 2) <= t0(q) for q in range(3, 399999, 2)) and all(t0(q + 4) <= t0(q) for q in range(4, 399997, 4)))

