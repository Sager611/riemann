# jerby_data.py -- jerbyconf: checks of Jerby's 250,000-point Zenodo table (GramPoints_CORE_1_250000.csv) against 2610.21952:
# class counts (l. 341-343), min mu (l. 353), confined-pair normalised viscosity s*mu = 2 pi mu / log(g/2pi) (l. 657-663) per range,
# trunk-slope ratio rho = |Z'(g_n)| (g_{n+1}-g_n) / A_act, A_act = max(|Z(g_{n-1})|, |Z(g_{n+1})|) (l. 1828-1851).
# data: curl -s -L -o /tmp/jc_zen.zip https://zenodo.org/api/records/23390338/files/Zonedo.zip/content; unzip into /tmp/jc_zen
import csv, math, sys, collections
fn = '/tmp/jc_zen/Zonedo/250000-GramPoints_Dataset/GramPoints_CORE_1_250000.csv'
def num(s):
    s = s.strip()
    if '\x60' in s: s = s.split('\x60')[0]
    if s.endswith('*I'): s = s.split('+')[0]
    return float(s)
G, Zv, Zp, MU, CL = {}, {}, {}, {}, {}
with open(fn) as f:
    r = csv.reader(f); h = next(r)
    for row in r:
        n = int(row[0]); G[n] = float(row[1]); Zv[n] = num(row[2]); Zp[n] = num(row[3]); MU[n] = float(row[4]); CL[n] = row[8]
N = max(G); print('rows', len(G), 'n range', min(G), N, 'labels', collections.Counter(CL.values()))
good = {n: (Zv[n] > 0) == (n % 2 == 0) for n in G}
iso = [n for n in G if not good[n] and good.get(n-1, None) is True and good.get(n+1, None) is True]
print('own count: good', sum(good.values()), 'isolated bad (both neighbours in table)', len(iso), ' GBG label but not own-isolated:',
      [n for n in G if CL[n] == 'GBG' and n not in set(iso)][:10], ' own-isolated but not GBG:', [n for n in iso if CL[n] != 'GBG'][:10])
print('boundary: n=1 good', good[1], ' n=250000 good', good[N], 'class', CL[N], ' n=249999', CL[N-1])
mism = [n for n in G if abs(abs(Zp[n]/Zv[n]) - MU[n]) > 1e-9*max(1, MU[n])]
print('mu_exact vs |Zp/Z| mismatches:', len(mism))
gbg = [n for n in G if CL[n] == 'GBG']
mn = min(gbg, key=lambda n: MU[n]); print('min mu over GBG', MU[mn], 'at n', mn)
s = lambda n: 2*math.pi/math.log(G[n]/(2*math.pi))
sm = {n: s(n)*MU[n] for n in gbg}
low = sorted(gbg, key=lambda n: sm[n])[:12]
print('lowest s*mu over GBG (all confined: no Rosser violation below n = 13,999,525):')
for n in low: print(f'  n={n} g={G[n]:.4f} mu={MU[n]:.6f} s*mu={sm[n]:.5f} Z[-1,0,+1]={Zv[n-1]:.4f},{Zv[n]:.4f},{Zv[n+1]:.4f}')
print('| n range | isolated bad | min s*mu | at n | min mu | # s*mu<2 | # s*mu<2.2 |')
print('|---|---|---|---|---|---|---|')
for a, b in [(1, 1000), (1000, 10000), (10000, 50000), (50000, 100000), (100000, 150000), (150000, 200000), (200000, 250001)]:
    sel = [n for n in gbg if a <= n < b]
    if not sel: continue
    m = min(sel, key=lambda n: sm[n]); mm = min(sel, key=lambda n: MU[n])
    print(f'| [{a}, {b}) | {len(sel)} | {sm[m]:.4f} | {m} | {MU[mm]:.4f} | {sum(1 for n in sel if sm[n] < 2)} | {sum(1 for n in sel if sm[n] < 2.2)} |')
def rho(n):
    A = max(abs(Zv[n-1]), abs(Zv[n+1])); return abs(Zp[n])*(G[n+1]-G[n])/A, A
print('trunk-slope ratio rho = |Z\'(g_n)|(g_{n+1}-g_n)/A_act:')
for lab, L in [('five lowest-viscosity (paper Fig.)', [8569, 45924, 64256, 39324, 96237]), ('five generic (paper Fig.)', [60677, 169832, 34454, 124634, 87302])]:
    print(' ', lab, ', '.join(f'n={n}: mu={MU[n]:.4f} rho={rho(n)[0]:.3f} A_act={rho(n)[1]:.3f} |Z|={abs(Zv[n]):.3f}' for n in L))
lo40 = sorted(gbg, key=lambda n: MU[n])[:40]
rr = sorted(rho(n)[0] for n in lo40); print(f'  40 lowest-mu isolated bad: rho min {rr[0]:.3f} median {rr[20]:.3f} max {rr[-1]:.3f}')
rall = sorted(rho(n)[0] for n in gbg if n + 1 in G); q = lambda p: rall[int(p*(len(rall)-1))]
print(f'  all {len(gbg)} isolated bad: rho min {rall[0]:.4f} 1% {q(0.01):.3f} median {q(0.5):.3f} 99% {q(0.99):.3f} max {rall[-1]:.1f}')
