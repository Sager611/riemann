
import json, math, random
random.seed(20261114)
res = json.load(open('Riemann/Deep/Ingham3/logs/cert_result.json'))['main']
cells = []
def scanlist(Ls, Le, w):
    L = Ls; out = []
    while L < Le:
        Lb = min(L + w, Le); out.append((L, Lb)); L = Lb
    return out
s1 = scanlist(7*math.log(100), 180, 0.25); s2 = scanlist(180, 240, 0.02)
print('k=7 scan cells', len(s1), len(s2))
cells += [dict(k=7, La=s1[0][0], Lb=s1[0][1], tag='k7first'), dict(k=7, La=s1[-1][0], Lb=s1[-1][1], tag='k7end1'),
          dict(k=7, La=s2[0][0], Lb=s2[0][1], tag='k7start2'), dict(k=7, La=s2[-1][0], Lb=s2[-1][1], tag='k7last')]
wi = min(range(len(s2)), key=lambda i: abs(s2[i][0] - 205.52))
cells.append(dict(k=7, La=s2[wi][0], Lb=s2[wi][1], tag='k7worst'))
for i in random.sample(range(len(s1)), 3): cells.append(dict(k=7, La=s1[i][0], Lb=s1[i][1], tag='k7rand1'))
for i in random.sample(range(len(s2)), 3): cells.append(dict(k=7, La=s2[i][0], Lb=s2[i][1], tag='k7rand2'))
W = {3: (0.25, 2400), 4: (0.1, 640), 5: (0.05, 400), 6: (0.05, 300)}
for k in (6, 5, 4, 3):
    L0 = res[str(k)]['L0']; w, Ls = W[k]
    sc = scanlist(L0, Ls, w)
    assert len(sc) == res[str(k)]['cells'], (k, len(sc), res[str(k)]['cells'])
    cells.append(dict(k=k, La=sc[0][0], Lb=sc[0][1], tag='k%dfirst' % k))
    i = random.randrange(len(sc)); cells.append(dict(k=k, La=sc[i][0], Lb=sc[i][1], tag='k%drand' % k))
# the k = 5 sliver: decimal 288.73 < float L0 = 288.73000000000001819 -> the cell [float below 288.73, L0]
cells.append(dict(k=5, La=math.nextafter(288.73, 0), Lb=288.73, tag='k5sliver'))
json.dump(cells, open('Riemann/Deep/Ingham3Audit/logs/cells.json', 'w'), indent=1)
print(len(cells), 'cells;', [(c['tag'], round(c['La'], 4)) for c in cells])
