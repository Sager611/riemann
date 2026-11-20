# ceilaudit6: pick pockets from a pk_dump json: the 4 thinnest, 2 median, 1 largest, and 2 pocket-free window centres.
import sys, json, math
pk = json.load(open(sys.argv[1])); pk = [p for p in pk if p[1] > 0]
pk.sort(key=lambda p: p[1])
n = len(pk); pick = pk[:4] + [pk[n // 2], pk[(3 * n) // 4]] + [pk[-1]]
out = [f'{p[0]}:{p[1]}' for p in pick]
# pocket-free centres: midpoints of the two largest gaps
c = sorted(p[0] for p in pk); gaps = sorted(((c[i + 1] - c[i], 0.5 * (c[i] + c[i + 1])) for i in range(len(c) - 1)), reverse=True)[:2]
out += [f'{g[1]:.2f}:0' for g in gaps]
print(f'# {n} pockets; areas min {pk[0][1]} median {pk[n//2][1]} max {pk[-1][1]}; sum A/2pi {sum(p[1] for p in pk)/(2*math.pi):.3f}', file=sys.stderr)
print(' '.join(out))
