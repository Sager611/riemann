# ceilaudit6: list pockets of a pk_dump json inside [tc - tau1, tc + tau1]; print count, areas, and the semiclassical sum.
import sys, json, math
pk = json.load(open(sys.argv[1])); tc = float(sys.argv[2]); t1 = float(sys.argv[3])
ins = [p for p in pk if abs(p[0] - tc) <= t1]
print(f'pockets in window: {len(ins)}; A/2pi: ' + ', '.join(f'{p[1]/(2*math.pi):.2e}@{p[0]:.0f}' for p in ins) + f'; sum A/2pi {sum(p[1] for p in ins)/(2*math.pi):.3f}')
