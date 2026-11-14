
import sys, json
sys.argv = ['indep.py']
exec(open('indep.py').read().split("if __name__ == '__main__':")[0])
D = json.load(open('/Users/adriansager/playground/riemann/.scratch/robinchaudit/boxes.json')); B = D['boxes']
def A(me): return arb(me[0])*arb(2)**me[1]
for s in D['summary']:
    if s['X'] != 'probe': continue
    r = B[s['last'] - 1]                         # the last (narrowest) failing box of the probe
    La, Lb = A(r['La']), A(r['Lb']); b0 = arb(r['b0'])
    m, k, lam, t = choose(r['kind'], La, Lb, b0, kmax=120, fine=True)
    print(f"{r['kind']} b0={r['b0']} box at 10^{float((La/LN10).mid()):.6f} (width {float((Lb-La).mid()):.2g}): RobinCH m {r['m']:.3g} (k={r['k']}); "
          f"own fine grid (4 best k, lam x 0.5..2.0 in 31 steps): best m {float(m.lower()):.3g} at k={k}, lam={lam:.4g}")
