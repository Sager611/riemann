
import sys
exec(open('indep.py').read().split("if __name__ == '__main__':")[0])
for kind, end in (('N', 47.799393), ('R', 46.810759)):
    out = []
    for d in (0.002, 0.005, 0.01, 0.02, 0.04):
        La = arb(end + d)*LN10; Lb = La + arb('1e-6')
        m, k, lam, t = choose(kind, La, Lb, arb(1), kmax=120, fine=True)
        out.append(f'+{d}: {float(m.lower()):.3g} (k={k})')
    print(f'b0=1 {kind}, own fine grid, point boxes (width 1e-6) at 10^({end} + d):', '; '.join(out), flush=True)
