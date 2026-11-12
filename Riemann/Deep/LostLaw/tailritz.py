# tailritz.py (lostlaw): rigorous upper bound for lambda_min(Q_c) (even sector) without an N2 x N2 solve: Rayleigh-Ritz on the
# orthonormal set {x/|x|} u {phi_j : N1 < j <= N2}, x = Galerkin ground state at N1 (x lies in span{phi_0..phi_N1}, so x is
# orthogonal to the tail). Entries of Q's Gram matrix from galerkin_fast's closed forms (same formulas, computed per entry, Arb).
# The final Rayleigh quotient of the Ritz vector is evaluated in Arb balls: an upper bound for lambda_min. Usage: tailritz.py c N1 N2 bits
import sys, os, json, time, random, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, acb, arb_mat, acb_mat, ctx
import galerkin_fast as gf
from galerkin_x import vm
from gal import dot, midv, col, tomat, matvec, gs, ritz

def arrays(c, N):
    L = arb(c).log(); pi = arb.pi(); cc = arb(c)
    om = [2*pi*j/L for j in range(N+1)]
    Nj = [(1/L).sqrt()] + [(2/L).sqrt()]*N
    tl = []; n = 0
    while True:
        an = arb(2*n) + arb(1)/2; t = cc**(-an); tl.append((an, t))
        if (t < arb(2)**(-(ctx.prec+40))) and n > 2: break
        n += 1
    tau0 = sum((2*t/an for an, t in tl), arb(0))
    hp0 = acb(arb(1)/4, 0).digamma().real - pi.log()
    base = hp0 + tau0
    S = [arb(0)]*(N+1); Cm = [arb(0)]*(N+1); Um = [arb(0)]*(N+1)
    for j in range(N+1):
        w = om[j]; z = acb(arb(1)/4, w/2)
        ps = z.digamma(); p1 = z.polygamma(1)
        T = -2*sum((t*an/(an*an + w*w) for an, t in tl), arb(0))
        Cm[j] = (ps.real - pi.log()) - hp0 - tau0 - T
        S[j] = ps.imag - 2*w*sum((t/(an*an + w*w) for an, t in tl), arb(0))
        um = arb(0)
        for an, t in tl:
            s = acb(an, -w); um += t*(L/s + 1/(s*s)).real
        Um[j] = p1.real/2 - 2*um
    P = [arb(0)]*(N+1); Dg = [arb(0)]*(N+1)
    n = 2
    while n < c:
        lam = vm(n)
        if lam is not None:
            t = arb(n).log(); wn = lam/arb(n).sqrt()
            for j in range(N+1):
                sj, cj = (om[j]*t).sin(), (om[j]*t).cos()
                P[j] += wn*sj
                Dg[j] += 2*wn*((1 - t/L) if j == 0 else ((L - t)*cj - sj/om[j])/L)
        n += 1
    sh = (L/4).sinh()
    a = [Nj[j]*(1 if j % 2 == 0 else -1)*sh/(om[j]**2 + arb(1)/4) for j in range(N+1)]
    return dict(L=L, pi=pi, om=om, Nj=Nj, base=base, S=S, Cm=Cm, Um=Um, P=P, Dg=Dg, a=a)

def entry(A, j, k):
    L, pi, om, Nj, S, P = A['L'], A['pi'], A['om'], A['Nj'], A['S'], A['P']
    if j == k:
        v = A['base'] + A['Um'][0]/L if j == 0 else A['base'] + A['Cm'][j] + A['Um'][j]/L + S[j]/(L*om[j])
        v = v - A['Dg'][j]
    else:
        sg = 1 if (j+k) % 2 == 0 else -1
        v = -(Nj[j]*Nj[k]/2)*sg*(S[k]-S[j])/(om[j]-om[k]) + (Nj[j]*Nj[k]/2)*sg*(S[k]+S[j])/(2*pi*(j+k)/L)
        v = v - Nj[j]*Nj[k]*sg*((P[k]-P[j])/(om[j]-om[k]) - (P[k]+P[j])/(2*pi*(j+k)/L))
    return v + 2*A['a'][j]*A['a'][k]

def main():
    c = int(sys.argv[1]); N1 = int(sys.argv[2]); N2 = int(sys.argv[3]); ctx.prec = int(sys.argv[4])
    t0 = time.time()
    M = gf.build(c, N1)
    random.seed(1)
    X0 = [[arb(random.uniform(-1, 1)) for _ in range(N1+1)] for _ in range(2)]
    Y = M.solve(tomat(X0), algorithm="approx"); X1 = gs([midv(col(Y, j)) for j in range(2)])
    Y = M.solve(tomat(X1), algorithm="approx"); X2 = gs([midv(col(Y, j)) for j in range(2)])
    Y = M.solve(tomat(X2), algorithm="approx"); X3 = gs([midv(col(Y, j)) for j in range(2)])
    del Y
    E1, x, _ = ritz(M, X3, [], [])
    x = [xi.mid() for xi in x]
    xx = dot(x, x)
    E1r = dot(x, matvec(M, x))/xx          # rigorous Rayleigh quotient at N1
    del M
    A = arrays(c, N2)
    tail = list(range(N1+1, N2+1)); m = len(tail)
    r = [sum((entry(A, j, k)*x[k] for k in range(N1+1)), arb(0)) for j in tail]   # <phi_j, A x>
    T = arb_mat(m, m)
    for p in range(m):
        for q in range(p, m):
            v = entry(A, tail[p], tail[q]); T[p, q] = v; T[q, p] = v
    # Ritz on orthonormal basis {x/|x|, phi_tail}: H = [[x^T M x/|x|^2, r/|x|], [r/|x|, T]]
    nx = xx.sqrt()
    H = arb_mat(m+1, m+1)
    H[0, 0] = E1r
    for p in range(m):
        H[0, p+1] = r[p]/nx; H[p+1, 0] = r[p]/nx
        for q in range(m): H[p+1, q+1] = T[p, q]
    Hm = arb_mat(m+1, m+1)
    for i in range(m+1):
        for j in range(m+1): Hm[i, j] = H[i, j].mid()
    # inverse iteration at shift 0 (H > 0, lambda_1/lambda_2 ~ 1e-7): no full eigendecomposition
    v = tomat([[arb(1)]*(m+1)])
    for _ in range(4):
        v = Hm.solve(v, algorithm="approx")
        vl = [v[i, 0].mid() for i in range(m+1)]; nv = dot(vl, vl).sqrt().mid()
        v = tomat([[(t/nv).mid() for t in vl]])
    w = [v[i, 0].mid() for i in range(m+1)]
    Hv = matvec(H, w)
    X = dot(w, Hv)/dot(w, w)
    print(json.dumps({'c': c, 'N1': N1, 'N2': N2, 'prec': ctx.prec, 'E1_N1_rig_upper': E1r.upper().str(10, radius=False),
                      'tail_rig_upper': X.upper().str(10, radius=False),
                      'gain_nats': float((E1r/X).log().mid()), 'wall': round(time.time() - t0, 1),
                      'rss': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}), flush=True)

if __name__ == '__main__':
    main()
