#!/usr/bin/env python3
# AlpAdj probes P1-P4 (PARI via cypari2). Run: uv run --no-project --with cypari2 python3 probes.py <P1|P2|P3|P4>
import sys, cypari2
pari = cypari2.Pari(); pari.allocatemem(120 * 10**6, silent=True)
gp = lambda s: pari(s)
which = sys.argv[1]
if which == 'P1':
    # weight-one newforms of small level: dihedral type, real-quadratic RM detection, Satake ratio at split l
    gp(r'''rmD(N, v) = {my(Ds = List()); for(D = 2, 4*N, if(isfundamental(D),
      my(ok = 1, cnt = 0); forprime(p = 3, #v - 1, if(N % p && kronecker(D, p) == -1, cnt++; if(v[p + 1] != 0, ok = 0; break)));
      if(ok && cnt > 5, listput(Ds, D)))); Vec(Ds)}''')
    gp(r'''scan(Nmax) = {my(out = List()); for(N = 1, Nmax, my(G = znstar(N, 1)); my(L = chargalois(G));
      for(i = 1, #L, my(chi = L[i]); if(zncharisodd(G, chi),
        my(mf = mfinit([N, 1, [G, chi]], 0)); if(mfdim(mf) > 0,
          my(B = mfeigenbasis(mf)); for(j = 1, #B, my(t = mfgaloistype(mf, B[j]));
            if(t > 0, my(v = mfcoefs(B[j], 400)); my(v1 = vector(#v, k, if(type(v[k]) == "t_POLMOD", lift(v[k]), v[k])));
              my(vz = vector(#v1, k, v1[k] == 0)); my(Ds = rmD(N, vector(#v1, k, if(vz[k], 0, 1))));
              if(#Ds, listput(out, [N, i, j, t, Ds])))))))); Vec(out)}''')
    res = gp('scan(150)')
    print('dihedral weight-one newforms with a real-quadratic D (a_p = 0 at all p <= 400 inert in Q(sqrt D)), N <= 150:')
    for r in res: print('  N =', r[0], 'char orbit', r[1], 'form', r[2], 'galoistype D_' + str(r[3]), 'real D =', r[4])
    if len(res):
        N = int(res[0][0]); D = int(res[0][4][0])
        print('first example: N =', N, 'D =', D)
        gp('G0 = znstar(%d, 1); L0 = chargalois(G0); chi0 = L0[%d]; mf0 = mfinit([%d, 1, [G0, chi0]], 0); f0 = mfeigenbasis(mf0)[%d]' % (N, int(res[0][1]), N, int(res[0][2])))
        print('  field of coefficients:', gp('mffields(mf0)'))
        k = 0
        for l in gp('primes([3, 200])'):
            l = int(l)
            if N % l and int(gp('kronecker(%d, %d)' % (D, l))) == 1:
                al = gp('mfcoef(f0, %d)' % l); ch = gp('chareval(G0, chi0, %d)' % l)
                # Satake polynomial X^2 - a_l X + chi(l); roots are roots of unity (finite image): ratio has finite order
                print('  l =', l, 'split in Q(sqrt %d): a_l =' % D, al, ' chi(l) = e(2 pi i *', ch, ')')
                k += 1
                if k >= 4: break
elif which == 'P2':
    # elliptic curves with rational l-torsion (residually reducible) and SPLIT multiplicative reduction at l
    found = {3: [], 5: [], 7: []}
    for a1 in range(-12, 13):
        for a3 in range(1, 25):
            if int(gp('E = ellinit([%d, 0, %d, 0, 0]); #E' % (a1, a3))) == 0: continue
            N = gp('ellglobalred(E)[1]')
            if int(gp('valuation(%s, 3)' % N)) == 1 and int(gp('ellap(E, 3)')) == 1 and len(found[3]) < 3:
                found[3].append((a1, a3, int(N), str(gp('E.j')), int(gp('#ellisomat(E)[2]'))))
    for t in range(-30, 31):
        if t in (0,): continue
        s = 'E = ellinit([1 - (%d), -(%d), -(%d), 0, 0])' % (t, t, t)
        if int(gp(s + '; #E')) == 0: continue
        N = gp('ellglobalred(E)[1]')
        if int(gp('valuation(%s, 5)' % N)) == 1 and int(gp('ellap(E, 5)')) == 1 and len(found[5]) < 3:
            found[5].append((t, int(N), str(gp('E.j'))))
    print('l = 3 (Hessian y^2 + a1 xy + a3 y = x^3, 3-torsion point (0,0)), v_3(N) = 1, a_3 = +1 (split):')
    for f in found[3]: print('  [a1, a3] =', f[:2], 'N =', f[2], 'j =', f[3], '#isogeny class =', f[4])
    print('l = 5 (Tate normal form E(t), 5-torsion point), v_5(N) = 1, a_5 = +1 (split):')
    for f in found[5]: print('  t =', f[0], 'N =', f[1], 'j =', f[2])
    if found[3]:
        a1, a3 = found[3][0][:2]
        gp('E = ellinit([%d, 0, %d, 0, 0])' % (a1, a3))
        print('check first l = 3 curve: CM?', gp('ellisCM = if(type(E.j) == "t_INT" && setsearch([0, 1728, -3375, 8000, -32768, 54000, 287496, -884736, -12288000, 16581375, -884736000, -147197952000, -262537412640768000], E.j), 1, 0)'),
              ' v_3(j) =', gp('valuation(E.j, 3)'), '(< 0: multiplicative)')
        print('  Tate curve at 3: rho|G_3 = (chi * ; 0 1) with unramified parts a_1 = a_2 = 1 (a_3 = +1) => psi_1/psi_2 = chi_cyc exactly => FreeRatios fails on G_3 cap G\'')
elif which == 'P3':
    # trace part: Z_l-rank of the maximal abelian pro-l extension unramified outside l (ray class groups mod l^k), CM fields
    for (pol, l, kmax, name) in [('x^2 + 6', 5, 6, 'Q(sqrt -6)'), ('x^2 + 14', 3, 6, 'Q(sqrt -14)'), ('x^4 + 1', 17, 3, 'Q(zeta_8)')]:
        gp('bnf = bnfinit(%s, 1)' % pol)
        h = int(gp('bnf.no')); vl = []
        for k in range(1, kmax + 1):
            c = gp('bnrclassno(bnf, %d^%d)' % (l, k)); vl.append(int(gp('valuation(%s, %d)' % (c, l))))
        diffs = [vl[i + 1] - vl[i] for i in range(len(vl) - 1)]
        print(name, 'l =', l, 'splitting:', gp('#idealprimedec(bnf, %d)' % l), 'primes; h =', h, 'v_l(h) =', int(gp('valuation(%d, %d)' % (h, l))),
              '; v_l #Cl(L, l^k), k = 1..%d:' % kmax, vl, 'increments', diffs)
elif which == 'P4':
    # Greenberg-Wiles count for the fixed-weight ordinary condition on ad~ over L+ (CHT 'odd' polarization)
    for n in (2, 3, 4):
        for d in (1, 2, 3, 4):
            loc = d * n * (n - 1) // 2      # sum_{v|l} [L+_v:Q_l] * dim(b / gr-diagonal unramified) - h0 = d*n(n-1)/2
            inf = d * n * (n - 1) // 2      # sum_{v|inf} h0(G_v, ad~) = d * n(n-1)/2 for the odd sign
            free = d * n                    # free weight adds Hom(O_v^x, Q_l) for each diagonal entry
            print('n=%d d=%d: fixed weight chi = %d ; free weight chi = %d' % (n, d, loc - inf, loc + free - inf))
