\\ serrenumK sanity checks (PARI/GP). Run: gp -q sanity.gp > logs/sanity_output.txt
\\ Theorem A/B mechanics over K = Q(sqrt 5), Q(sqrt 2) with curves defined over K (base changes and twists by non-rational d).
default(parisize, 64000000);
radlog(nf, I) = my(F = idealfactor(nf, I)[,1]); sum(i = 1, #F, log(F[i].p^F[i].f));
chi(nf, D, pr) = my(m = nfmodprinit(nf, pr), r = nfmodpr(nf, D, m)); if (r == 0, 0, if (issquare(r), 1, -1));
deg1(nf, pmax) = my(L = List()); forprime (p = 3, pmax, foreach (idealprimedec(nf, p), pr, if (pr.f == 1, listput(L, pr)))); Vec(L);
twist(nf, E, d) = my(c4 = E.c4, c6 = E.c6); ellinit([0, 0, 0, -27*c4*d^2, -54*c6*d^3], nf);
{
fields = [[y^2 - 5, "Q(sqrt5)", 124.201, 651.302, (1 + y)/2], [y^2 - 2, "Q(sqrt2)", 124.201, 680.540, 1 + y]];
base = [[0,-1,1,-10,-20], [1,0,1,4,-6], [1,1,1,-10,-10], [1,-1,1,-1,-14], [0,1,1,-9,-15], [0,0,1,-1,0], [1,-1,1,-3,3], [0,1,1,-23,-50]];
foreach (fields, F,
  my(nf = nfinit(F[1]), a = F[3], b = F[4], u = F[5], P = deg1(nf, 3000), Es = List(), worst = [0, 0], npairs = 0, ndiff = 0);
  foreach (base, c, my(E = ellinit(c, nf)); listput(Es, E); listput(Es, twist(nf, E, Mod(u, nf.pol))); listput(Es, twist(nf, E, Mod(y, nf.pol))));
  my(N = vector(#Es, i, ellglobalred(Es[i])[1]), V = vector(#Es, i, vector(#P, k, if (idealval(nf, idealmul(nf, 2, N[i]), P[k]) > 0, oo, ellap(Es[i], P[k])))));
  for (i = 1, #Es, for (j = i + 1, #Es,
    my(m = idealmul(nf, idealmul(nf, 2, N[i]), N[j]), k0 = 0);
    for (k = 1, #P, if (V[i][k] != oo && V[j][k] != oo && V[i][k] != V[j][k], k0 = k; break));
    npairs++;
    if (k0, ndiff++; my(p0 = P[k0].p, rho = radlog(nf, m), r = log(p0) - 32*log(a*rho + b)); if (p0 > worst[1], worst = [p0, r]))));
  print(F[2], ": ", #Es, " curves over K (8 base changes, 16 twists by d = ", u, ", ", y, "), ", npairs, " pairs, ", ndiff, " with a separating degree-1 good prime <= 3000;");
  print("   largest least separating prime p0 = ", worst[1], ", log(p0) - 32 log(a rho + b) = ", worst[2]);
  \\ Theorem B mechanics: twists chi_D, D in {u, y, -1, -u} (units / sqrt): least good degree-1 p with chi_D(p) = -1 and a_p != 0
  my(mx = 0);
  for (i = 1, #Es, foreach ([Mod(u, nf.pol), Mod(y, nf.pol), Mod(-1, nf.pol), Mod(-u, nf.pol)], D,
    for (k = 1, #P, if (V[i][k] != oo && V[i][k] != 0 && chi(nf, D, P[k]) == -1, mx = max(mx, P[k].p); break))));
  print("   twist step: least good degree-1 p with chi_D(p) = -1, a_p != 0 is <= ", mx, " for all ", #Es, " curves and 4 characters");
  \\ non-surjective primes from isogenies over K (base changes of 11a1, 17- and 37-isogeny curves)
  foreach ([[[0,-1,1,-10,-20], 5], [[1,0,1,-190891,-36002922], 17], [[1,1,1,-8,6], 37]], c,
    my(E = ellinit(c[1], nf), M = iferr(ellisomat(E)[2], e, 0));
    my(degs = if (M === 0, "ellisomat n/a", Set(concat(Vec(M)))));
    my(rho = radlog(nf, idealmul(nf, 2, ellglobalred(E)[1])));
    print("   ", c[1], ": isogeny degrees over K ", degs, "; known non-surjective l = ", c[2], " <= 2(a rho + b)^16 = 10^", 16*log(a*rho + b)/log(10) + log(2)/log(10)));
  \\ Hensel-Serre (4.4) on L = K(E[2]) for the base curves: log d_L <= m log d_K + (m-1) rho + m n_K log m, m = [L:K]
  my(worstH = -oo, cnt = 0);
  foreach (base, c, my(E = ellinit(c, nf), f = elldivpol(ellinit(c), 2), g = nfsplitting(f), Lp = polredbest(polcompositum(subst(nf.pol, y, x), g)[1]));
    my(m = poldegree(Lp) / 2, rho = radlog(nf, idealmul(nf, 2, ellglobalred(E)[1])), lhs = log(abs(nfdisc(Lp))));
    my(rhs = m*log(abs(nf.disc)) + (m - 1)*rho + m*2*log(m)); worstH = max(worstH, lhs - rhs); cnt++);
  print("   Hensel-Serre (4.4) on ", cnt, " fields K(E[2]): max(log d_L - bound) = ", worstH, " (must be <= 0)");
);
}
