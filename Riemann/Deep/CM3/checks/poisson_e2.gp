\\ poisson_e2.gp -- CM2 Lemma 9 (Poisson in t over O_F, local densities of Lemma 10) for mu with SQUARE factors
\\ (e = v_p(n2 n2') = 2: the pairs with common factors that Lemma P keeps), F = Q(sqrt5), l = O_F, c = 1, S0 empty.
\\ brute:  2 * sum_{t in O_F, t^2-4mu << 0} sum_{y | f(t)} Phi(N(t^2-4mu)/(N(y)^2 X')),  (t^2-4mu) = f(t)^2 d_{K_t/F}  (rnfdisc)
\\ main :  (2^{1-n}/sqrt(d_F)) (X'/sqrt(M0)) sum_y A_mu(y) K_2(N(y)^2 v),  v = X'/(4^n M0),
\\         A_mu(p^k) = #{x mod p^{2k} : C_k(x, mu)} by brute force (C_k of CM2 Lemma 1(b), dyadic prime included).
default(parisizemax, 256*10^6);
alarm(1500);
\p 30
F = bnfinit(y^2 - y - 1, 1); rho = real(polroots(F.pol)); rho = [vecmax(rho), vecmin(rho)];
emb(a, v) = subst(lift(a), y, rho[v]);
sm(t) = { if (t <= 0, return(0)); if (t >= 1, return(1)); my(z = 1/t - 1/(1-t)); if (z > 2000, return(0)); if (z < -2000, return(1)); 1/(1 + exp(z)) };
Phi(u) = if (u <= 1/4 || u >= 2, 0, if (u < 1/2, sm(4*u - 1), if (u <= 1, 1, sm(2 - u))));
\\ g(lam) = int_{-1}^{1} Phi(lam (1 - s^2)) ds ;  K2(v) = v^{-1} int_{-1}^{1} g((1-tau^2)/v) dtau
g(lam) = if (lam <= 1/4, 0, 2*intnum(s = 0, 1, Phi(lam*(1 - s^2))));
gs(lam) = { \\ = g(lam) = int_0^lam Phi(w)/sqrt(lam(lam-w)) dw; with w = lam - r^2: (2/sqrt(lam)) int Phi(lam - r^2) dr
  if (lam <= 1/4, return(0));
  my(r0 = sqrt(max(0, lam - 2)), r1 = sqrt(lam - 1/4));
  (2/sqrt(lam)) * intnum(r = r0, r1, Phi(lam - r^2))
};
K2(v) = { if (v >= 4, return(0));
  \\ tau in [0,1): (1 - tau^2)/v ranges over (0, 1/v]; need > 1/4 i.e. tau^2 < 1 - v/4
  my(tm = sqrt(max(0, 1 - v/4)));
  2*intnum(tau = 0, tm, gs((1 - tau^2)/v)) / v
};
\\ totally positive generator of a prime above p (degree one)
tpgen(pr) = { my(g = nfbasistoalg(F, bnfisprincipal(F, pr, 1)[2]), u = [1, -1, y, -y, y^2, -y^2, y^-1, -y^-1]);
  for (i = 1, #u, my(h = Mod(g*u[i], F.pol)); if (emb(h,1) > 0 && emb(h,2) > 0, return(h))); error("tp") };
nrm(a) = norm(Mod(a, F.pol));
\\ reps of O_F / p^m: a*1 + b*y with a, b mod (p-part): use HNF of p^m
reps(I) = { my(H = idealhnf(F, I), L = List()); for (a = 0, H[1,1]-1, for (b = 0, H[2,2]-1, listput(L, a + b*y))); Vec(L) };
inid(z, I) = { my(H = idealhnf(F, I)); denominator(H^(-1) * nfalgtobasis(F, z)) == 1 };
\\ A_{mu,p}(k): # x mod p^{2k} with C_k(x, mu): exists a: p^k | x - 2a, p^{2k} | a^2 - x a + mu
ACACHE = Map();
Aloc(pr, k, mu) = { if (k == 0, return(1));
  my(key = [idealhnf(F, pr), k, mu], val); if (mapisdefined(ACACHE, key, &val), return(val));
  my(P2 = idealpow(F, pr, 2*k), Pk = idealpow(F, pr, k), X = reps(P2), cnt = 0);
  for (i = 1, #X, my(x = X[i], ok = 0);
    for (j = 1, #X, my(a = X[j]); if (inid(x - 2*a, Pk) && inid(a^2 - x*a + mu, P2), ok = 1; break));
    cnt += ok);
  mapput(ACACHE, key, cnt); cnt };
runtest(name, mu, xlist) = {
  my(M0 = nrm(mu), B1 = 2*sqrt(emb(mu,1)), B2 = 2*sqrt(emb(mu,2)), tl = List());
  \\ t = a + b y : sigma_v(t) = a + b rho_v
  my(bmax = floor((B1 + B2)/(rho[1] - rho[2])) + 1);
  for (b = -bmax, bmax, my(lo = ceil(max(-B1 - b*rho[1], -B2 - b*rho[2])), hi = floor(min(B1 - b*rho[1], B2 - b*rho[2])));
    for (a = lo, hi, my(t = Mod(a + b*y, F.pol), D = t^2 - 4*mu);
      if (emb(D,1) < 0 && emb(D,2) < 0,
        my(dd = rnfdisc(F, x^2 - lift(t)*x + lift(mu))[1], ff = idealfactor(F, idealdiv(F, D, dd)), nys = [1]);
        for (i = 1, #ff~, my(e = ff[i,2]); if (e % 2, error("odd")); my(np = idealnorm(F, ff[i,1]), new = List());
          for (s = 1, #nys, for (j = 0, e/2, listput(new, nys[s]*np^j))); nys = Vec(new));
        listput(tl, [abs(nrm(D)), nys]))));
  tl = Vec(tl);
  printf("%s: M0 = %d = %s, #t = %d\n", name, M0, factor(M0)[,1]~, #tl);
  \\ local densities for all y with N(y) <= Ymax
  for (ix = 1, #xlist, my(Xp = xlist[ix], v = Xp/(16*M0), brute = 0., main = 0.);
    for (i = 1, #tl, for (j = 1, #tl[i][2], brute += Phi(tl[i][1]/(tl[i][2][j]^2*Xp))));
    brute *= 2;
    my(Ymax = floor(2/sqrt(v)) + 1, Lid = ideallist(F, Ymax), nterm = 0);
    for (Ny = 1, Ymax, for (q = 1, #Lid[Ny], my(yy = Lid[Ny][q], fa = idealfactor(F, yy), A = 1);
      if (Ny^2*v >= 4, next);
      for (i = 1, #fa~, A *= Aloc(fa[i,1], fa[i,2], mu));
      if (A, nterm++; main += A * K2(Ny^2*v))));
    main *= (1/(2*sqrt(5))) * Xp/sqrt(M0);
    printf("   X'/sqrt(M0) = %.1f: brute %.4f  main %.4f  ratio %.4f  (%d y-terms)\n", Xp/sqrt(M0), brute, main, brute/main, nterm));
};
{
printf("K2 Mellin check at z=1/2: %.6f vs %.6f\n", intnum(v = 0, 4, K2(v)*v^(-1/2)), Pi^2*intnum(u = 1/4, 2, Phi(u)*u^(-1/2)));
P11 = idealprimedec(F, 11); P19 = idealprimedec(F, 19); P29 = idealprimedec(F, 29); P31 = idealprimedec(F, 31);
p11 = tpgen(P11[1]); p19 = tpgen(P19[1]); p29 = tpgen(P29[1]); p31 = tpgen(P31[1]); p19b = tpgen(P19[2]);
mu1 = p11^2 * p19 * p29;           \\ e = 2 at p11
mu2 = p11^2 * p19^2 * p29;         \\ e = 2 at p11 and p19
mu3 = p11 * p19 * p29 * p31;       \\ squarefree control
mu4 = p11^2 * p19 * p19b;          \\ e = 2 at p11; the two primes above 19
runtest("mu1 = p11^2 p19 p29", mu1, [20, 40] * sqrt(nrm(mu1)));
runtest("mu2 = p11^2 p19^2 p29", mu2, [20, 40] * sqrt(nrm(mu2)));
runtest("mu3 = p11 p19 p29 p31 (control)", mu3, [20, 40] * sqrt(nrm(mu3)));
runtest("mu4 = p11^2 p19 p19'", mu4, [20, 40, 80, 160] * sqrt(nrm(mu4)));
runtest("mu1 again, larger X'", mu1, [80, 160] * sqrt(nrm(mu1)));
}

