\\ Numeric sanity checks for SerreNum (PARI/GP 2.19; not part of the proof).
\\ Run: gp -q sanity.gp > logs/sanity_output.txt
default(parisize, "120M");
LOG10 = log(10);
\\ certified bounds (cert.py): p <= (a log R + b)^32, C_E <= 2 (a log rad(2N) + b)^16
bnd(a, b, R) = (a*log(R) + b)^32;
radi(n) = factorback(factor(n)[,1]);
\\ ---------------------------------------------------------------- 1. curves of small conductor
curves = List(); seen = Map();
{
for (a1 = 0, 1, for (a3 = 0, 1, for (a2 = -1, 1, for (a4 = -12, 12, for (a6 = -12, 12,
  my(E, Em, N);
  iferr(E = ellinit([a1, a2, a3, a4, a6]), err, next);
  if (#E == 0 || E.disc == 0, next);
  Em = ellminimalmodel(E); N = ellglobalred(Em)[1];
  if (N > 300, next);
  my(key = [N, Em.j, Em.a1, Em.a2, Em.a3, Em.a4, Em.a6]);
  if (mapisdefined(seen, key), next);
  mapput(seen, key, 1);
  listput(curves, [N, Em]))))));
}
print("1. curves with conductor <= 300 from small Weierstrass coefficients: ", #curves);
P = primes(60);
apvec(c) = vector(#P, i, if (c[1] % P[i] == 0, oo, ellap(c[2], P[i])));
\\ group curves by their a_p at good p <= 281; curves in different groups are certainly not isogenous
PP = primes(60);
cls = List(); clsrep = List();
{
for (i = 1, #curves,
  my(c = curves[i], v = vector(#PP, k, if (c[1] % PP[k] == 0, 0, ellap(c[2], PP[k]))), found = 0);
  for (j = 1, #cls, if (cls[j][1] == c[1] && cls[j][2] == v, found = 1; break));
  if (!found, listput(cls, [c[1], v]); listput(clsrep, c)));
}
print("   classes with distinct a_p (good p <= 281) among them: ", #cls, " (distinct a_p => not isogenous)");
\\ least good prime with a_p(E) != a_p(E') for every pair of classes; compare with Theorem 1
worst = [0]; worstratio = -oo; npairs = 0; nmod2 = 0; worstmod2 = 0;
{
for (i = 1, #clsrep, for (j = i + 1, #clsrep,
  my(E = clsrep[i][2], F = clsrep[j][2], N1 = clsrep[i][1], N2 = clsrep[j][1], M = 2*N1*N2, R = radi(M), p0 = 0, p2 = 0);
  forprime (p = 3, 10^5, if (M % p == 0, next);
    my(x = ellap(E, p), y = ellap(F, p));
    if (!p2 && (x - y) % 2 != 0, p2 = p);
    if (x != y, p0 = p; break));
  if (!p0, error("no separating prime below 10^5 for ", [N1, N2]));
  npairs++;
  if (p0 > worst[1], worst = [p0, N1, N2]);
  my(r = log(p0)/LOG10 - log(bnd(17, 100, R))/LOG10);
  if (r > worstratio, worstratio = r);
  if (p2, nmod2++; if (p2 > worstmod2, worstmod2 = p2))));
}
print("   non-isogenous pairs: ", npairs, "; largest least separating good prime: ", worst[1], " (conductors ", worst[2], ", ", worst[3], ")");
print("   max over pairs of log10(p / (17 log R + 100)^32) = ", precision(worstratio, 10), "  (must be < 0)");
print("   pairs separated mod 2 below 10^5: ", nmod2, ", largest least prime with a_p(E) != a_p(E') mod 2: ", worstmod2);
\\ ---------------------------------------------------------------- 2. twist pairs (Theorem 2 mechanism)
print("2. twists E_D, D | 2N_E squarefree, D != 1: least p not dividing 2DN_E with chi_D(p) = -1 and a_p(E) != 0");
tw = [[[0,-1,1,-10,-20], "11a1"], [[0,0,1,-1,0], "37a1"], [[1,0,1,-190891,-36002922], "j=-17*373^3/2^17"], [[1,1,1,-8,6], "j=-7*11^3"]];
{
for (t = 1, #tw,
  my(E = ellinit(tw[t][1]), N = ellglobalred(E)[1], R = radi(2*N), mx = 0, Ds = []);
  fordiv (2*N, d, if (!issquarefree(d), next);
    foreach ([d, -d], D, if (D == 1, next);
      my(p0 = 0);
      forprime (p = 3, 10^5, if ((2*D*N) % p == 0, next);
        if (kronecker(D, p) == -1 && ellap(E, p) != 0, p0 = p; break));
      if (p0 > mx, mx = p0)));
  printf("   %-18s N = %-6d rad(2N) = %-5d max over D of least p = %-4d  C_E bound 2(11.5 log rad(2N)+68)^16 = 10^%.2f\n",
    tw[t][2], N, R, mx, log(2*(11.5*log(R) + 68)^16)/LOG10));
}
\\ ---------------------------------------------------------------- 3. non-surjective primes <= 37 lie below C_E
print("3. known non-surjective primes (rational isogenies, Zywina 1508.07660 examples) against the C_E bound");
{
foreach ([[[0,-1,1,-10,-20], 5], [[1,0,1,-190891,-36002922], 17], [[1,1,1,-8,6], 37]], c,
  my(E = ellinit(c[1]), N = ellglobalred(E)[1], M = ellisomat(E)[2], degs = Set(concat(Vec(M))), R = radi(2*N));
  printf("   N = %-6d isogeny degrees %s; non-surjective at %d; bound 2(11.5 log %d + 68)^16 = 10^%.2f > %d\n",
    N, degs, c[2], R, log(2*(11.5*log(R) + 68)^16)/LOG10, c[2]));
}
\\ ---------------------------------------------------------------- 4. Hensel-Serre bound on 2-division fields
print("4. Hensel-Serre: log Delta_L <= (n-1) log rad(d_L) + n log n, rad(d_L) | 2N, for L = Q(E[2], E'[2])");
okHS = 1; nHS = 0; maxratio = 0;
{
for (i = 1, min(#clsrep, 40), for (j = i + 1, min(#clsrep, 40),
  my(E = clsrep[i][2], F = clsrep[j][2], N = 2*clsrep[i][1]*clsrep[j][1], f, g, h, K, d, n);
  f = elldivpol(E, 2); g = elldivpol(F, 2);
  h = polredbest(nfsplitting(f*g/content(f*g)));   \\ splitting field of both 2-division polynomials
  if (poldegree(h) < 2, next);
  K = nfinit(h); d = abs(K.disc); n = poldegree(h);
  nHS++;
  if (N % radi(d) != 0, okHS = 0; print("   rad(d) does not divide 2N: ", [clsrep[i][1], clsrep[j][1]]));
  my(rhs = (n - 1)*log(radi(d)) + n*log(n));
  if (log(d) > rhs, okHS = 0; print("   Hensel-Serre FAILS ", [clsrep[i][1], clsrep[j][1], n]));
  if (log(d)/rhs > maxratio, maxratio = log(d)/rhs)));
}
print("   fields checked: ", nHS, "; all OK: ", okHS, "; max log Delta / bound = ", precision(maxratio, 6));
quit;
