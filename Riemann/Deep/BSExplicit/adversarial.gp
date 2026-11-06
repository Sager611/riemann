\\ bsexpl adversarial checks (PARI/GP 2.19, single thread). Run: perl -e 'alarm 1500; exec @ARGV' gp -q adversarial.gp > adversarial_output.txt
default(nbthreads, 1); default(parisizemax, 800000000); default(realprecision, 38);
c4 = 2*exp(21/8 - (2/log(3))/8*psi(1/2)); cqd = (3 - 2*sqrt(2))/2; cbs = 2*cqd/c4;
NMAX = 20000;
\\ quartic fields, r2 = 2, |d| <= NMAX, containing k = Q(sqrt(dk))
containsk(P, dk) = { my(S = nfsubfields(P, 2)); for (i = 1, #S, if (nfdisc(S[i][1]) == dk, return(1))); 0; }
fams = [["V4"], ["C4"], ["D4"]];
{
for (ik = 1, 2,
  my(dk = [-4, -3][ik], kpol = [x^2 + 1, x^2 + x + 1][ik], Lk = lfuncreate(kpol), allF = List(), mind = 10^9, nreal = 0, worst = 10.);
  for (ig = 1, 3, my(G = fams[ig][1], Ls = nflist(G, [1, NMAX], 2));
    for (j = 1, #Ls, if (containsk(Ls[j], dk), listput(allF, [G, Ls[j]]))));
  print("k = Q(sqrt(", dk, ")): ", #allF, " quartic F/k with d_F <= ", NMAX);
  for (j = 1, #allF,
    my(G = allF[j][1], P = allF[j][2], dF = abs(nfdisc(P)), LF = lfuncreate(P), LX = lfundiv(LF, Lk), vals, sg, zs = List());
    mind = min(mind, dF);
    \\ L(sigma, chi_{F/k}) on [1/2, 1): sign scan (a real zero of odd order changes sign), plus value at 1
    vals = vector(26, i, real(lfun(LX, 1/2 + (i-1)/50)));
    vals = concat(vals, [real(lfun(LX, 0.99)), real(lfun(LX, 0.999)), real(lfun(LX, 1))]);
    for (i = 2, #vals, if (vals[i-1]*vals[i] < 0, listput(zs, i)));
    if (#zs, nreal++; print("  REAL ZERO sign change: ", G, " ", P, " d_F=", dF, " at index ", zs));
    worst = min(worst, vecmin(vals));
    if (j <= 6 || #zs, print("  ", G, " d_F=", dF, " P=", P, "  min_{[1/2,1]} L(sigma,chi) = ", vecmin(vals), "  L(1,chi) = ", vals[#vals],
      "  forbidden zone [1 - c/log d_F, 1) = [", 1 - cqd/log(dF), ", 1)")));
  print("  min d_F = ", mind, "; fields with a sign change on [1/2,1): ", nreal, "; min of L(sigma,chi) over all samples: ", worst));
}
\\ D4: Galois closure M, discriminant bound d_M <= (d_F/d_k)^4, double-zero structure zeta_M = zeta_k L(chi_F)^2 L(chi_F3)
{
my(tests = [[x^2+1, -4], [x^2+x+1, -3]], cnt);
for (ik = 1, 2, my(kpol = tests[ik][1], dk = tests[ik][2], Ls = nflist("D4", [1, NMAX], 2), Lk = lfuncreate(kpol)); cnt = 0;
  for (j = 1, #Ls, if (cnt >= 4, break); my(P = Ls[j]); if (!containsk(P, dk), next); cnt++;
    my(dF = abs(nfdisc(P)), M = polredbest(nfsplitting(P)), dM = abs(nfdisc(M)), Q4 = nfsubfields(M, 4), F3 = 0, LM, LF, L3, lhs, rhs);
    for (i = 1, #Q4, my(R = Q4[i][1]); if (containsk(R, dk) && #nfsubfields(R, 2) == 3, F3 = R));
    LM = lfuncreate(M); LF = lfuncreate(P); L3 = lfuncreate(F3);
    lhs = lfun(LM, 2); rhs = lfun(Lk, 2) * (lfun(LF, 2)/lfun(Lk, 2))^2 * (lfun(L3, 2)/lfun(Lk, 2));
    print("D4 over Q(sqrt(", dk, ")): d_F=", dF, " [M:Q]=", poldegree(M), " d_M=", dM, " (d_F/d_k)^4=", (dF/abs(dk))^4,
      " ok=", dM <= (dF/abs(dk))^4, " log d_M/log d_F=", log(dM)/log(dF),
      "  zeta_M(2)/(zeta_k L(chi_F)^2 L(chi_F3))(2) - 1 = ", lhs/rhs - 1)));
}
\\ Stark Lemma 3 / (9)-(10) sanity on one D4 octic M: sum over computed zeros with multiplicity vs 1/(s-1) + log(d_M)/2
{
my(P = 0, dk = -4, kpol = x^2+1);
my(Ls = nflist("D4", [1, NMAX], 2)); for (j = 1, #Ls, if (containsk(Ls[j], dk), P = Ls[j]; break));
my(M = polredbest(nfsplitting(P)), dM = abs(nfdisc(M)), LM = lfuncreate(M), Z = lfunzeros(LM, 20), s);
for (t = 1, 3, s = [1.05, 1.2, 1.6][t];
  my(S = sum(i = 1, #Z, 2*real(1/(s - (1/2 + I*Z[i])))));   \\ zeros rho, conj(rho) on the line, |gamma| <= 20 (lfunzeros lists gamma >= 0 with multiplicity)
  print("Stark p.140 check, M = closure of ", P, ", d_M=", dM, ", s=", s, ": partial zero sum (|gamma|<=20) = ", S, " < 1/(s-1) + log(d_M)/2 = ", 1/(s-1) + log(dM)/2, "  ", S < 1/(s-1) + log(dM)/2));
}
\\ Stark Lemma 4 / Theorem B sanity: Galois K over k (biquadratic F, octic closures M); bnfinit is GRH-conditional (sanity only)
{
my(cands = List());
for (ik = 1, 2, my(dk = [-4, -3][ik], c = 0, V = nflist("V4", [1, NMAX], 2), D = nflist("D4", [1, NMAX], 2));
  for (j = 1, #V, if (c >= 3, break); if (containsk(V[j], dk), c++; listput(cands, [dk, V[j]])));
  c = 0; for (j = 1, #D, if (c >= 2, break); if (containsk(D[j], dk), c++; listput(cands, [dk, polredbest(nfsplitting(D[j]))]))));
for (j = 1, #cands, my(dk = cands[j][1], K = cands[j][2], b = bnfinit(K, 1), n = poldegree(K), d = abs(b.disc), h = b.no, R = b.reg, w = b.tu[1],
    r2 = b.sign[2], r1 = b.sign[1], kap = 2^r1*(2*Pi)^r2*h*R/(w*sqrt(d)), lowk = 1/(4*c4*log(d)), bs = cbs*sqrt(d)/((2*Pi)^(n/2)*log(d)));
  print("K over Q(sqrt(", dk, ")), [K:Q]=", n, " d_K=", d, " h=", h, " R=", R, " w=", w, ": kappa=", kap, " > 1/(4 c4 log d)=", lowk, " ", kap > lowk,
    ";  hR=", h*R, " > c_BS d^(1/2)/((2pi)^(n/2) log d)=", bs, " ", h*R > bs));
}
quit;
