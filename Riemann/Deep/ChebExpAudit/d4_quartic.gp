\\ chebaudit: D4 quartic fields L, |d_L| <= 10^4: K = L (least degree-one prime) and K = k (quadratic subfield; split, inert),
\\ witnesses p !| poldisc(f) via factormod patterns, against bounds A (0.484, exponent 32) and B (0.62, exponent 1/delta_2).
default(parisize, "64M");
d2 = 1 - (777 - sqrt(1081))/852;
nlin(f, p) = my(F = factormod(f, p)); sum(i = 1, #F~, poldegree(F[i,1]) == 1);
nquad(f, p) = my(F = factormod(f, p)); sum(i = 1, #F~, poldegree(F[i,1]) == 2);
V = nflist("D4", [1, 10^4]); mxA = 0; mxB = 0; bad = 0; ncase = 0;
{
  for (i = 1, #V, my(f = V[i], DL = abs(nfdisc(f)), D = poldisc(f), S = nfsubfields(f, 2), w = [0, 0, 0], dk);
    if (#S != 1, error("D4 subfields"));
    dk = nfdisc(S[1][1]);
    forprime (p = 2, 10^6, if (D % p,
      if (!w[1] && nlin(f, p) > 0, w[1] = p);
      if (kronecker(dk, p) == 1, if (!w[2] && nlin(f, p) > 0, w[2] = p); if (!w[3] && nquad(f, p) > 0, w[3] = p));
      if (vecmin(w) > 0, break)));
    for (j = 1, 3, ncase++; mxA = max(mxA, exp(log(w[j])/32)/log(3*DL)); mxB = max(mxB, exp(d2*log(w[j]))/log(3*DL));
      if (w[j] > (0.484*log(3*DL))^32 || w[j] > (0.62*log(3*DL))^(1/d2), bad++)));
  printf("D4 quartic fields |d| <= 10^4: %d; (L, K, C) cases: %d; violations of A(0.484) or B(0.62): %d; max N^(1/32)/log 3D = %.4f; max N^(d2)/log 3D = %.4f\n", #V, ncase, bad, mxA, mxB);
}
