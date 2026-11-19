\\ chebaudit: S3 closures L of all S3 cubic fields F with |d_F| <= 3000, K = Q, classes {1}, transpositions, 3-cycles.
\\ Delta_L = d_F^2 |d_k| (k = Q(sqrt d_F)); witnesses p !| poldisc(f) by the factorisation pattern of f mod p.
default(parisize, "64M");
d2 = 1 - (777 - sqrt(1081))/852;
V = nflist("S3", [1, 3000]); mxA = 0; mxB = 0; bad = 0; ncl = 0;
{
  for (i = 1, #V, my(f = V[i], dF = nfdisc(f), dk = coredisc(dF), DL = dF^2*abs(dk), D = poldisc(f), w = [0, 0, 0], r);
    forprime (p = 2, 10^6, if (D % p,
      r = #polrootsmod(f, p); my(j = if (r == 3, 1, r == 1, 2, 3));
      if (!w[j], w[j] = p); if (vecmin(w) > 0, break)));
    for (j = 1, 3, ncl++; my(a = exp(log(w[j])/32)/log(3*DL), b = exp(d2*log(w[j]))/log(3*DL));
      mxA = max(mxA, a); mxB = max(mxB, b);
      if (w[j] > (0.484*log(3*DL))^32 || w[j] > (0.62*log(3*DL))^(1/d2), bad++)));
  printf("S3 cubic fields |d| <= 3000: %d; classes: %d; violations of A(0.484) or B(0.62): %d; max N^(1/32)/log 3D_L = %.4f; max N^(d2)/log 3D_L = %.4f\n", #V, ncl, bad, mxA, mxB);
}
