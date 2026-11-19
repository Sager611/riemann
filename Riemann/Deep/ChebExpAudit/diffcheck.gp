\\ chebaudit: for the rows where ChebExp's least prime (PARI idealprimedec method, as direct.gp) is below the own
\\ factormod witness, check that it divides poldisc(f) (so Dedekind-Kummer could not see it)
default(parisize, "64M");
nlin(f, p) = my(F = factormod(f, p)); sum(i = 1, #F~, poldegree(F[i,1]) == 1);
nquad(f, p) = my(F = factormod(f, p)); sum(i = 1, #F~, poldegree(F[i,1]) == 2);
lp1(nf) = forprime(p = 2, 10^6, my(dec = idealprimedec(nf, p)); for (i = 1, #dec, if (dec[i].f == 1, return(p))));
w1(f) = my(D = poldisc(f)); forprime(p = 2, 10^6, if (D % p && nlin(f, p) > 0, return(p)));
pariq(k, rnf, p) =
{
  my(dec = idealprimedec(k, p), sp = 0, in = 0, DD);
  for (j = 1, #dec, if (dec[j].f == 1, DD = rnfidealprimedec(rnf, dec[j]);
    if (#DD == 2, sp = 1); if (#DD == 1 && DD[1].e == 1, in = 1)));
  [sp, in];
}
cntL = 0; okL = 0; cntQ = 0; okQ = 0;
dorow(G, f) =
{
  my(nf = nfinit(f), D = poldisc(f), a = lp1(nf), b = w1(f));
  if (a < b, cntL++; if (D % a == 0, okL++, print("K=L ", G, " ", f, " ", a, " ", b)));
  if (poldegree(f) != 4, return);
  my(S = nfsubfields(nf, 2));
  for (s = 1, #S,
    my(k = nfinit(subst(S[s][1], x, y)), P = nffactor(k, f)[1, 1], rnf = rnfinit(k, P), dk = k.disc, ps = 0, pi = 0, qs = 0, qi = 0, t);
    forprime (p = 2, 10^5,
      t = pariq(k, rnf, p);
      if (t[1] && !ps, ps = p); if (t[2] && !pi, pi = p);
      if (D % p && kronecker(dk, p) == 1, if (!qs && nlin(f, p) > 0, qs = p); if (!qi && nquad(f, p) > 0, qi = p));
      if (ps && pi && qs && qi, break));
    if (ps < qs, cntQ++; if (D % ps == 0, okQ++, print("split ", G, " ", f, " ", dk, " ", ps, " ", qs)));
    if (pi < qi, cntQ++; if (D % pi == 0, okQ++, print("inert ", G, " ", f, " ", dk, " ", pi, " ", qi))));
}
foreach (["C3", "S3", "C4", "V4", "D4", "S4"], G, V = nflist(G, [1, 257]); for (i = 1, #V, dorow(G, V[i])));
printf("K=L rows with PARI least < own witness: %d, of which the PARI prime divides poldisc(f): %d\n", cntL, okL);
printf("K=quad rows with PARI least < own witness: %d, of which the PARI prime divides poldisc(f): %d\n", cntQ, okQ);
