\\ Lemma generator check (Banks 2610.31402): tau in O_k with F = k(sqrt tau) and
\\ |iota(tau)| <= C_k d_F^{1/n_k} for all embeddings; ratio R = max|iota(tau)| / d_F^{1/n_k} must stay bounded.
default(parisize, 64000000);
setrand(3);
test(pol, ntr) =
{
  my(bnf = bnfinit(pol, 1), nk = poldegree(pol), dk = abs(bnf.disc), Rmax = 0., Rmin = 1e100, cnt = 0, big = 0);
  for (t = 1, ntr,
    my(tau0 = vector(nk, i, random(2001) - 1000)~);
    tau0 = nfalgtobasis(bnf, tau0);
    if (tau0 == 0 || #nfroots(bnf, x^2 - nfbasistoalg(bnf, tau0)) > 0, next);
    my(fa = idealfactor(bnf, tau0), a = 1, b = 1);
    for (i = 1, #fa~,
      a = idealmul(bnf, a, idealpow(bnf, fa[i,1], fa[i,2] % 2));
      b = idealmul(bnf, b, idealpow(bnf, fa[i,1], fa[i,2] \ 2)));
    my(c = idealred(bnf, b));                                   \\ small ideal in the class of b
    my(pr = bnfisprincipal(bnf, idealdiv(bnf, b, c)));          \\ (gamma) = b / c
    if (pr[1] != 0, error("not principal"));
    my(g = pr[2]);
    my(tau1 = nfeltdiv(bnf, tau0, nfeltmul(bnf, g, g)));
    \\ unit balancing: minimize max |iota(tau1 u^2)| over u in the unit group (search small exponents)
    my(fu = bnf.fu, best = tau1, bestv = vecmax(abs(nfeltembed(bnf, tau1))));
    if (#fu > 0,
      my(u = nfalgtobasis(bnf, fu[1]));
      for (m = -60, 60,
        my(cand = nfeltmul(bnf, tau1, nfeltpow(bnf, u, 2*m)), v = vecmax(abs(nfeltembed(bnf, cand))));
        if (v < bestv, bestv = v; best = cand)));
    my(dF = dk^2 * idealnorm(bnf, rnfdisc(bnf, x^2 - nfbasistoalg(bnf, best))[1]));
    my(R = bestv / dF^(1/nk));
    cnt++; Rmax = max(Rmax, R); Rmin = min(Rmin, R);
    if (dF > 10^8, big++));
  print("k = Q[y]/(", pol, "): d_k=", dk, " h=", bnf.no, " fields=", cnt, " (d_F>1e8: ", big, ")  max R=", Rmax, "  min R=", Rmin);
}
test(y^2 - 10, 300);
test(y^2 + 5, 300);
test(y^3 - y^2 - 3*y + 1, 150);
quit;

