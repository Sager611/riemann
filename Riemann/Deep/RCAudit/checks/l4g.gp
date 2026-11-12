\\ Lemma 4 audit, any D = f^2 dK: S = sum_{chi ring class char != 1} |L(1/2+iu,chi) M_X(1/2+iu,chi)|^2,
\\ L(s,chi) = L(s,chi*) E_chi(s) (imprimitive at p | f), M_X over ideals prime to f.
chiv(v, cyc, dl) = exp(2 * Pi * I * sum(k = 1, #cyc, v[k] * dl[k] / cyc[k]));
lem4g(dK, f, Xs, us, out) =
{
  my(D = dK * f^2, bnf = bnfinit(x^2 - dK, 1), bnr = bnrinit(bnf, f), cyc = bnr.cyc, Hs, ring = List(), Xm = vecmax(Xs),
     IL, ids = List(), res, k, pr = List(), cache = Map(), T0 = getabstime());
  Hs = [bnrisprincipal(bnr, a, 0) | a <- [1..f], gcd(a, f) == 1];
  forvec(v = vector(#cyc, i, [0, cyc[i] - 1]), if (sum(j = 1, #Hs, abs(chiv(v, cyc, Hs[j]) - 1)) < 1e-9, listput(ring, v)));
  if (#ring != qfbclassno(D), error("ring class count mismatch ", D));
  IL = ideallist(bnf, Xm);
  for (nn = 1, Xm, if (gcd(nn, f) == 1, foreach(IL[nn], J0, my(fa = idealfactor(bnf, J0), mu);
       mu = if (#fa~ == 0, 1, if (vecmax(fa[,2]) > 1, 0, (-1)^(#fa~)));
       if (mu, listput(ids, [nn, mu, bnrisprincipal(bnr, J0, 0)])))));
  foreach(factor(f)[,1], p, foreach(idealprimedec(bnf, p), P, listput(pr, P)));
  res = vector(#Xs * #us);
  for (t = 1, #ring,
    my(v = ring[t], f0, bnr0, chi0, L, key);
    if (v == 0 * v, next);
    f0 = bnrconductor(bnr, v); key = Str(f0);
    if (mapisdefined(cache, key), bnr0 = mapget(cache, key), bnr0 = bnrinit(bnf, f0); mapput(cache, key, bnr0));
    chi0 = bnrmap(bnrmap(bnr, bnr0), v);
    L = lfuncreate([bnr0, chi0]);
    for (iu = 1, #us,
      my(s = 1/2 + I * us[iu], Lv = lfun(L, s), E = 1.);
      for (j = 1, #pr, if (idealval(bnf, f0, pr[j]) == 0,
          E *= 1 - chiv(chi0, bnr0.cyc, bnrisprincipal(bnr0, pr[j], 0)) * idealnorm(bnf, pr[j])^(-s)));
      for (ix = 1, #Xs,
        my(M = 0.);
        for (j = 1, #ids, if (ids[j][1] <= Xs[ix], M += ids[j][2] * chiv(v, cyc, ids[j][3]) * ids[j][1]^(-s)));
        k = (iu - 1) * #Xs + ix; res[k] += abs(Lv * E * M)^2)));
  write(out, [D, dK, f, #ring, Xs, us, res, (getabstime() - T0) / 1000.]);
}
