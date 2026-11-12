\\ cmls.gp (cmaudit, 2026-11-12): per-field tests of CMQuartic Lemma 1 (relative determinant, nontrivial ideals c),
\\ Lemma 2 (O_F-primitive ideals per class), Lemma 3 (exact sharp large sieve constant = top eigenvalue of the dual Gram
\\ over Cl_K^, content weight N(g)^{1+2eps}), plus CM2 Lemma 7 (boxes in ideal lattices of F, skewed boxes).
\\ K = F(sqrt(-alpha)), alpha >> 0. Output: one line per field and test.
\\ run with: gp -q -s 64000000 -D parisizemax=120000000 (any default(parisize*) inside a file discards the later definitions in gp 2.19)

\\ ---------- field construction ----------
mkK(Fpol, alpha) = {
  my(nfF = nfinit(Fpol), R = rnfequation(nfF, x^2 + alpha, 1), z = R[1], a = R[2], k = R[3], pr, P, mp, ay, bet, bnf, nf,
     aut, c = 0, fp, Drel, dF = abs(nfF.disc), dK, D, jcl = []);
  pr = polredbest(z, 1); P = pr[1]; mp = pr[2];
  ay = Mod(subst(lift(a), x, lift(mp)), P);
  bet = Mod(lift(mp), P) - k*ay;
  if (bet^2 + Mod(subst(lift(Mod(alpha, Fpol)), y, lift(ay)), P) != 0, error("beta"));
  bnf = bnfinit(P, 1); nf = bnf.nf;
  aut = nfgaloisconj(nf);
  for (i = 1, #aut, my(s = aut[i]);
    if (s != x && Mod(subst(lift(ay), x, s), P) == ay, c = s));
  if (c == 0, error("no c"));
  if (Mod(subst(lift(bet), x, c), P) != -bet, error("c(beta)"));
  fp = Mod(subst(deriv(Fpol, y), y, lift(ay)), P);
  Drel = idealdiv(nf, nf.diff, lift(fp));
  dK = abs(nf.disc); D = dK / dF^2;
  if (idealnorm(nf, Drel) != D, error("Drel norm"));
  \\ classes of j(Cl_F): extend generators of Cl_F
  my(bnfF = bnfinit(Fpol, 1));
  if (bnfF.no > 1,
    forprime (p = 2, 200, my(dec = idealprimedec(nfF, p));
      for (i = 1, #dec, my(v = bnfisprincipal(bnfF, dec[i], 0));
        if (v != 0*v,
          my(g = Mod(subst(lift(nfbasistoalg(nfF, dec[i][2])), y, lift(ay)), P), J = idealhnf(nf, p, lift(g)));
          jcl = concat(jcl, [bnfisprincipal(bnf, J, 0)]))));
  );
  [bnf, nf, c, ay, bet, Drel, dK, D, dF, nfF, bnfF, jcl]
}

capply(K, e) = { my(nf = K[2], c = K[3], a = nfbasistoalg(nf, e)); nfalgtobasis(nf, Mod(subst(lift(a), x, c), nf.pol)) };

\\ ---------- Lemma 1 with nontrivial ideals ----------
lemma1(K, ntr) = {
  my(bnf = K[1], nf = K[2], c = K[3], Drel = K[6], D = K[8], ps = List(), ok = 0, bad = 0, minr = 10^100, maxid = 0, nprinc = 0);
  forprime (p = 2, 60, my(dec = idealprimedec(nf, p)); for (i = 1, #dec, if (idealnorm(nf, dec[i]) < 4000, listput(ps, dec[i]))));
  for (t = 1, ntr,
    my(Jd = idealhnf(nf, 1), k = 1 + random(3));
    for (j = 1, k, Jd = idealmul(nf, Jd, ps[1 + random(#ps)]));
    my(v = bnfisprincipal(bnf, Jd, 0)); if (v == 0*v, nprinc++);
    my(Jc = nfgaloisapply(nf, c, Jd), T = idealmul(nf, idealmul(nf, Jd, Jc), Drel), NI = idealnorm(nf, Jd));
    for (r = 1, 6,
      my(b1 = Jd * vectorv(4, i, random(21) - 10), b2 = Jd * vectorv(4, i, random(21) - 10));
      my(Dl = nfeltmul(nf, b1, capply(K, b2)) - nfeltmul(nf, capply(K, b1), b2));
      if (Dl == 0*Dl, next);
      my(nd = abs(nfeltnorm(nf, Dl)));
      if (idealhnf(nf, idealadd(nf, T, Dl)) == T, ok++, bad++);
      my(rr = nd / (NI^2 * D)); if (denominator(rr) != 1, bad++); minr = min(minr, rr);
      \\ identity prod_v 2|s b1||s b2||sin| = |N Delta|^{1/2}
      my(e1 = nfeltembed(nf, b1), e2 = nfeltembed(nf, b2), pr = 1.);
      for (v = 1, #e1, pr *= 2 * abs(imag(e1[v] * conj(e2[v]))));
      maxid = max(maxid, abs(pr / sqrt(nd) - 1))));
  [ok, bad, minr, maxid, nprinc]
}

\\ ---------- prime data, DFS over ideals of norm <= Nmax ----------
\\ returns list of [minnorm, localList] per rational prime; local entries [norm, classvec, gnorm, zprim]
primedata(K, Nmax) = {
  my(bnf = K[1], nf = K[2], c = K[3], Drel = K[6], nfF = K[10], out = List(), cyc = bnf.cyc, r = #cyc);
  forprime (p = 2, Nmax,
    my(fm = floor(log(Nmax + 0.5) / log(p)), dec = idealprimedec(nf, p, fm), orb = List(), used = vector(#dec), decF = idealprimedec(nfF, p), ftype);
    if (#dec == 0, next);
    ftype = if(#decF == 2, "s", if(decF[1][3] == 2, "r", "i"));
    for (i = 1, #dec, if (used[i], next);
      my(Pc = nfgaloisapply(nf, c, dec[i]), j = 0, cl = if(r, bnfisprincipal(bnf, dec[i], 0), vectorv(0)), q = idealnorm(nf, dec[i]));
      for (k = 1, #dec, if (k != i && !used[k] && idealval(nf, Pc, dec[k]) > 0, j = k; break));
      used[i] = 1;
      if (j, used[j] = 1; listput(orb, ["s", q, cl, if(r, bnfisprincipal(bnf, dec[j], 0), vectorv(0))]),
        if (idealval(nf, Drel, dec[i]) > 0, listput(orb, ["r", q, cl]), listput(orb, ["i", sqrtint(q), cl]))));
    \\ local factor lists per orbit: [norm, cl, cexp, q]
    my(LL = vector(#orb));
    for (o = 1, #orb, my(t = orb[o][1], q = orb[o][2], cl = orb[o][3], L = List());
      if (t == "s", my(cl2 = orb[o][4]);
          for (a = 0, logint(Nmax, q), for (b = 0, logint(Nmax, q) - a, listput(L, [q^(a+b), a*cl + b*cl2, min(a,b), q]))),
        t == "i", for (a = 0, logint(Nmax, q) \ 2, listput(L, [q^(2*a), a*cl, a, q])),
        for (a = 0, logint(Nmax, q), listput(L, [q^a, a*cl, a\2, q])));
      LL[o] = Vec(L));
    \\ combine orbits at p
    my(comb = List([[1, vectorv(r), 1, vector(#orb)]]));
    for (o = 1, #orb, my(nc = List());
      for (u = 1, #comb, for (w = 1, #LL[o], my(e = LL[o][w], nm = comb[u][1] * e[1]);
        if (nm <= Nmax, my(ce = comb[u][4]); ce[o] = e[3]; listput(nc, [nm, comb[u][2] + e[2], comb[u][3] * e[4]^e[3], ce]))));
      comb = nc);
    my(loc = List());
    for (u = 1, #comb, my(e = comb[u]); if (e[1] == 1, next);
      my(zc);
      if (ftype == "s", zc = if(#orb == 2, min(e[4][1], e[4][2]), 0),
          ftype == "i", zc = e[4][1],
          zc = e[4][1] \ 2);
      listput(loc, [e[1], e[2], e[3], zc == 0]));
    if (#loc, loc = vecsort(Vec(loc), 1); listput(out, [loc[1][1], loc])));
  vecsort(Vec(out), 1)
}

\\ globals for DFS
G_pd = 0; G_N = 0; G_T = 0; G_cyc = 0; G_rad = 0; G_A = 0; G_cnt = 0;
rec(i, nm, cl, g, zp) = {
  my(idx = 1, b = 1);
  for (k = 1, #G_cyc, idx += (cl[k] % G_cyc[k]) * G_rad[k]);
  while (nm > G_T[b], b++);
  G_A[1][idx, b] += 1. / g^2; G_A[2][idx, b] += exp(-1.2 * log(g)); G_A[3][idx, b] += (g == 1); G_A[4][idx, b] += 1; G_A[5][idx, b] += zp;
  G_cnt++;
  for (j = i, #G_pd, my(e = G_pd[j]);
    if (nm * e[1] > G_N, break);
    my(L = e[2]);
    for (w = 1, #L, my(l = L[w]); if (nm * l[1] > G_N, break);
      rec(j + 1, nm * l[1], cl + l[2], g * l[3], zp && l[4])));
}

lemma23(K, fac) = {
  my(bnf = K[1], dK = K[7], D = K[8], dF = K[9], h = bnf.no, cyc = bnf.cyc, sq = sqrt(dK), T, Nmax, rad, A, res = List(), t0 = getabstime());
  T = vector(#fac, k, floor(fac[k] * sq)); Nmax = T[#T];
  rad = vector(#cyc); my(rr = 1); for (k = 1, #cyc, rad[k] = rr; rr *= cyc[k]);
  G_pd = primedata(K, Nmax); my(t1 = getabstime() - t0);
  G_N = Nmax; G_T = T; G_cyc = cyc; G_rad = rad; G_cnt = 0;
  G_A = vector(5, i, matrix(h, #T, a, b, 0.));
  rec(1, 1, vectorv(#cyc), 1, 1);
  \\ cumulative
  for (i = 1, 5, for (b = 2, #T, for (a = 1, h, G_A[i][a, b] += G_A[i][a, b-1])));
  \\ j(Cl_F) cosets
  my(J = K[12], cos = vector(h, a, a));
  if (#J, \\ union-find style: for each class, coset rep = min index over orbit under adding jcl generators
    for (a = 1, h, my(orbit = Set([a]), changed = 1);
      while (changed, changed = 0;
        foreach (orbit, u, for (s = 1, #J, my(vv = vectorv(#cyc), uu = u - 1);
          forstep (k = #cyc, 1, -1, vv[k] = uu \ rad[k]; uu = uu % rad[k]);
          vv += J[s]; my(id2 = 1); for (k = 1, #cyc, id2 += (vv[k] % cyc[k]) * rad[k]);
          if (!setsearch(orbit, id2), orbit = setunion(orbit, [id2]); changed = 1))));
      cos[a] = orbit[1]));
  my(jsize = vecsum(vector(h, a, cos[a] == cos[1])));
  for (b = 1, #T,
    my(N = T[b], mx = vector(5, i, vecmax(G_A[i][, b])), mean5 = vecsum(G_A[1][, b]) / h, Bor, Bnew, mxZ = 0.);
    if (jsize > 1, my(cs = Map()); for (a = 1, h, my(v = 0.); mapisdefined(cs, cos[a], &v); mapput(cs, cos[a], v + G_A[1][a, b]));
      foreach (Mat(cs)[, 2], v, mxZ = max(mxZ, v)), mxZ = mx[1]);
    Bnew = (1 + 4 * Pi * sqrt(N / sqrt(D)))^2;
    listput(res, [N, N / sq, h * mx[1] / (N + sq), h * mx[2] / (N + sq), mx[3], mx[3] / Bnew, h * mx[4] / (N + sq), h * mx[5] / (N + sq), mx[1] / mean5, (h / jsize) * mxZ / (N + sq)]));
  [Vec(res), G_cnt, t1, getabstime() - t0, jsize]
}

\\ ---------- CM2 Lemma 7: boxes in ideal lattices of F (n = 2) ----------
boxtest(Fpol, ntr) = {
  my(nfF = nfinit(Fpol), rts = real(polroots(Fpol)), worst = 0., ps = List(), cnt = 0);
  forprime (p = 2, 40, my(dec = idealprimedec(nfF, p)); for (i = 1, #dec, listput(ps, dec[i])));
  for (t = 1, ntr,
    my(Jd = idealhnf(nfF, 1)); for (j = 1, 1 + random(3), Jd = idealmul(nfF, Jd, ps[1 + random(#ps)]));
    my(NI = idealnorm(nfF, Jd), B = matrix(2, 2), Bi, sk = 10.^(random(13) - 6), A = [0.3, 1, 4, 16][1 + random(4)]);
    \\ embedding of the Z-basis of I: columns = elements
    for (j = 1, 2, my(e = nfbasistoalg(nfF, Jd[, j])); for (v = 1, 2, B[v, j] = subst(lift(e), y, rts[v])));
    Bi = B^(-1);
    my(Z1 = sqrt(A * NI * sk), Z2 = sqrt(A * NI / sk), c1 = (random(10^6) / 10.^6 - 0.5) * 50 * Z1, c2 = (random(10^6) / 10.^6 - 0.5) * 50 * Z2);
    \\ lattice point t = B*k (t0 = 0, shift absorbed in c); |t_v - c_v| <= Z_v
    my(k1lo = -1, k1hi = 1, n = 0);
    \\ bounds on k1 = Bi[1,]*t, t in box
    my(m1 = Bi[1,1]*c1 + Bi[1,2]*c2, w1 = abs(Bi[1,1])*Z1 + abs(Bi[1,2])*Z2);
    for (k1 = ceil(m1 - w1), floor(m1 + w1),
      \\ t_v = B[v,1]k1 + B[v,2]k2 ; need |t_v - c_v| <= Z_v for v = 1, 2: intervals in k2
      my(lo = -10.^30, hi = 10.^30);
      for (v = 1, 2, my(a0 = B[v,1]*k1 - (if(v == 1, c1, c2)), s = B[v,2], Zv = if(v == 1, Z1, Z2), l1 = (-Zv - a0) / s, l2 = (Zv - a0) / s);
        lo = max(lo, min(l1, l2)); hi = min(hi, max(l1, l2)));
      if (hi >= lo, n += max(0, floor(hi) - ceil(lo) + 1)));
    my(bd = 2 * (1 + 8 * Z1 * Z2 / NI));
    worst = max(worst, n / bd); cnt++);
  [cnt, worst]
}

\\ ---------- explicit eigenvalue check for one small field ----------
eigcheck(K, N) = {
  my(bnf = K[1], h = bnf.no, cyc = bnf.cyc, ids, w, cls, G, chars, lam);
  ids = ideallist(bnf.nf, N);
  my(L = List());
  for (m = 1, N, for (i = 1, #ids[m], my(Jd = ids[m][i], cl = bnfisprincipal(bnf, Jd, 0));
    \\ content: largest ideal g of F with gO_K | I : norm of I + Ic ... gO_K = I + Ic? (not in general) -> compute via valuations
    my(fa = idealfactor(bnf.nf, Jd), gn = 1);
    \\ group primes by c-orbits
    my(done = vector(#fa~));
    for (a = 1, #fa~, if (done[a], next); done[a] = 1;
      my(P = fa[a, 1], ea = fa[a, 2], Pc = nfgaloisapply(bnf.nf, K[3], P), eb = 0, isc = 0);
      if (idealval(bnf.nf, Pc, P) > 0, \\ c-fixed
        if (idealval(bnf.nf, K[6], P) > 0, gn *= idealnorm(bnf.nf, P)^(ea \ 2), gn *= sqrtint(idealnorm(bnf.nf, P))^ea),
        for (b = 1, #fa~, if (!done[b] && idealval(bnf.nf, Pc, fa[b, 1]) > 0, done[b] = 1; eb = fa[b, 2]));
        gn *= idealnorm(bnf.nf, P)^min(ea, eb)));
    listput(L, [cl, 1. / gn^2])));
  chars = List(); forvec (e = vector(#cyc, k, [0, cyc[k] - 1]), listput(chars, e));
  my(nc = #chars, V = matrix(nc, #L));
  for (s = 1, nc, for (t = 1, #L, my(ph = 0.); for (k = 1, #cyc, ph += chars[s][k] * L[t][1][k] / cyc[k]);
    V[s, t] = exp(2 * Pi * I * ph) * sqrt(L[t][2])));
  G = V * conj(V~);
  my(RG = matconcat([real(G), -imag(G); imag(G), real(G)]));
  lam = vecmax(qfjacobi(RG)[1]);  \\ real symmetric 2h x 2h form of the Hermitian G (mateigen hangs on repeated eigenvalues)
  my(byc = Map()); for (t = 1, #L, my(v = 0., key = L[t][1]); mapisdefined(byc, key, &v); mapput(byc, key, v + L[t][2]));
  my(mxc = 0.); foreach (Mat(byc)[, 2], v, mxc = max(mxc, v));
  [h, #L, lam, h * mxc, lam / (h * mxc) - 1]
}
