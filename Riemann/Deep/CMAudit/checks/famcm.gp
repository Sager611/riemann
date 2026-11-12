\\ famcm.gp (cmaudit): CM3 Theorem 1 family sieve test over F = Q(sqrt5) (h_F = 1, S_0 empty).
\\ Family: all CM K/F with X/2 < N(d_{K/F}) <= X (one K per squarefree ideal a = (alpha), alpha >> 0).
\\ b on admissible n (squarefree, all prime factors of norm >= 19), detector weights as fampos.gp (sigma, V = sqrt(X)/2, Y = X);
\\ LHS = sum_K h_K sum_C |S_C|^2 = sum_K sum_{chi in Cl_K^} |sum_n b(n) lambda_chi(n)|^2.
sstep(t) = if(t <= 0, 0., if(t >= 1, 1., my(a = exp(-1/t), c = exp(-1/(1-t))); a/(a+c)));
phi2(v) = if(v <= 0.25 || v >= 4, 0., my(l = log(v)/log(2)); if(l < -1, sstep(l+2), if(l > 1, sstep(2-l), 1.)));
famrun(X, Nl, sig, OUT) = {
  my(bnfF = bnfinit(y^2 - y - 1, 1), nfF = bnfF.nf, eps = nfbasistoalg(nfF, bnfF.fu[1]), rts = real(polroots(nfF.pol)), Nm = 4 * vecmax(Nl), V = max(2, floor(sqrt(X)/2)), Y = X, l2 = 4.^(1 - sig), t0 = getabstime());
  \\ admissible F-primes with norm <= Nm, and admissible n
  my(il = ideallist(nfF, Nm), Pl = List(), pidx = Map(), adm = List());
  forprime (p = 2, Nm, my(dec = idealprimedec(nfF, p)); for (i = 1, #dec, my(q = idealnorm(nfF, dec[i])); if (q >= 19 && q <= Nm && p != 2, listput(Pl, dec[i]); mapput(pidx, [p, i], #Pl))));
  for (m = 1, Nm, for (i = 1, #il[m], my(I = il[m][i], fa = idealfactor(nfF, I), ok = 1, ix = List());
    for (k = 1, #fa~, my(pr = fa[k, 1]); if (fa[k, 2] > 1 || idealnorm(nfF, pr) < 19 || pr.p == 2, ok = 0; break);
      my(dec = idealprimedec(nfF, pr.p), j = 0); for (u = 1, #dec, if (idealval(nfF, idealhnf(nfF, pr), dec[u]) > 0, j = u));
      listput(ix, mapget(pidx, [pr.p, j])));
    if (ok, listput(adm, [m, Vec(ix)]))));
  adm = Vec(adm);
  \\ coefficient matrix B: rows = sequences (killed, pole) per N, columns = admissible n
  my(K2 = 2 * #Nl, B = matrix(K2, #adm), NB = vector(K2));
  for (k = 1, #Nl, my(N = Nl[k]);
    for (i = 1, #adm, my(n = adm[i][1], pn = vector(#adm[i][2], u, idealnorm(nfF, Pl[adm[i][2][u]])), a = 0);
      forsubset (#pn, S, my(d = prod(u = 1, #S, pn[S[u]])); if (d <= V, a += (-1)^#S));
      my(pre = a * n^(-sig), w1 = exp(-n/Y) * phi2(n/N), w2 = l2 * exp(-4.*n/Y) * phi2(4.*n/N));
      B[2*k-1, i] = pre * (w1 - w2); B[2*k, i] = pre * w1));
  for (k = 1, K2, NB[k] = sum(i = 1, #adm, B[k, i]^2));
  \\ family
  my(sq = il, nK = 0, famh = 0, LHS = vector(K2), TRIV = vector(K2), DIAG = vector(K2));
  for (m = 1, X, for (i = 1, #sq[m], my(I = sq[m][i], fa = idealfactor(nfF, I));
    if (#fa~ && vecmax(fa[, 2]) > 1, next);
    my(g = nfbasistoalg(nfF, bnfisprincipal(bnfF, I)[2]));
    \\ make totally positive: try unit multiples
    my(found = 0); foreach ([1, -1, eps, -eps], u, my(h = g * u); if (subst(lift(h), y, rts[1]) > 0 && subst(lift(h), y, rts[2]) > 0, g = h; found = 1; break));
    if (!found, error("sign"));
    my(al = lift(g), rd = rnfdisc(nfF, x^2 + al), Nd = idealnorm(nfF, rd[1]));
    if (Nd <= X/2 || Nd > X, next);
    nK++;
    my(R = rnfequation(nfF, x^2 + al, 1), pr = polredbest(R[1], 1), P = pr[1], ay = Mod(subst(lift(R[2]), x, lift(pr[2])), P), bnf = bnfinit(P, 1), nfK = bnf.nf, h = bnf.no, cyc = bnf.cyc, rad = vector(#cyc), rr = 1);
    for (k = 1, #cyc, rad[k] = rr; rr *= cyc[k]);
    famh += h;
    \\ type and class per admissible F-prime: 0 inert, [1, cl] split (cl of one prime above), [2, cl] ramified
    my(typ = vector(#Pl), done = Map());
    for (t = 1, #Pl, my(pf = Pl[t], p = pf.p, decK);
      if (!mapisdefined(done, p, &decK), decK = idealprimedec(nfK, p); mapput(done, p, decK));
      my(gK = subst(lift(nfbasistoalg(nfF, pf[2])), y, lift(ay)), above = List());
      for (u = 1, #decK, if (p == 0 || nfeltval(nfK, gK, decK[u]) > 0 || (pf[2] == 0), listput(above, decK[u])));
      if (#above == 0, \\ generator p alone (inert in F): all K-primes above p lie above pf
        above = List(decK));
      if (#above == 2, typ[t] = [1, if(#cyc, bnfisprincipal(bnf, above[1], 0), vectorv(0))],
        if (above[1][3] == 2 * pf[3], typ[t] = [2, if(#cyc, bnfisprincipal(bnf, above[1], 0), vectorv(0))], typ[t] = 0)));
    my(S = matrix(K2, max(h, 1)));
    for (i = 1, #adm, my(ix = adm[i][2], spl = List(), base = vectorv(#cyc), ok = 1);
      for (u = 1, #ix, my(tt = typ[ix[u]]); if (tt == 0, ok = 0; break); if (tt[1] == 1, listput(spl, tt[2]), base += tt[2]));
      if (!ok, next);
      my(s = #spl, col = B[, i]);
      for (k = 1, K2, DIAG[k] += h * col[k]^2 * 2^s);
      forvec (e = vector(s, u, [0, 1]), my(v = base); for (u = 1, s, v += if(e[u], spl[u], -spl[u]));
        my(id = 1); for (k = 1, #cyc, id += (v[k] % cyc[k]) * rad[k]);
        for (k = 1, K2, S[k, id] += col[k])));
    for (k = 1, K2, my(row = S[k, ]); LHS[k] += h * sum(c = 1, h, row[c]^2); TRIV[k] += vecsum(row)^2)));
  for (k = 1, #Nl,
    my(a = 2*k-1, b2 = 2*k);
    write(OUT, Strprintf("X=%d fields=%d sum_h=%d N=%d  killed: R1=%.3f RX=%.4f R1diag=%.3f triv=%.2f   pole: R1=%.3f RX=%.4f R1diag=%.3f triv=%.2f  ms=%d",
      X, nK, famh, Nl[k], LHS[a]/(famh*NB[a]), LHS[a]/(X^1.5*NB[a]), LHS[a]/DIAG[a], TRIV[a]/LHS[a], LHS[b2]/(famh*NB[b2]), LHS[b2]/(X^1.5*NB[b2]), LHS[b2]/DIAG[b2], TRIV[b2]/LHS[b2], getabstime()-t0)));
}
