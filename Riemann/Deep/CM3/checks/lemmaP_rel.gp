\\ lemmaP_rel.gp -- numerical check of the relative Lemma P (Riemann/Deep/CM3/NOTE.md, section 2)
\\ for CM fields K = F(sqrt(-beta)) over real quadratic F, orders O_f = O_F + f O_K (no coprimality to f).
\\ For each (K, f):  Q'_{K,f}(b) := sum_{n,n'} b(n) conj(b(n')) sum_{c | (n,n')} eta_f(c) sum_{u in R(m)} #{alpha in c r O_f : N alpha = m u},
\\ m = nu_n nu_n' (nu_n generates n r, times a random unit), counted by direct enumeration of alpha (no class groups), against
\\ sum over squarefree e | f of omega(K, f/e) * sum_{C in G_{K,f/e}} |S_C(b_e)|^2  (ring class groups from bnrinit / image of (O_F/f')^*).
\\ Also checked: pairs with different f-parts give 0 (step a); per pair, the count in c r O_f equals the count in e c r O_{f/e} (step b);
\\ |G_{K,f}| against h_K |(O_K/f)^*| / (|(O_F/f)^*| [O_K^*:O_f^*]).
default(parisizemax, 256*10^6);
alarm(2400);
setrand(20261021);

\\ ---------- base field: real quadratic, zk = [1, y], y^2 = QS*y + QR ----------
qmul(a, b) = [a[1]*b[1] + QR*a[2]*b[2], a[1]*b[2] + a[2]*b[1] + QS*a[2]*b[2]];
qnorm(a) = a[1]^2 + QS*a[1]*a[2] - QR*a[2]^2;
qtr(a) = 2*a[1] + QS*a[2];
qemb(a, v) = a[1] + a[2]*RHO[v];
qpm(a) = Mod(a[1] + a[2]*y, FPOL);
qconj(a) = [a[1] + QS*a[2], -a[2]];
istp(a) = qnorm(a) > 0 && qtr(a) > 0;
col2(v) = if (type(v) == "t_INT", [v, 0], [v[1], v[2]]);

setupF(polF) = {
  BNF = bnfinit(polF, 1); FPOL = polF;
  if (BNF.zk != [1, y], error("zk"));
  my(c = Vec(polF)); QS = -c[2]; QR = -c[3];
  my(rt = real(polroots(polF))); RHO = [vecmax(rt), vecmin(rt)];
  EPS = col2(nfalgtobasis(BNF, BNF.fu[1]));
  EPSI = qconj(EPS) * qnorm(EPS);   \\ inverse (norm +-1)
  UREPS = [[1,0], [-1,0], EPS, -EPS];
};
qunitpow(k) = { my(r = [1,0], e = if (k >= 0, EPS, EPSI)); for (i = 1, abs(k), r = qmul(r, e)); r };

\\ ---------- K ----------
setupK(beta) = {
  BETA = beta; NBETA = qnorm(beta);
  RNF = rnfinit(BNF, x^2 + (beta[1] + beta[2]*y));
  KB = bnfinit(RNF.polabs, 1);
  LATT = Map(); PRK = Map();
};
emb(a) = rnfeltreltoabs(RNF, Mod(a[1] + a[2]*y, FPOL));
emb2(z) = { my(w = nfbasistoalg(BNF, z)); if (type(w) != "t_POLMOD", w = Mod(w, FPOL)); rnfeltreltoabs(RNF, w) };
embid(id) = {
  my(H = idealhnf(BNF, id));
  if (H == matid(2), return(idealhnf(KB, 1)));
  my(t = idealtwoelt(BNF, H));
  idealadd(KB, idealhnf(KB, emb2(t[1])), idealhnf(KB, emb2(t[2])))
};
\\ inverse of a Z-basis matrix of the lattice l + l f O_K  (= l O_f), in KB coordinates
lattinv(l, f) = {
  my(key = [idealhnf(BNF, l), idealhnf(BNF, f)], M);
  if (mapisdefined(LATT, key, &M), return(M));
  my(H = key[1], cols = List(), J = embid(idealmul(BNF, l, f)));
  for (j = 1, 2, listput(cols, nfalgtobasis(KB, emb2(H[,j]))));
  for (j = 1, 4, listput(cols, J[,j]));
  M = mathnf(Mat(Vec(cols)));
  if (matsize(M) != [4,4], error("lattice rank"));
  M = M^(-1); mapput(LATT, key, M); M
};
inlatt(Mi, v) = denominator(Mi * v) == 1;
\\ primes of K above a prime pr of F: [eta, list of primes P with N_{K/F}P = pr]
primesabove(pr) = {
  my(key = idealhnf(BNF, pr), r);
  if (mapisdefined(PRK, key, &r), return(r));
  my(fa = idealfactor(KB, embid(pr)));
  if (#fa~ == 2, r = [1, [fa[1,1], fa[2,1]]], if (fa[1,2] == 2, r = [0, [fa[1,1]]], r = [-1, []]));
  mapput(PRK, key, r); r
};
\\ all alpha in K with N_{K/F} alpha = mu (mu >> 0) and Tr alpha in l; KB coordinate vectors
enumalpha(mu, l) = {
  my(H = idealhnf(BNF, l), c1 = [H[1,1], H[2,1]], c2 = [H[1,2], H[2,2]]);
  my(B1 = 2*sqrt(qemb(mu,1))*(1+1e-12) + 1e-9, B2 = 2*sqrt(qemb(mu,2))*(1+1e-12) + 1e-9);
  my(E = [qemb(c1,1), qemb(c2,1); qemb(c1,2), qemb(c2,2)], Ei = E^(-1));
  my(umax = floor(abs(Ei[1,1])*B1 + abs(Ei[1,2])*B2) + 1, res = List(), fm = 4*mu);
  for (u = -umax, umax,
    my(l1 = (-B1 - E[1,1]*u)/E[1,2], h1 = (B1 - E[1,1]*u)/E[1,2], l2 = (-B2 - E[2,1]*u)/E[2,2], h2 = (B2 - E[2,1]*u)/E[2,2]);
    if (l1 > h1, [l1, h1] = [h1, l1]); if (l2 > h2, [l2, h2] = [h2, l2]);
    for (v = ceil(max(l1, l2)), floor(min(h1, h2)),
      my(t = u*c1 + v*c2, d0 = fm - qmul(t, t));
      if (d0 == [0,0], listput(res, nfalgtobasis(KB, emb(t/2))),
        if (istp(d0) && issquare(qnorm(d0)*NBETA),
          my(rr = nfroots(BNF, x^2 - qpm(d0)/qpm(BETA)));
          if (#rr,
            my(cc = rr[1], tt = qpm(t));
            listput(res, nfalgtobasis(KB, rnfeltreltoabs(RNF, Mod(tt/2 + cc/2*x, RNF.pol))));
            listput(res, nfalgtobasis(KB, rnfeltreltoabs(RNF, Mod(tt/2 - cc/2*x, RNF.pol)))))))));
  Vec(res)
};
\\ omega(K, f) = |mu(O_f)| * q(O_f)
omegaof(f) = {
  my(Mi = lattinv(1, f), tu = KB.tu, z = tu[2], nm = 0, q = 1);
  for (k = 0, tu[1]-1, if (inlatt(Mi, nfalgtobasis(KB, z^k)), nm++));
  my(fu = KB.fu[1], Nfu = rnfeltnorm(RNF, rnfeltabstorel(RNF, fu)));
  if (#nfroots(BNF, x^2 - Nfu) == 0,   \\ N(fu) not in U^2: Hasse index 2 for O_K
    for (k = 0, tu[1]-1, if (inlatt(Mi, nfalgtobasis(KB, z^k*fu)), q = 2; break)));
  [nm, q, nm*q]
};
\\ [O_K^* : O_f^*] = [mu_K : mu(O_f)] * kmin
unitindex(f) = {
  my(Mi = lattinv(1, f), tu = KB.tu, z = tu[2], fu = KB.fu[1], nm = 0);
  for (k = 0, tu[1]-1, if (inlatt(Mi, nfalgtobasis(KB, z^k)), nm++));
  for (kk = 1, 400, for (k = 0, tu[1]-1, if (inlatt(Mi, nfalgtobasis(KB, z^k*fu^kk)), return((tu[1]/nm)*kk))));
  error("unit index > 400")
};
hnfreduce(v, L) = { my(w = v); forstep (i = #w, 1, -1, w -= floor(w[i]/L[i,i]) * L[,i]); w };

\\ ---------- one test configuration ----------
\\ plist: primes of F for the support; Aclass: class vector; rr: reference ideal r (in A^{-1}); fl: list of conductors (prime to rr)
NTEST = 0; NBAD = 0; NBLK = 0; NBLKNZ = 0; NBLKNZE = 0; NOFF = 0; NSTEPB = 0; NGOK = 0; NG = 0;
runconfig(name, plist, Aclass, rr, Nmax, fl, betas) = {
  my(k = #plist, supp = List());
  \\ squarefree products of the primes in plist, in class A, norm <= Nmax
  forvec (sel = vector(k, i, [0,1]),
    my(id = idealhnf(BNF, 1), nm = 1, mask = 0);
    for (i = 1, k, if (sel[i], id = idealmul(BNF, id, plist[i]); nm *= idealnorm(BNF, plist[i]); mask += 2^(i-1)));
    if (nm > Nmax, next);
    if (Vec(bnfisprincipal(BNF, id, 0)) != Vec(Aclass), next);
    my(g = bnfisprincipal(BNF, idealmul(BNF, id, rr), 1));
    if (g[1] != 0*g[1], error("not principal"));
    my(nu = qmul(col2(g[2]), qmul(qunitpow(random(5)-2), [(-1)^random(2), 0])));
    my(bval = (random(7)-3) + I*(random(7)-3));
    listput(supp, [id, mask, nm, nu, bval]));
  supp = Vec(supp); my(S = #supp);
  printf("%s: |supp| = %d, Nmax = %d\n", name, S, Nmax);
  for (ib = 1, #betas,
    setupK(betas[ib]);
    my(etas = vector(k, i, primesabove(plist[i])[1]));
    \\ enumerate the alpha lists once per pair and unit
    my(alist = Map());
    for (i = 1, S, for (j = i, S,
      my(m = qmul(supp[i][4], supp[j][4]), L = List());
      for (iu = 1, 4, my(mu = qmul(m, UREPS[iu])); if (istp(mu), listput(L, enumalpha(mu, rr))));
      mapput(alist, [i,j], Vec(L))));
    for (jf = 1, #fl,
      my(f = fl[jf], fmask = 0, fa = idealfactor(BNF, f));
      for (i = 1, k, if (idealval(BNF, f, plist[i]) > 0, fmask += 2^(i-1)));
      my(lhs = Map(), offbad = 0, stepbbad = 0);
      for (i = 1, S, for (j = i, S,
        my(gm = bitand(supp[i][2], supp[j][2]), ei = bitand(supp[i][2], fmask), ej = bitand(supp[j][2], fmask));
        my(w = if (i == j, norm(supp[i][5]), 2*real(supp[i][5]*conj(supp[j][5]))));
        my(tot = 0, AL = mapget(alist, [i,j]));
        \\ e-ideal and f' for step (b)
        my(eid = idealhnf(BNF, 1)); for (t = 1, k, if (bitand(ei, 2^(t-1)), eid = idealmul(BNF, eid, plist[t])));
        my(fp = idealdiv(BNF, f, eid));
        forsubset (k, cs,   \\ c = product over subset cs of plist, must divide gm
          my(cm = 0); for (t = 1, #cs, cm += 2^(cs[t]-1));
          if (bitand(cm, gm) != cm, next);
          my(eta = 1, cid = idealhnf(BNF, 1));
          for (t = 1, #cs, my(pi = cs[t]); eta *= if (bitand(fmask, 2^(pi-1)), 0, etas[pi]); cid = idealmul(BNF, cid, plist[pi]));
          if (eta == 0, next);
          my(Mi = lattinv(idealmul(BNF, cid, rr), f), cnt = 0);
          for (iu = 1, #AL, for (a = 1, #AL[iu], if (inlatt(Mi, AL[iu][a]), cnt++)));
          if (ei == ej,
            my(Mi2 = lattinv(idealmul(BNF, idealmul(BNF, eid, cid), rr), fp), cnt2 = 0);
            for (iu = 1, #AL, for (a = 1, #AL[iu], if (inlatt(Mi2, AL[iu][a]), cnt2++)));
            NSTEPB++; if (cnt2 != cnt, stepbbad++));
          tot += eta*cnt);
        if (ei != ej, NOFF++; if (tot != 0, offbad++), 
          my(old = 0); mapisdefined(lhs, ei, &old); mapput(lhs, ei, old + w*tot))));
      \\ right side, block by block
      my(keys = if (#lhs, Vec(Mat(lhs)[,1]), []), ok = (offbad == 0 && stepbbad == 0), totl = 0, totr = 0, desc = "");
      for (ik = 1, #keys,
        my(e = keys[ik], eid = idealhnf(BNF, 1));
        for (t = 1, k, if (bitand(e, 2^(t-1)), eid = idealmul(BNF, eid, plist[t])));
        my(fp = idealdiv(BNF, f, eid), md = embid(fp), bnr = bnrinit(KB, md));
        my(gens = idealstar(BNF, fp, 2).gen, Hs = List());
        for (t = 1, #gens, listput(Hs, bnrisprincipal(bnr, idealhnf(KB, emb2(gens[t])), 0)));
        my(r = #bnr.cyc, Lm);
        if (r == 0, Lm = matrix(0,0), Lm = if (#Hs, mathnf(concat(matdiagonal(bnr.cyc), Mat(Vec(Hs)))), mathnf(matdiagonal(bnr.cyc))));
        my(Gsz = if (r == 0, 1, prod(t = 1, r, Lm[t,t])));
        \\ sanity: class number formula for the ring class group
        my(ui = unitindex(fp), phK = idealstar(KB, md).no, phF = idealstar(BNF, fp).no);
        NG++; if (Gsz * phF * ui == KB.no * phK, NGOK++);
        my(Sc = Map());
        for (i = 1, S, if (bitand(supp[i][2], fmask) != e, next);
          my(opts = List([idealhnf(KB, 1)]), dead = 0);
          for (t = 1, k, if (bitand(supp[i][2], 2^(t-1)) && !bitand(e, 2^(t-1)),
            my(pa = primesabove(plist[t])[2]); if (#pa == 0, dead = 1; break);
            my(no = List()); for (s = 1, #opts, for (q = 1, #pa, listput(no, idealmul(KB, opts[s], pa[q])))); opts = no));
          if (dead, next);
          for (s = 1, #opts,
            my(lab = if (r == 0, [], hnfreduce(bnrisprincipal(bnr, opts[s], 0), Lm)), old = 0);
            mapisdefined(Sc, lab, &old); mapput(Sc, lab, old + supp[i][5])));
        my(ssum = 0, vals = if (#Sc, Mat(Sc)[,2], []));
        for (t = 1, #vals, ssum += norm(vals[t]));
        my(om = omegaof(fp), rhs = om[3]*ssum, lv = mapget(lhs, e));
        NBLK++; if (lv != 0, NBLKNZ++; if (e != 0, NBLKNZE++));
        if (lv != rhs, ok = 0);
        totl += lv; totr += rhs;
        desc = concat(desc, Str(" e", e, ":", lv, "=", rhs, "(w", om[1], ",q", om[2], ",|G|", Gsz, ")")));
      NTEST++; if (!ok, NBAD++);
      printf("  beta=%s h_K=%d f=%s N=%d: %s total %s=%s off-bad=%d stepb-bad=%d%s\n",
        BETA, KB.no, idealfactor(BNF, f)[,2]~, idealnorm(BNF, f), if (ok, "OK", "FAIL"), totl, totr, offbad, stepbbad, desc)));
};

\\ ================= F = Q(sqrt 5) =================
{
setupF(y^2 - y - 1);
P5 = idealprimedec(BNF, 5)[1]; P3 = idealprimedec(BNF, 3)[1];
P11 = idealprimedec(BNF, 11); P19 = idealprimedec(BNF, 19); P29 = idealprimedec(BNF, 29)[1]; P31 = idealprimedec(BNF, 31)[1];
P2 = idealprimedec(BNF, 2)[1];
plist = [P5, P3, P11[1], P11[2], P19[1], P19[2], P29, P31];
fl = [P11[1], idealpow(BNF, P11[1], 2), idealmul(BNF, P11[1], P19[1]), P3, idealpow(BNF, P3, 2),
      idealmul(BNF, P11[1], P2), P5, idealmul(BNF, P11[1], P11[2]), idealmul(BNF, idealpow(BNF, P11[1], 2), P5)];
betas = [[1,0], [3,0], [2,1], [7,0], [5,2], [11,3]];
runconfig("Q(sqrt5)", plist, [], idealhnf(BNF, 1), 200, fl, betas);
}


\\ ================= F = Q(sqrt 3): U+ != U^2, |R(m)| = 2, Hasse index 2 for Q(zeta_12) =================
{
setupF(y^2 - 3);
Q3 = idealprimedec(BNF, 3)[1]; Q11 = idealprimedec(BNF, 11); Q13 = idealprimedec(BNF, 13); Q23 = idealprimedec(BNF, 23)[1];
Q5 = idealprimedec(BNF, 5)[1]; Q2 = idealprimedec(BNF, 2)[1];
plist = [Q3, Q11[1], Q11[2], Q13[1], Q13[2], Q23, Q5];
fl = [Q11[1], idealpow(BNF, Q11[1], 2), idealmul(BNF, Q11[1], Q13[1]), Q3, idealmul(BNF, Q2, Q11[1]), Q5, idealmul(BNF, Q11[1], Q11[2])];
betas = [[1,0], [2,0], [2,1], [5,0], [7,2]];
runconfig("Q(sqrt3)", plist, [], idealhnf(BNF, 1), 160, fl, betas);
}


\\ ================= F = Q(sqrt 10): h_F = 2 (non-principal e, reference ideal r) =================
{
setupF(y^2 - 10);
R3 = idealprimedec(BNF, 3); R5 = idealprimedec(BNF, 5)[1]; R13 = idealprimedec(BNF, 13); R31 = idealprimedec(BNF, 31); R7 = idealprimedec(BNF, 7)[1];
R2 = idealprimedec(BNF, 2)[1];
print("Q(sqrt10): cl=", BNF.cyc, " classes 3:", [bnfisprincipal(BNF, R3[1], 0), bnfisprincipal(BNF, R3[2], 0)], " 5:", bnfisprincipal(BNF, R5, 0),
      " 13:", [bnfisprincipal(BNF, R13[1], 0)], " 31:", [bnfisprincipal(BNF, R31[1], 0)]);
\\ (i) A = nontrivial class, r = R3[1] (in A^{-1} = A); supports and conductors prime to r
Anon = bnfisprincipal(BNF, R3[1], 0);
plist = [R3[2], R5, R13[1], R13[2], R31[1], R31[2], R7];
fl = [R3[2], idealpow(BNF, R3[2], 2), R5, idealmul(BNF, R3[2], R5), R13[1], idealmul(BNF, R13[1], R3[2]), idealmul(BNF, R2, R3[2])];
betas = [[1,0], [2,0], [3,0], [7,1]];
runconfig("Q(sqrt10) A!=1", plist, Anon, R3[1], 600, fl, betas);
\\ (ii) A = trivial class, r = O_F
plist = [R3[1], R3[2], R5, R13[1], R13[2], R31[1], R7];
fl = [R3[2], idealpow(BNF, R3[2], 2), R5, idealmul(BNF, R3[2], R5), idealmul(BNF, R3[1], R3[2]), idealmul(BNF, R13[1], R3[2])];
runconfig("Q(sqrt10) A=1", plist, 0*Anon, idealhnf(BNF, 1), 500, fl, betas);
}



{
printf("SUMMARY: %d (K,f) tests, %d failed; %d blocks (%d nonzero, %d nonzero with e != 1); %d off-diagonal pairs (all must be 0); %d step-(b) comparisons; ring class number formula %d/%d\n",
  NTEST, NBAD, NBLK, NBLKNZ, NBLKNZE, NOFF, NSTEPB, NGOK, NG);
}
\\ note: NBAD counts tests where any block, the off-diagonal check or step (b) failed


