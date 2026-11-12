\\ genus.gp (cmaudit): CMQuartic Lemma 4 count check: |Cl_K/(1-c)Cl_K| (= #genus characters chi^c = chi) vs |G(K)|/2,
\\ G(K) = characters of Cl_F(d oo) with order <= 2 on each O_p^x (p | d), computed independently from F's narrow ray class group.
read("Riemann/Deep/CMAudit/checks/cmls.gp");
findalpha2(Fpol, T, mode) = {
  my(nfF = nfinit(Fpol), r = real(polroots(Fpol)));
  for (a = ceil(sqrt(T)), ceil(sqrt(T)) + 3000, for (b = 1, 5,
    my(al = a + b*y, N = abs(norm(Mod(al, Fpol))));
    if (a + b*r[1] <= 0 || a + b*r[2] <= 0 || N < T, next);
    if (mode == "p" && !isprime(N), next);
    if (mode == "c" && (isprime(N) || !issquarefree(N) || omega(N) != 2 || N % 2 == 0), next);
    my(rd = rnfdisc(nfF, x^2 + al), dd = rd[1], fa = idealfactor(nfF, dd), ok = 1);
    for (i = 1, #fa~, if (fa[i,1].f != 1 || fa[i,2] != 1 || fa[i,1].p == 2, ok = 0));
    if (ok && idealnorm(nfF, dd) == N, return([al, dd]))));
  error("no alpha");
}
genuscheck(Fpol, T, mode) = {
  my(R = findalpha2(Fpol, T, mode), al = R[1], dd = R[2], K = mkK(Fpol, al), bnf = K[1], nf = K[2], c = K[3], cyc = bnf.cyc, r = #cyc, lhs, rhs);
  if (r == 0, lhs = 1,
    my(M = matrix(r, r));
    for (i = 1, r, my(v = bnfisprincipal(bnf, nfgaloisapply(nf, c, bnf.gen[i]), 0)); v[i] -= 1; M[, i] = v);
    lhs = abs(matdet(mathnf(concat(M, matdiagonal(cyc))))));
  my(bnfF = bnfinit(Fpol, 1), bnr = bnrinit(bnfF, [dd, [1, 1]]), rc = bnr.cyc, fa = idealfactor(bnfF, dd), H = List(), Nd = idealnorm(bnfF, dd));
  for (i = 1, #fa~, my(pr = fa[i, 1], p = pr.p, g = znprimroot(p), xs = vector(#fa~, j, if(j == i, lift(g), 1)), xx = idealchinese(bnfF, fa, xs));
    xx = nfbasistoalg(bnfF, xx);
    my(k = 0); while (!(subst(lift(xx), y, real(polroots(Fpol))[1]) + k*Nd > 0 && subst(lift(xx), y, real(polroots(Fpol))[2]) + k*Nd > 0), k++);
    xx += k*Nd;
    listput(H, bnrisprincipal(bnr, lift(xx)^2 % Fpol, 0)));
  my(rr = #rc);
  rhs = if (rr == 0, 1, abs(matdet(mathnf(concat(Mat(Vec(H)), matdiagonal(rc))))));
  [Fpol, al, K[7], idealnorm(bnfF, dd), bnf.no, cyc, bnfF.no, bnrinit(bnfF, [1, [1,1]]).no, lhs, rhs, lhs == rhs / 2]
}
OUT = "Riemann/Deep/CMAudit/checks/out_genus.txt";
{
  foreach ([[y^2 - y - 1, [500, 3000, 2*10^4, 10^5]], [y^2 - 2, [500, 5000, 5*10^4]], [y^2 - 3, [500, 5000, 5*10^4]], [y^2 - 10, [500, 5000, 3*10^4]]], FT,
    foreach (FT[2], T, foreach (["p", "c"], md,
      my(res = iferr(genuscheck(FT[1], T, md), E, [FT[1], T, md, "ERR", E]));
      write(OUT, res))));
  write(OUT, "DONE");
}
