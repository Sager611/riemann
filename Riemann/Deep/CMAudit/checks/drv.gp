\\ drv.gp (cmaudit): production run of cmls.gp. usage: gp ... drv.gp with OUT, FIELDS set below by the caller file.
read("Riemann/Deep/CMAudit/checks/cmls.gp");
findalpha(Fpol, T, mode) = {
  my(nfF = nfinit(Fpol), r = real(polroots(Fpol)));
  for (a = ceil(sqrt(T)), 10^7, for (b = 1, 3,
    my(al = a + b*y, N = abs(norm(Mod(al, Fpol))));
    if (a + b*r[1] <= 0 || a + b*r[2] <= 0 || N < T, next);
    if (mode == "p" && !isprime(N), next);
    if (mode == "c" && (isprime(N) || !issquarefree(N) || omega(N) != 2 || N % 2 == 0), next);
    if (mode == "e" && (!isprime(N)), next);
    my(rd = rnfdisc(nfF, x^2 + al), nd = idealnorm(nfF, rd[1]));
    if ((mode == "p" || mode == "c") && nd == N, return(al));
    if (mode == "e" && nd > N && nd % 2 == 0, return(al))));
}
f3(v) = Strprintf("%.3g", v);
runfield(OUT, Fpol, T, mode, fac) = {
  my(t0 = getabstime(), al = findalpha(Fpol, T, mode), K = mkK(Fpol, al), bnf = K[1], h = bnf.no, kap);
  kap = (2*Pi)^2 * h * bnf.reg / (bnf.tu[1] * sqrt(K[7]));
  my(L1 = lemma1(K, 25), R = lemma23(K, fac));
  write(OUT, "FIELD F=", Fpol, " alpha=", al, " mode=", mode, " dK=", K[7], " D=", K[8], " h=", h, " cyc=", bnf.cyc, " kappa=", f3(kap),
    " jsize=", R[5], " nideals=", R[2], " L1[ok,bad,minNDelta/(Nc^2D),maxIdErr,nprincipal]=", [L1[1], L1[2], L1[3], f3(L1[4]), L1[5]], " ms=", getabstime() - t0);
  for (i = 1, #R[1], my(e = R[1][i]);
    write(OUT, "  N/sqrtdK=", f3(e[2]), " N=", e[1], " K05/(N+sqrtdK)=", f3(e[3]), " K01=", f3(e[4]), " maxpiC=", e[5], " maxpiC/Bnew=", f3(e[6]),
      " U/(N+sqrtdK)=", f3(e[7]), " Zprim/(N+sqrtdK)=", f3(e[8]), " max/mean=", f3(e[9]), " KZ/(N+sqrtdK)=", f3(e[10]), " K05/((N+sqrtdK)kappa)=", f3(e[3]/kap)));
}
