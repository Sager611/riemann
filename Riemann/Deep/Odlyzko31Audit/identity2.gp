\\ odlaudit: Poitou's explicit formula (P) as an identity, independent PARI/GP implementation.
\\ (P) log|d| = n(g + log 8pi - Ic) + r1(pi/2 - Ir) - E + sum_rho Phi(rho) + 2 sum_{P,m} log NP NP^(-m/2) F(m log NP)
\\ Run: gp -q -f identity.gp > logs/identity.log
default(realprecision, 50);
default(parisizemax, 120000000);
SG(t) = (1 - t/2)*cos(Pi*t/2) + sin(Pi*t/2)/Pi;
SH(t) = (2 - t)*(1 + cos(Pi*t)/2)/3 + sin(Pi*t)/(2*Pi);
\\ kernel parameters set globally: shape sh (0 = G, 1 = H), b, w
Sx(t) = if(sh, SH(t), SG(t));
F(x) = my(y = abs(x)); if(y >= 2*b, 0, Sx(y/b)/cosh(w*y));
GL = intnumgaussinit(60);
\\ Fhat(t) = int F(x) e^{itx} dx = 2 int_0^X F cos(tx), Gauss-Legendre on 24 pieces (F analytic on [0, X])
Fhat(t) = my(X = 2*b, h = X/24); 2*sum(k = 0, 23, intnumgauss(x = k*h, (k+1)*h, Sx(x/b)/cosh(w*x)*cos(t*x), GL));
kconsts() = {
  my(X = 2*b, Ic, Ir, E);
  Ic = intnum(x = 0, X, (1 - F(x))/(2*sinh(x/2))) + intnum(x = X, [+oo, 1/2], 1/(2*sinh(x/2)));
  Ir = intnum(x = 0, X, (1 - F(x))/(2*cosh(x/2))) + intnum(x = X, [+oo, 1/2], 1/(2*cosh(x/2)));
  E = 4*intnum(x = 0, X, F(x)*cosh(x/2));
  [Ic, Ir, E];
}
\\ prime-ideal sum 2 sum log NP NP^{-m/2} F(m log NP); pol = 0 means K = Q
psum(pol) = {
  my(X = 2*b, s = 0., nf, dec, Np);
  if(pol != 0, nf = nfinit(pol));
  forprime(p = 2, floor(exp(X)) + 1,
    dec = if(pol == 0, [1], idealprimedec(nf, p));
    for(j = 1, #dec,
      Np = if(pol == 0, p, p^(dec[j].f));
      my(m = 1); while(m*log(Np) < X, s += log(Np)*Np^(-m/2)*F(m*log(Np)); m++)));
  2*s;
}
check(name, pol, Ldata, absd, n, r1, T) = {
  my(c = kconsts(), Ic = c[1], Ir = c[2], E = c[3], zs, Z1 = 0., Z2 = 0., P, gam, rhs, t0 = getabstime());
  zs = lfunzeros(Ldata, T);
  for(k = 1, #zs, my(v = 2*Fhat(zs[k])); if(zs[k] <= T/2, Z1 += v); Z2 += v);
  P = psum(pol);
  gam = n*(Euler + log(8*Pi) - Ic) + r1*(Pi/2 - Ir);
  rhs = gam - E + Z2 + P;
  printf("%s (n=%d r1=%d |d|=%d) kernel %s b=%s w=%s: Ic=%.12f Ir=%.12f E=%.12f\n", name, n, r1, absd, if(sh, "H", "G"), b, w, Ic, Ir, E);
  printf("   Gamma terms %.12f  -E %.12f  zeros(%d, gamma<=%d) %.12f [half-height %.12f]  primes %.12f\n", gam, -E, #zs, T, Z2, Z1, P);
  printf("   log|d| = %.15f  rhs = %.15f  residual %.3e  (tail change T/2->T: %.2e; %d ms)\n", log(absd), rhs, log(absd) - rhs, Z2 - Z1, getabstime() - t0);
}
\\ odlaudit: rerun of the Q(sqrt-3) identity with the zeros of zeta_K taken as the union of the zeros of zeta and of L(s, chi_-3)
\\ (lfunzeros on the product misses close pairs of zeros from the two factors). Run: gp -q -f identity2.gp > logs/identity2.log
check2(name, pol, Ls, absd, n, r1, T) = {
  my(c = kconsts(), Ic = c[1], Ir = c[2], E = c[3], zs = [], Z1 = 0., Z2 = 0., P, gam, rhs, t0 = getabstime());
  for(j = 1, #Ls, zs = concat(zs, lfunzeros(Ls[j], T, 32)));
  for(k = 1, #zs, my(v = 2*Fhat(zs[k])); if(zs[k] <= T/2, Z1 += v); Z2 += v);
  P = psum(pol);
  gam = n*(Euler + log(8*Pi) - Ic) + r1*(Pi/2 - Ir);
  rhs = gam - E + Z2 + P;
  printf("%s (n=%d r1=%d |d|=%d) kernel %s b=%s w=%s: Ic=%.12f Ir=%.12f E=%.12f\n", name, n, r1, absd, if(sh, "H", "G"), b, w, Ic, Ir, E);
  printf("   Gamma terms %.12f  -E %.12f  zeros(%d, gamma<=%d) %.12f [half-height %.12f]  primes %.12f\n", gam, -E, #zs, T, Z2, Z1, P);
  printf("   log|d| = %.15f  rhs = %.15f  residual %.3e  (tail change T/2->T: %.2e; %d ms)\n", log(absd), rhs, log(absd) - rhs, Z2 - Z1, getabstime() - t0);
}
sh = 1; b = 3/2; w = 15/32;
check2("Q(sqrt-3) [zeta x L(chi_-3)]", x^2 + x + 1, [lfuncreate(1), lfuncreate(-3)], 3, 2, 0, 400);
sh = 0; b = 2; w = 3733819666132076/10000000000000000;
check2("Q(sqrt-3) [zeta x L(chi_-3)]", x^2 + x + 1, [lfuncreate(1), lfuncreate(-3)], 3, 2, 0, 400);
quit;
