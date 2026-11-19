\\ chebaudit: numerical check of the Hadamard identity (3.2) of ChebExp DERIVATION for K = Q(sqrt d)
\\ zeros of zeta_K = zeros of zeta and of L(s, chi_d), taken separately (PARI misses close pairs on the product)
default(parisize, "128M"); default(realprecision, 38);
T = 300; Z1 = lfuncreate(1); zz = lfunzeros(Z1, T);
chk(d, s) =
{
  my(Ld = lfuncreate(d), zl = lfunzeros(Ld, T), q = abs(d), sumz, tail, rhs, r1, r2, gam, zp);
  sumz = sum(i = 1, #zz, 2*(s - 1/2)/((s - 1/2)^2 + zz[i]^2)) + sum(i = 1, #zl, 2*(s - 1/2)/((s - 1/2)^2 + zl[i]^2));
  tail = (s - 1/2)/Pi*((log(T/(2*Pi)) + 1)/T + (log(q*T/(2*Pi)) + 1)/T);
  if (d < 0, r1 = 0; r2 = 1, r1 = 2; r2 = 0);
  gam = -r1*(log(Pi) - psi(s/2))/2 - r2*(log(2*Pi) - psi(s));
  zp = lfun(Z1, s, 1)/zeta(s) + lfun(Ld, s, 1)/lfun(Ld, s);
  rhs = log(q)/2 + 1/s + 1/(s - 1) + gam + zp;
  printf("d=%d s=%s zeros<=T: %d+%d  zero side=%.8f  rhs (3.2)=%.8f  diff=%.2e\n", d, s, #zz, #zl, sumz + tail, rhs, sumz + tail - rhs);
}
chk(-3, 2); chk(-3, 1.3); chk(5, 2); chk(5, 1.3); chk(-4, 1.5); chk(12, 1.5);
