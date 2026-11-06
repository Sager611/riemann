\\ Littlewood polynomials: exhaustive irreducibility counts (monic, coefficients +-1)
\\ columns: n, 2^n, #reducible, #reducible with f(1)f(-1)=0, #other reducible, frac reducible, n*frac other, ms
lw(n) = {
  my(red = 0, pm = 0, other = 0, f, v, s1, s2);
  for (b = 0, 2^n - 1,
    v = vector(n + 1, i, if (i == n + 1, 1, if (bittest(b, i - 1), 1, -1)));
    s1 = vecsum(v); s2 = sum(i = 1, n + 1, if (i % 2, v[i], -v[i]));
    f = Pol(Vecrev(v));
    if (!polisirreducible(f),
      red++;
      if (s1 == 0 || s2 == 0, pm++, other++)));
  [n, 2^n, red, pm, other];
}
{
  for (n = NMIN, NMAX,
    my(t = getabstime(), r = lw(n));
    printf("%d %d %d %d %d %.6f %.6f %d\n", r[1], r[2], r[3], r[4], r[5], r[3] / r[2], r[5] * n / r[2], getabstime() - t));
}
