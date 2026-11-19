\\ Falsification near the Corollary 2 thresholds: fundamental -d just above d0, worst-like (kronecker(-d,p) = -1 for
\\ p = 2..31) and plain consecutive ones; h by qfbclassno and quadclassunit (compared); bound h > H and
\\ h >= 0.0350 sqrt(d)/max(3, loglog d).  Run under memcap with a time limit.
default(parisize, 100000000);
run(d0, H, need, worst) = {
  my(d = ceil(d0) - 1, found = 0, mn = [10^9, 0, 0], mh = [10^18, 0], h, h2, ok, r);
  while (found < need,
    d++;
    if (worst && d % 8 != 3, next);
    if (worst, ok = 1; forprime(p = 3, 31, if (kronecker(-d, p) != -1, ok = 0; break)); if (!ok, next));
    if (!isfundamental(-d), next);
    h = qfbclassno(-d); h2 = quadclassunit(-d).no; found++;
    if (h != h2, print("  MISMATCH d=", d, " qfbclassno=", h, " quadclassunit=", h2));
    r = h * max(3, log(log(d))) / sqrt(d);
    if (r < mn[1], mn = [r, d, h]);
    if (h < mh[1], mh = [h, d]);
    if (h <= H || r < 0.0350, print("  COUNTEREXAMPLE d=", d, " h=", h)));
  printf("d0=%.4g H=%d %s: %d discs, min h = %d (d=%d), min h*max(3,loglog d)/sqrt(d) = %.5f (d=%d, h=%d) vs 0.0350\n", d0, H, if(worst,"worst-like","consecutive"), found, mh[1], mh[2], mn[1], mn[2], mn[3]);
}
{
  my(T = [[6.973e7, 100], [7.292e9, 1000], [2.957e10, 2000], [7.627e11, 10^4], [7.976e13, 10^5], [3.047e10, 2000], [8.137e13, 10^5]]);
  for (i = 1, #T, run(T[i][1], T[i][2], 200, 0); run(T[i][1], T[i][2], 25, 1));
}
quit;
