\\ Side result (e) where its bound is smallest (loglog q = 3, bound 45): L(1,chi_-d) = pi h/sqrt(d) for fundamental -d
\\ with d >= e^(e^3) and kronecker(-d,p) = +1 for p <= 31 (largest-L-like), vs 15 loglog d.
default(parisize, 100000000);
{
  my(d = ceil(exp(exp(3))) - 1, n = 0, mx = 0, md = 0, h, Lv, ok);
  while (n < 25, d++;
    if (d % 8 != 7, next);
    ok = 1; forprime(p = 3, 31, if (kronecker(-d, p) != 1, ok = 0; break)); if (!ok, next);
    if (!isfundamental(-d), next);
    h = qfbclassno(-d); if (h != quadclassunit(-d).no, print("MISMATCH ", d));
    Lv = Pi*h/sqrt(d); n++; if (Lv > mx, mx = Lv; md = d));
  printf("loglog d >= 3: max L(1,chi_-d) over 25 largest-L-like d = %.4f (d = %d) vs 15 loglog d = %.2f\n", mx, md, 15*log(log(md)));
}
quit;
