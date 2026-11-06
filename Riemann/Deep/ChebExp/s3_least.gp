\\ least primes in the three classes of Gal(L/Q) = S_3, L = Galois closure of a cubic field K3
\\ output lines: disc(K3) Delta_L p_id p_transposition p_3cycle
V = concat(nflist("S3", [1, 3000], 0), nflist("S3", [1, 3000], 1));
{
for (i = 1, #V,
  f = V[i]; K = nfinit(f); d = K.disc; DL = abs(nfdisc(nfsplitting(f)));
  pid = 0; ptr = 0; pc = 0;
  forprime (p = 2, 10^7,
    if (d % p == 0, next);
    dec = idealprimedec(K, p); fs = vecsort(vector(#dec, j, dec[j].f));
    if (fs == [1,1,1] && !pid, pid = p);
    if (fs == [1,2] && !ptr, ptr = p);
    if (fs == [3] && !pc, pc = p);
    if (pid && ptr && pc, break));
  print(d, " ", DL, " ", pid, " ", ptr, " ", pc));
}
