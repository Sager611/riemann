\\ Lemma frob check (Banks 2610.31402): theta^{Np} = sigma(theta) mod Q_{p,sigma},
\\ N(Q_{p,sigma}) = Np^{|G|/|C_p|}.  Case 1: k = Q, K = splitting field of x^3-2 (G = S_3).
\\ Case 2: k = Q(sqrt-3), K the same field, G = Gal(K/k) = C_3, Np = p or p^2.
default(parisize, 64000000);
setrand(1);
P6 = subst(polredabs(polcompositum(x^3-2, x^2+3)[1]), x, y);
K = nfinit(P6); Gal = galoisinit(K);
n = poldegree(P6);
clsize(g) = #Set([ g2^-1*g*g2 | g2 <- Gal.group ]);
randelt() = vector(n, i, random(41) - 20)~;
{
bad = 0; checked = 0; normbad = 0; nprimes = 0;
forprime(p = 5, 400,
  nprimes++;
  dec = idealprimedec(K, p);
  frobs = vector(#dec, i, idealfrobenius(K, Gal, dec[i]));
  for (i = 1, #dec,
    sig = galoispermtopol(Gal, frobs[i]);
    for (t = 1, 6,
      th = randelt();
      d = nfeltpow(K, th, p) - nfgaloisapply(K, sig, th);
      if (d != 0 && nfeltval(K, d, dec[i]) < 1, bad++);
      checked++));
  fs = Set(frobs);
  for (j = 1, #fs,
    nrm = prod(i = 1, #dec, if (frobs[i] == fs[j], idealnorm(K, dec[i]), 1));
    c = clsize(fs[j]);
    if (nrm != p^(6/c), normbad++; print("norm mismatch p=", p, " ", nrm, " vs ", p^(6/c)))));
print("case1 S3 over Q: primes=", nprimes, " congruence checks=", checked, " failures=", bad, " norm-formula failures=", normbad);
}
{
r = nfroots(K, x^2 + x + 1)[1];
bad2 = 0; chk2 = 0; notrel = 0; normbad2 = 0;
forprime(p = 5, 400,
  f = if (p % 3 == 1, 1, 2);
  Np = p^f;
  dec = idealprimedec(K, p);
  for (i = 1, #dec,
    phi = idealfrobenius(K, Gal, dec[i]);
    sig = galoispermtopol(Gal, phi^f);
    if (nfgaloisapply(K, sig, r) != r, notrel++);
    for (t = 1, 6,
      th = randelt();
      d = nfeltpow(K, th, Np) - nfgaloisapply(K, sig, th);
      if (d != 0 && nfeltval(K, d, dec[i]) < 1, bad2++);
      chk2++));
  \\ G = C_3 abelian: |C_p| = 1, Q_{p,sigma} = all primes of K above the prime of k: norm Np^3 per prime of k
  nk = if (f == 1, 2, 1);
  if (prod(i = 1, #dec, idealnorm(K, dec[i])) != Np^(3 * nk), normbad2++));
print("case2 C3 over Q(sqrt-3): congruence checks=", chk2, " failures=", bad2, " Frobenius not fixing k: ", notrel, " norm failures: ", normbad2);
}
quit;

