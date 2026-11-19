\\ chebaudit: independent witnesses for ChebExp's direct range B (|disc L| <= 257), via factormod and
\\ cyclotomic cosets only (no idealprimedec / idealfrobenius / rnfidealprimedec). Field list = nflist (the stated input).
\\ Output lines: G n DL kind my_witness bound_ok
default(parisize, "64M");
d2 = 1 - (777 - sqrt(1081))/852;
bnd(DL) = (0.62*log(3*DL))^(1/d2);
nlin(f, p) = my(F = factormod(f, p)); sum(i = 1, #F~, poldegree(F[i,1]) == 1);
nquad(f, p) = my(F = factormod(f, p)); sum(i = 1, #F~, poldegree(F[i,1]) == 2);
nfac(f, p) = #factormod(f, p)~;
w1(f) = my(D = poldisc(f)); forprime(p = 2, 10^6, if (D % p && nlin(f, p) > 0, return(p))); error("w1");
COND = Map([49, 7; 81, 9; 169, 13; 125, 5; 144, 12; 225, 15; 256, 8]);
abel(G, n, DL, f) =
{
  my(m = mapget(COND, DL), D = poldisc(f), H = List(), U = [a | a <- [1..m-1], gcd(a, m) == 1], key, seen = Map(), cnt = 0, out = List());
  foreach (U, a, forprime (p = 2, 10^6, if (p % m == a && D % p, if (nfac(f, p) == n, listput(H, a)); break)));
  H = Vec(H);
  key = (a -> vecmin([a*h % m | h <- H]));
  forprime (p = 2, 10^6, if (m % p, my(k = key(p % m)); if (!mapisdefined(seen, k), mapput(seen, k, p); cnt++; listput(out, p)); if (cnt == #U/#H, break)));
  if (#U/#H != n, error("coset count ", DL));
  \\ consistency: number of factors of f mod p = n / order of p's coset, for p !| D up to 3000
  forprime (p = 2, 3000, if (D % p && m % p, my(o = 1, a = p % m); while (setsearch(Set(H), lift(Mod(a, m)^o)) == 0, o++); if (nfac(f, p) != n/o, error("pattern ", DL, " ", p))));
  foreach (out, p, print(G, " ", n, " ", DL, " K=Q coset ", p, " ", p <= bnd(DL)));
}
quad(G, n, DL, nf, f) =
{
  my(S = nfsubfields(nf, 2), D = poldisc(f), dk, ps, pi);
  for (s = 1, #S,
    dk = nfdisc(S[s][1]); ps = 0; pi = 0;
    forprime (p = 2, 10^6, if (D % p && kronecker(dk, p) == 1,
      if (!ps && nlin(f, p) > 0, ps = p); if (!pi && nquad(f, p) > 0, pi = p); if (ps && pi, break)));
    print(G, " ", n, " ", DL, " K=quad[", dk, "] split ", ps, " ", ps <= bnd(DL));
    print(G, " ", n, " ", DL, " K=quad[", dk, "] inert ", pi, " ", pi <= bnd(DL)));
}
{
  foreach (["C3", "S3", "C4", "V4", "D4", "A4", "S4"], G,
    V = nflist(G, [1, 257]);
    for (i = 1, #V,
      f = V[i]; nf = nfinit(f); DL = abs(nf.disc); n = poldegree(f); w = w1(f);
      print(G, " ", n, " ", DL, " K=L id ", w, " ", w <= bnd(DL));
      if (G == "C3" || G == "C4" || G == "V4", abel(G, n, DL, f));
      if (n == 4, quad(G, n, DL, nf, f))));
}
