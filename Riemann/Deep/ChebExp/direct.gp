\\ Direct range: every number field L with |disc L| <= D0 = 257, every subfield K with L/K Galois, every class C.
\\ Minkowski: degree >= 5 forces |disc| >= 258 (degree 5, r2 = 2: ((pi/4)^2 5^5/5!)^2 = 258.0), so only
\\ degrees 2, 3, 4 occur; PARI nflist lists them completely by Galois group.
\\ Output: one line per (L, K, class): G deg(L) Delta_L K class least-p
D0 = 257;
groups = ["C2", "C3", "S3", "C4", "V4", "D4", "A4", "S4"];

lp1(nf) =
{
  forprime (p = 2, 10^6,
    my(dec = idealprimedec(nf, p));
    for (i = 1, #dec, if (dec[i].f == 1, return(p))));
  error("lp1");
}

galQ(G, n, DL, nf) =
{
  my(gal = galoisinit(nf), found = Map(), cnt = 0, fr);
  forprime (p = 2, 10^6,
    if (DL % p == 0, next);
    fr = idealfrobenius(nf, gal, idealprimedec(nf, p)[1]);
    if (!mapisdefined(found, fr),
      mapput(found, fr, p); cnt++;
      print(G, " ", n, " ", DL, " K=Q elt", Vec(fr), " ", p));
    if (cnt == n, return(1)));
  error("galQ");
}

quadsub(G, n, DL, nf, f) =
{
  my(S = nfsubfields(nf, 2), k, P, rnf, psplit, pinert, dec, pr, D);
  for (s = 1, #S,
    k = nfinit(subst(S[s][1], x, y)); P = nffactor(k, f)[1, 1]; rnf = rnfinit(k, P);
    psplit = 0; pinert = 0;
    forprime (p = 2, 10^6,
      dec = idealprimedec(k, p);
      for (j = 1, #dec,
        pr = dec[j];
        if (pr.f == 1,
          D = rnfidealprimedec(rnf, pr);
          if (#D == 2 && !psplit, psplit = p);
          if (#D == 1 && D[1].e == 1 && !pinert, pinert = p)));
      if (psplit && pinert, break));
    print(G, " ", n, " ", DL, " K=quad[", S[s][1], "] split ", psplit);
    print(G, " ", n, " ", DL, " K=quad[", S[s][1], "] inert ", pinert));
}

{
  for (gi = 1, #groups,
    G = groups[gi]; V = nflist(G, [1, D0]);
    for (i = 1, #V,
      f = V[i]; nf = nfinit(f); DL = abs(nf.disc); n = poldegree(f);
      print(G, " ", n, " ", DL, " K=L id ", lp1(nf));
      if (G == "C2" || G == "C3" || G == "C4" || G == "V4", galQ(G, n, DL, nf));
      if (n == 4, quadsub(G, n, DL, nf, f))));
}
