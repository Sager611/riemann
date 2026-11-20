f=nfsplitting(x^3-2); print("x^3-2: n_M=", poldegree(f), " Delta_M=", nfdisc(f), " Delta_k=", nfdisc(x^3-2));
g=nfsplitting(x^3-x^2-2*x-20); print("x^3-x^2-2x-20: n_M=", poldegree(g), " Delta_M=", nfdisc(g), " Delta_k=", nfdisc(x^3-x^2-2*x-20));
foreach([x^2-x-1, x^3-x^2-2*x+1], h, print(h, " disc ", nfdisc(h), " galois? ", #nfgaloisconj(h) == poldegree(h)));
quit;
