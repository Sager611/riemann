default(parisize, "120M");
\\ PARI: qfbclassno(-d, 1) for d in file D (only |D| <= 10^9)
v = readvec(getenv("DF")); for(i = 1, #v, my(d = v[i]); print(d, " ", qfbclassno(-d, 1)));
