\\ PARI quadclassunit(-d).no (GRH-conditional; consistency check of hcount at large d only)
default(parisize, "64M");
v = readvec(getenv("DF")); for(i = 1, #v, my(d = v[i]); print(d, " ", quadclassunit(-d).no));
