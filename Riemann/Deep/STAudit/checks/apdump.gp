default(parisize, 64000000);
E = ellinit([0,-1,1,-10,-20]);
f = "ap.txt";
v = List();
forprime(p = 2, 2000000, listput(v, Str(p, " ", ellap(E, p))));
writefile = 0;
s = strjoin(Vec(v), "\n");
write(f, s);
print("done ", #v);
