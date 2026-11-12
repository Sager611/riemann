default(parisize, 120000000);
default(realprecision, 12);
E = ellinit([0,-1,1,-10,-20]);
Vga=[0,0,1,1];
clgtab(D) = {
  my(q=quadclassunit(D), cyc=q[2], gen=q[3], h=q[1], M=Map(), ng=#cyc, f, one);
  one = Qfb(1, D%2, (D%2-D)/4);
  if (ng==0, mapput(M, Vec(one), []); return([1, [], M]));
  forvec(e = vector(ng, j, [0, cyc[j]-1]),
    f = one;
    for(j=1, ng, f = qfbcomp(f, qfbpow(gen[j], e[j])));
    f = qfbred(f);
    mapput(M, Vec(f), e));
  [h, cyc, M];
}
chival(chi, cyc, e) = if(#cyc==0, 1, exp(2*Pi*I*sum(j=1,#cyc, chi[j]*e[j]/cyc[j])));
cls(M, D, p) = mapget(M, Vec(qfbred(qfbprimeform(D,p))));
locF(D, M, cyc, chi, p, X) = {
  my(k=kronecker(D,p), ap=ellap(E,p), c);
  if (p==11,
    if (k==1, c=chival(chi,cyc,cls(M,D,p)); return((1-ap*c*X)*(1-ap*conj(c)*X)));
    if (k==-1, return(1-ap^2*X^2));
    return(1));
  if (k==1, c=chival(chi,cyc,cls(M,D,p)); return((1-ap*c*X+p*c^2*X^2)*(1-ap*conj(c)*X+p*conj(c)^2*X^2)));
  if (k==-1, return(1-(ap^2-2*p)*X^2+p^2*X^4));
  c=chival(chi,cyc,cls(M,D,p)); return(1-ap*c*X+p*c^2*X^2);
}
DL = [0, 0.05, 0.15]; TL = [0, 1.5];
f = "msq2.out";
Dmax = eval(getenv("DMAX")); Dmin = eval(getenv("DMIN"));
main() = {
for (a = Dmin, Dmax, my(D = -a);
  if (!isfundamental(D) || D % 11 == 0, next);
  my(T = clgtab(D), cyc = T[2], h = T[1], Lmax, sums = vector(6), t0 = getabstime(), feqw = 0);
  Lmax = lfuncost(lfuncreate([1,0,Vga,2,121*D^2,1]), [2])[1];
  my(chis = if(#cyc==0, [[]], my(v=List()); forvec(e = vector(#cyc, j, [0, cyc[j]-1]), listput(v, e)); Vec(v)));
  for (ic = 1, #chis, my(chi = chis[ic], an, L, Li, k = 0);
    an = direuler(p=2, Lmax, 1/locF(D, T[3], cyc, chi, p, X));
    L = lfuncreate([an, 0, Vga, 2, 121*D^2, -kronecker(D,11)]);
    if (ic == 2 || #chis == 1, feqw = lfuncheckfeq(L));
    Li = lfuninit(L, [1.1, 0.1, 2]);
    for (i = 1, 3, for (j = 1, 2, k++; sums[k] += abs(lfun(Li, 1 + DL[i] + I*TL[j]))^2)));
  write(f, D, " ", h, " ", Lmax, " ", feqw, " ", sums, " ", getabstime()-t0));
}
main();
