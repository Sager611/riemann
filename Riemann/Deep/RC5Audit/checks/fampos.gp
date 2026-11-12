\\ fampos.gp (rc5audit): D > 0 analogue of fam.c. Family D in (X/2, X] fundamental, narrow class group Cl+(K) via bnrinit(bnf,[1,[1,1]]).
\\ LHS = sum_D h+ sum_C |S_C|^2 (= sum_D sum_{chi in Cl+^} |sum_n b(n) lambda_chi(n)|^2), TRIV = sum_D |sum_C S_C|^2.
\\ b on admissible n (squarefree, prime factors >= 19); detector weights as fam.c (sigma, V = sqrt(X)/2, Y = X, kill = ell-term).
default(parisize, 120000000);
sstep(t) = if(t <= 0, 0., if(t >= 1, 1., my(a = exp(-1/t), c = exp(-1/(1-t))); a/(a+c)));
phi2(v) = if(v <= 0.25 || v >= 4, 0., my(l = log(v)/log(2)); if(l < -1, sstep(l+2), if(l > 1, sstep(2-l), 1.)));
run(X, Nl, sig) = {
  my(V = max(2, floor(sqrt(X)/2)), nmax = 4*vecmax(Nl), Y = X, l2 = 4.^(1-sig), adm = List(), K = 2*#Nl, B, NB = vector(K), LHS = vector(K), TRIV = vector(K), fam = 0, famR = 0, nD = 0);
  for (n = 1, nmax, if (issquarefree(n) && gcd(n, 510510) == 1, listput(adm, n)));
  adm = Vec(adm);
  my(pf = vector(#adm, i, factor(adm[i])[,1]~));
  B = matrix(K, #adm);
  for (k = 1, #Nl, my(N = Nl[k]);
    for (i = 1, #adm, my(n = adm[i], v = n/N, a = sumdiv(n, d, if(d <= V, moebius(d), 0)), pre = a * n^(-sig), w1 = exp(-n/Y)*phi2(v), w2 = l2*exp(-4.*n/Y)*phi2(4*v));
      B[2*k-1, i] = pre*(w1 - w2); B[2*k, i] = pre*w1));
  for (k = 1, K, NB[k] = sum(i = 1, #adm, B[k,i]^2));
  forstep (D = floor(X/2) + 1, X, 1, if (!isfundamental(D) || issquare(D), next);
    my(bnf = bnfinit(quadpoly(D), 1), bnr = bnrinit(bnf, [1, [1,1]]), cyc = bnr.cyc, hp = if(#cyc, vecprod(cyc), 1), rad, dl = Map(), S, hR);
    nD++; fam += hp;
    hR = 2 * bnf.no * bnf.reg; famR += hR;
    rad = vector(#cyc); my(r = 1); for (j = 1, #cyc, rad[j] = r; r *= cyc[j]);
    S = matrix(K, hp);
    for (i = 1, #adm, my(ps = pf[i], ok = 1, vs = List());
      for (j = 1, #ps, my(p = ps[j], kr = kronecker(D, p)); if (kr == -1, ok = 0; break);
        my(v); if (!mapisdefined(dl, p, &v), v = if(#cyc, Col(bnrisprincipal(bnr, idealprimedec(bnf, p)[1], 0)), vectorv(0)); mapput(dl, p, v));
        listput(vs, [v, kr]));
      if (!ok, next);
      my(m = #vs, cnt = 1 << m);
      for (s = 0, cnt - 1, my(vv = vectorv(#cyc), skip = 0);
        for (j = 1, m, my(e = if(bittest(s, j-1), -1, 1)); if (vs[j][2] == 0 && e == -1, skip = 1; break); vv += e * vs[j][1]);
        if (skip, next);
        my(idx = 1); for (j = 1, #cyc, idx += (vv[j] % cyc[j]) * rad[j]);
        for (k = 1, K, S[k, idx] += B[k, i])));
    for (k = 1, K, my(t = 0, q = 0); for (c = 1, hp, t += S[k,c]; q += S[k,c]^2); LHS[k] += hp*q; TRIV[k] += t^2));
  printf("# X=%d nD=%d sum_h+=%d sum_h+R+=%.1f X^1.5=%.1f\n", X, nD, fam, famR, X^1.5);
  for (k = 1, #Nl, for (kk = 0, 1, my(j = 2*k-1+kk);
    printf("N=%d kill=%d NB=%.4e LHS=%.4e TRIV=%.4e R_fam=%.4f R_famR=%.5f R_X=%.5f triv%%=%.1f\n", Nl[k], 1-kk, NB[j], LHS[j], TRIV[j], LHS[j]/(fam*NB[j]), LHS[j]/(famR*NB[j]), LHS[j]/(X^1.5*NB[j]), 100*TRIV[j]/LHS[j])));
}
run(1000, [32, 64, 128, 256, 512], 5/6);
run(2000, [32, 64, 128, 256, 512, 1024], 5/6);
run(4000, [64, 128, 256, 512, 1024, 2048], 5/6);
quit;
