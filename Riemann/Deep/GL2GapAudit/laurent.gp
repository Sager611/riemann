\\ gl2gapaudit: Laurent coefficients of D(s) = sum tau(n)^2 n^(-11) n^(-s) at s = 1 (normalised), Delta, by smoothed sums
default(parisize, 120000000);
N = 160000;
t0 = getabstime();
V = mfcoefs(mfDelta(), N);  \\ V[n+1] = tau(n)
print("tau ms ", getabstime()-t0, "  tau(2) = ", V[3], " tau(3) = ", V[4]);
a2 = vector(N, n, V[n+1]^2 / n^11 * 1.);
kill(V);
w = 'w;
lauG(X) = exp(lngamma(1+w+O(w^6)) + w*log(X))/w;  \\ Gamma(w) X^w
res = vector(3);
{
for (i = 1, 3, X = 1000*2^(i-1);
  S0 = 0.; S1 = 0.; S2 = 0.; S3 = 0.;
  for (n = 1, N, e = a2[n]*exp(-n/X); S0 += e; e /= n; S1 += e; ln = log(n); S2 += e*ln; S3 += e*ln^2);
  G = lauG(X);
  res[i] = [X, S0, S1, S2, S3, polcoef(G,0), polcoef(G,1), polcoef(G,2)];
  print("X=", X, " S0/X=", S0/X, " S1=", S1, " S2=", S2, " S3=", S3);
);
}
\\ c from S1 differences (S1 = c A0 + d + O(X^-3/4) tiny)
c12 = (res[2][3]-res[1][3])/(res[2][6]-res[1][6]); c23 = (res[3][3]-res[2][3])/(res[3][6]-res[2][6]);
print("c (1000->2000) = ", c12, "   c (2000->4000) = ", c23);
c = c23;
for (i = 1, 3, r = res[i]; d = r[3] - c*r[6]; e = c*r[7] - r[4]; f = (r[5] - 2*c*r[8])/2; print("X=", r[1], "  d=", d, "  e=", e, "  f=", f));
\\ cross-check c = L(1, sym^2 Delta)/zeta(2) via Petersson norm route is not used; S0/X -> c as a rough check
