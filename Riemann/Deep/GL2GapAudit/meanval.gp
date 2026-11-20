\\ gl2gapaudit: numerical test of Thm mixed (1410.7765) for Delta (level 1, weight 12), mu, nu <= 2, on [T, 2T].
\\ usage: gp -q -D T=500 meanval.gp ; L-values by lfun (PARI), derivatives by 8th-order central differences (6th-order as error check).
\\ run as: gp -q -s 120000000 <wrapper defining T that reads this file>
if (type(T) != "t_INT", T = 500);
h = 1/50;
\\ ---- Laurent coefficients of D(1+z) = c/z + d + e z + f z^2 (D = sum |a(n)|^2 n^-s), smoothed sums with 1/X corrections
N = 160000; V = mfcoefs(mfDelta(), N); a2 = vector(N, n, V[n+1]^2 / n^11 * 1.); kill(V);
w = 'w; lauG(X) = exp(lngamma(1+w+O(w^6)) + w*log(X))/w;
Xs = [1000, 2000, 4000]; S = matrix(3, 4); A = matrix(3, 3);
{
for (i = 1, 3, X = Xs[i]; s0 = 0.; s1 = 0.; s2 = 0.; s3 = 0.;
  for (n = 1, N, e = a2[n]*exp(-n/X); s0 += e; e /= n; s1 += e; ln = log(n); s2 += e*ln; s3 += e*ln^2);
  S[i,] = [s0, s1, s2, s3]; G = lauG(X); A[i,] = [polcoef(G,0), polcoef(G,1), polcoef(G,2)]);
}
\\ S0/X = c + D0/X + D1/X^2 (3 equations)
M3 = matrix(3, 3, i, j, Xs[i]^(1-j)); sol = matsolve(M3, vector(3, i, S[i,1]/Xs[i])~); c = sol[1];
\\ S1 = c A0 + d + k/X ; S2 = c A1 - e + k/X ; S3 = 2 c A2 + 2 f + k/X : two-point fits (X = 2000, 4000), X = 1000 as check
fit(y1, y2) = (2*y2 - y1);  \\ value at 1/X -> 0 from X and 2X
d = fit(S[2,2] - c*A[2,1], S[3,2] - c*A[3,1]); d1 = fit(S[1,2] - c*A[1,1], S[2,2] - c*A[2,1]);
e = -fit(S[2,3] - c*A[2,2], S[3,3] - c*A[3,2]); e1 = -fit(S[1,3] - c*A[1,2], S[2,3] - c*A[2,2]);
f = fit(S[2,4] - 2*c*A[2,3], S[3,4] - 2*c*A[3,3])/2; f1 = fit(S[1,4] - 2*c*A[1,3], S[2,4] - 2*c*A[2,3])/2;
printf("Laurent: c = %.9f (D(0) = %.5f), d = %.7f [%.7f], e = %.6f [%.6f], f = %.5f [%.5f]  ([.] = fit from X = 1000, 2000)\n", c, sol[2], d, d1, e, e1, f, f1);
kill(a2);
\\ ---- grid of L(6 + i t), t in [T - 5h, 2T + 5h]
Li = lfuninit(lfunetaquo([1,24]), [6, 0, 2*T + 1]);
n0 = 5; K = T/h + 2*n0 + 1; t1 = T - n0*h;
t0 = getabstime(); Lv = vector(K, k, lfun(Li, 6 + I*(t1 + (k-1)*h))); printf("grid: %d values in %.1f s\n", K, (getabstime()-t0)/1000.);
c1_8 = [1/280, -4/105, 1/5, -4/5, 0, 4/5, -1/5, 4/105, -1/280]; c2_8 = [-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560];
c1_6 = [0, -1/60, 3/20, -3/4, 0, 3/4, -3/20, 1/60, 0]; c2_6 = [0, 1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90, 0];
fd(cf, k, p) = sum(j = 1, 9, cf[j]*Lv[k + j - 5]) / h^p;
\\ L^(mu)(s): L' = -i dL/dt, L'' = -d^2L/dt^2
Iint(o) = {
  my(I3 = matrix(3, 3), wt, D);
  for (k = n0 + 1, K - n0, wt = if (k == n0 + 1 || k == K - n0, h/2, h);
    D = [Lv[k], -I*fd(if (o == 8, c1_8, c1_6), k, 1), -fd(if (o == 8, c2_8, c2_6), k, 2)];
    for (m = 1, 3, for (n = 1, 3, I3[m, n] += wt * D[m] * conj(D[n]))));
  I3 };
I8 = Iint(8); I6 = Iint(6);
\\ ---- predictions
LT = log(T);
Plead(n) = (-1)^n * 2^(n+1)/(n+1) * c * T * LT^(n+1);
m(n, l) = c*(-1)^n*(2*l)^(n+1)/(n+1) + d*if (n == 0, 2, (-2*l)^n) + if (n >= 2, -e*n*(-2*l)^(n-1) + f*(2*(n == 2) + n*(n-1)*(-2*l)^(n-2)), 0);
Prec(n) = intnum(t = T, 2*T, m(n, log(t/(2*Pi))));
mf(n, l) = c*(-1)^n*(2*l)^(n+1)/(n+1) + d*if (n == 0, 2, (-2*l)^n);   \\ without e, f: size of the truncation
Precd(n) = intnum(t = T, 2*T, mf(n, log(t/(2*Pi))));
printf("T = %d, h = %s, log T = %.4f\n", T, h, LT);
{
for (mu = 0, 2, for (nu = mu, 2, n = mu + nu; v = I8[mu+1, nu+1];
  printf("mu=%d nu=%d: I = %.6f %+.6f i | FD8-FD6 %.1e | I/Plead = %.5f | I/Prec = %.6f | (c,d only: %.6f)\n", mu, nu, real(v), imag(v),
    abs(v - I6[mu+1, nu+1])/abs(v), real(v)/Plead(n), real(v)/Prec(n), real(v)/Precd(n))));
}
printf("check Im I(1,0) = -(1/2)[|L|^2]: %.6f vs %.6f\n", imag(I8[2,1]), -(abs(Lv[K-n0])^2 - abs(Lv[n0+1])^2)/2);
\\ g = e^{it log T} L: G_j = sum_{k,l<=j} C(j,k)C(j,l) LT^{2j-k-l} I(k,l); leading 2 c T LT^{2j+1}/(2j+1)
{
for (j = 0, 2, Gj = sum(k = 0, j, sum(l = 0, j, binomial(j,k)*binomial(j,l)*LT^(2*j-k-l)*I8[k+1,l+1]));
  Pj = sum(k = 0, j, sum(l = 0, j, binomial(j,k)*binomial(j,l)*LT^(2*j-k-l)*Prec(k+l)));
  if (j == 0, G0 = real(Gj));
  printf("j=%d: int|g^(j)|^2 = %.4f (Im %.1e) | / (2cT L^{2j+1}/(2j+1)) = %.5f | / recipe = %.6f | ratio G_j/(G_0 L^{2j}) = %.5f vs 1/(2j+1) = %.5f\n",
    j, real(Gj), imag(Gj), real(Gj)/(2*c*T*LT^(2*j+1)/(2*j+1)), real(Gj)/real(Pj), real(Gj)/(G0*LT^(2*j)), 1/(2*j+1.)));
}
