\\ gl2gapaudit: zero-counting normalisation check for Delta on [T, 2T]: sign changes of Z(t) = L(6+it) e^{i phi(t)},
\\ phi(t) = Im lngamma(6+it) - t log(2 pi) (Lambda(s) = (2pi)^{-s} Gamma(s) L(s) is real on Re s = 6), vs N(2T) - N(T) = (T/pi) log(2T/(pi e)).
\\ run: gp -q -s 64000000 zeros.gp (T set below)
T = 200; h = 1/50;
Li = lfuninit(lfunetaquo([1,24]), [6, 0, 2*T + 1]);
K = T/h + 1;
Z = vector(K, k, my(t = T + (k-1)*h); real(lfun(Li, 6 + I*t) * exp(I*(imag(lngamma(6 + I*t)) - t*log(2*Pi)))));
mx = vector(K, k, my(t = T + (k-1)*h); abs(imag(lfun(Li, 6 + I*t) * exp(I*(imag(lngamma(6 + I*t)) - t*log(2*Pi)))))); 
zs = List(); for (k = 1, K - 1, if (Z[k] == 0 || Z[k]*Z[k+1] < 0, listput(zs, T + (k-1)*h - Z[k]*h/(Z[k+1]-Z[k]))));
n = #zs; pred = T/Pi*log(2*T/(Pi*exp(1)));
printf("T = %d: max |Im Z| = %.2e (phase check); sign changes on [T,2T] = %d; (T/pi)log(2T/(pi e)) = %.2f\n", T, vecmax(mx), n, pred);
gmax = 0; for (i = 1, n - 1, g = (zs[i+1] - zs[i]) * log(zs[i]) / Pi; if (g > gmax, gmax = g; gi = zs[i]));
gmax2 = 0; for (i = 1, n - 1, g = (zs[i+1] - zs[i]) * log(zs[i]/(2*Pi)) / Pi; if (g > gmax2, gmax2 = g));
printf("largest gap / (pi/log gamma) = %.4f (at gamma = %.3f); with log(gamma/2pi): %.4f; mean gap * log(T)/pi = %.4f\n", gmax, gi, gmax2, (zs[n]-zs[1])/(n-1)*log(1.5*T)/Pi);
quit;
