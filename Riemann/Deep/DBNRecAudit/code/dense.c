// dbnrecaudit dense.c: float evaluation of a finite cosine sum on a uniform grid (falsification sampling only; nothing certified).
// input file: N, then N lines "w lo hi phi" (phi = phase at delta = 0, radians; w = log(n)/2 = d phi/d delta).
// mode a: S(delta) = sum lo*cos(phi + w delta)            (lo = signed actual amplitude; hi ignored)
// mode b: S(delta) = sum (c >= 0 ? hi*c : lo*c), c = cos(phi + w delta)  (the bound function Q+ with intervals [lo, hi])
// usage: dense FILE mode delta0 delta1 npts   -> prints max, argmax, min, mean
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); char mode = argv[2][0]; double d0 = atof(argv[3]), d1 = atof(argv[4]); long np = atol(argv[5]);
  int N; if (fscanf(f, "%d", &N) != 1) return 1;
  double *w = malloc(N*sizeof(double)), *lo = malloc(N*sizeof(double)), *hi = malloc(N*sizeof(double)), *ph = malloc(N*sizeof(double));
  for (int i = 0; i < N; i++) if (fscanf(f, "%lf %lf %lf %lf", &w[i], &lo[i], &hi[i], &ph[i]) != 4) return 2;
  double mx = -1e300, mn = 1e300, amx = 0, sum = 0;
  for (long k = 0; k < np; k++) {
    double dl = np > 1 ? d0 + (d1 - d0)*(double)k/(double)(np - 1) : d0, s = 0;
    if (mode == 'a') for (int i = 0; i < N; i++) s += lo[i]*cos(ph[i] + w[i]*dl);
    else for (int i = 0; i < N; i++) { double c = cos(ph[i] + w[i]*dl); s += c >= 0 ? hi[i]*c : lo[i]*c; }
    if (s > mx) { mx = s; amx = dl; } if (s < mn) mn = s; sum += s;
  }
  printf("%.12f %.12f %.12f %.12f\n", mx, amx, mn, sum/np);
  return 0;
}

