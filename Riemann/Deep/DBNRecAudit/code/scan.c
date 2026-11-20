// dbnrecaudit scan.c: independent anchor scan (falsification only). Ranks every anchor k of a shell by
//   P(k) = max over a few in-anchor offsets of sum over a few groups g of G_g(phase_g(k, offset)),  G_g tabulated on 4096 phase bins (centre values),
// phases as uint64 turns p0 + k dK + off (own mpmath data, scan_wrap.py). Keeps the top NTOP anchors. Not a bound: a search for alignments.
// usage: scan FILE NTOP  (FILE: NGr K NOFF, NOFF base values (a group held constant over anchors, e.g. b = 2), then per group: p0 dK, NOFF offsets (uint64), 4096 table values)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); int ntop = atoi(argv[2]);
  int NGr, NOFF; long long K; if (fscanf(f, "%d %lld %d", &NGr, &K, &NOFF) != 3) return 1;
  uint64_t *p0 = malloc(NGr*8), *dK = malloc(NGr*8), *off = malloc(NGr*NOFF*8); double *tab = malloc(NGr*4096*sizeof(double)), *base = malloc(NOFF*sizeof(double));
  for (int o = 0; o < NOFF; o++) if (fscanf(f, "%lf", &base[o]) != 1) return 5;
  for (int g = 0; g < NGr; g++) {
    unsigned long long a, b; if (fscanf(f, "%llu %llu", &a, &b) != 2) return 2; p0[g] = a; dK[g] = b;
    for (int o = 0; o < NOFF; o++) { unsigned long long c; if (fscanf(f, "%llu", &c) != 1) return 3; off[g*NOFF + o] = c; }
    for (int i = 0; i < 4096; i++) if (fscanf(f, "%lf", &tab[g*4096 + i]) != 1) return 4;
  }
  double *hv = malloc(ntop*sizeof(double)); long long *hk = malloc(ntop*8); int *ho = malloc(ntop*sizeof(int)); int nh = 0;
  for (int i = 0; i < ntop; i++) hv[i] = -1e300;
  double hmin = -1e300; int imin = 0;
  for (long long k = 0; k < K; k++) {
    double best = -1e300; int bo = 0;
    for (int o = 0; o < NOFF; o++) {
      double s = base[o];
      for (int g = 0; g < NGr; g++) { uint64_t ph = p0[g] + (uint64_t)k*dK[g] + off[g*NOFF + o]; s += tab[g*4096 + (ph >> 52)]; }
      if (s > best) { best = s; bo = o; }
    }
    if (best > hmin) {           // replace the current minimum of the top list
      hv[imin] = best; hk[imin] = k; ho[imin] = bo; if (nh < ntop) nh++;
      hmin = 1e300; for (int i = 0; i < ntop; i++) if (hv[i] < hmin) { hmin = hv[i]; imin = i; }
    }
  }
  for (int i = 0; i < ntop; i++) if (hv[i] > -1e299) printf("%lld %d %.9f\n", hk[i], ho[i], hv[i]);
  return 0;
}

