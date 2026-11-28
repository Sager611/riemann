// tabucol.c: k-colouring by tabu search (Hertz & de Werra).  Input (stdin, binary): int32 n, int32 m, int32 k, int32 iters,
// int32 seed, n x int32 initial colours, m x (int32 u, int32 v).  Output (stdout, binary): int32 conflicts, n x int32 colours.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned long long rs = 88172645463325252ULL;
static unsigned rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (unsigned)(rs >> 11); }
int main(void) {
  int h[5]; if (fread(h, 4, 5, stdin) != 5) return 1;
  int n = h[0], m = h[1], k = h[2], iters = h[3]; rs ^= (unsigned long long)h[4] * 2654435761ULL;
  int *col = malloc(4 * (size_t)n), *U = malloc(4 * (size_t)m), *V = malloc(4 * (size_t)m);
  if (fread(col, 4, n, stdin) != (size_t)n) return 1;
  for (int e = 0; e < m; e++) { int uv[2]; if (fread(uv, 4, 2, stdin) != 2) return 1; U[e] = uv[0]; V[e] = uv[1]; }
  int *deg = calloc(n + 1, 4), *ptr = malloc(4 * (size_t)(n + 1)), *nbr = malloc(8 * (size_t)m);
  for (int e = 0; e < m; e++) { deg[U[e]]++; deg[V[e]]++; }
  ptr[0] = 0; for (int v = 0; v < n; v++) ptr[v + 1] = ptr[v] + deg[v];
  int *fill = calloc(n, 4);
  for (int e = 0; e < m; e++) { nbr[ptr[U[e]] + fill[U[e]]++] = V[e]; nbr[ptr[V[e]] + fill[V[e]]++] = U[e]; }
  int *gam = calloc((size_t)n * k, 4); long long *tabu = calloc((size_t)n * k, 8);
  for (int v = 0; v < n; v++) for (int e = ptr[v]; e < ptr[v + 1]; e++) gam[(size_t)v * k + col[nbr[e]]]++;
  long long conf = 0; for (int v = 0; v < n; v++) conf += gam[(size_t)v * k + col[v]]; conf /= 2;
  long long best = conf; int *bestcol = malloc(4 * (size_t)n); memcpy(bestcol, col, 4 * (size_t)n);
  for (long long it = 0; it < iters && conf > 0; it++) {
    int bd = 1 << 30, bv = -1, bc = -1, cnt = 0;
    for (int v = 0; v < n; v++) {
      int g = gam[(size_t)v * k + col[v]]; if (!g) continue;
      for (int c = 0; c < k; c++) { if (c == col[v]) continue;
        int d = gam[(size_t)v * k + c] - g;
        if (tabu[(size_t)v * k + c] > it && !(conf + d < best)) continue;
        if (d < bd) { bd = d; bv = v; bc = c; cnt = 1; } else if (d == bd && rnd() % (++cnt) == 0) { bv = v; bc = c; } }
    }
    if (bv < 0) continue;
    int old = col[bv]; col[bv] = bc; conf += bd;
    for (int e = ptr[bv]; e < ptr[bv + 1]; e++) { int w = nbr[e]; gam[(size_t)w * k + old]--; gam[(size_t)w * k + bc]++; }
    tabu[(size_t)bv * k + old] = it + 10 + rnd() % 10 + (long long)(0.6 * conf);
    if (conf < best) { best = conf; memcpy(bestcol, col, 4 * (size_t)n); }
  }
  int c32 = (int)best; fwrite(&c32, 4, 1, stdout); fwrite(bestcol, 4, n, stdout); return 0;
}
