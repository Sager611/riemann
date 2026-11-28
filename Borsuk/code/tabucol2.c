// tabucol2.c: k-colouring with hard "far" edges (weight W, endpoints must differ) and soft "near" edges (weight 1,
// endpoints should agree).  Tabu search on total cost.  stdin (int32): n, mf, mn, k, iters, seed, W, n init colours,
// mf far edges (u,v), mn near edges (u,v).  stdout (int32): far conflicts, near disagreements, n colours (best found:
// fewest far conflicts, then fewest near disagreements).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static unsigned long long rs = 88172645463325252ULL;
static unsigned rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return (unsigned)(rs >> 11); }
static int *mkcsr(int n, int m, int *U, int *V, int **ptr_out) {
  int *deg = calloc(n + 1, 4), *ptr = malloc(4 * (size_t)(n + 1)), *nb = malloc(8 * (size_t)m + 4), *fill = calloc(n, 4);
  for (int e = 0; e < m; e++) { deg[U[e]]++; deg[V[e]]++; }
  ptr[0] = 0; for (int v = 0; v < n; v++) ptr[v + 1] = ptr[v] + deg[v];
  for (int e = 0; e < m; e++) { nb[ptr[U[e]] + fill[U[e]]++] = V[e]; nb[ptr[V[e]] + fill[V[e]]++] = U[e]; }
  free(deg); free(fill); *ptr_out = ptr; return nb;
}
int main(void) {
  int h[7]; if (fread(h, 4, 7, stdin) != 7) return 1;
  int n = h[0], mf = h[1], mn = h[2], k = h[3], iters = h[4], W = h[6]; rs ^= (unsigned long long)h[5] * 2654435761ULL;
  int *col = malloc(4 * (size_t)n); if (fread(col, 4, n, stdin) != (size_t)n) return 1;
  int *FU = malloc(4 * (size_t)mf + 4), *FV = malloc(4 * (size_t)mf + 4), *NU = malloc(4 * (size_t)mn + 4), *NV = malloc(4 * (size_t)mn + 4);
  for (int e = 0; e < mf; e++) { int uv[2]; if (fread(uv, 4, 2, stdin) != 2) return 1; FU[e] = uv[0]; FV[e] = uv[1]; }
  for (int e = 0; e < mn; e++) { int uv[2]; if (fread(uv, 4, 2, stdin) != 2) return 1; NU[e] = uv[0]; NV[e] = uv[1]; }
  int *fp, *np_; int *fn = mkcsr(n, mf, FU, FV, &fp), *nn = mkcsr(n, mn, NU, NV, &np_);
  long long *gam = calloc((size_t)n * k, 8), *tabu = calloc((size_t)n * k, 8);
  for (int v = 0; v < n; v++) {
    for (int e = fp[v]; e < fp[v + 1]; e++) gam[(size_t)v * k + col[fn[e]]] += W;
    int nd = np_[v + 1] - np_[v];
    for (int c = 0; c < k; c++) gam[(size_t)v * k + c] += nd;
    for (int e = np_[v]; e < np_[v + 1]; e++) gam[(size_t)v * k + col[nn[e]]] -= 1;
  }
  long long cost = 0; for (int v = 0; v < n; v++) cost += gam[(size_t)v * k + col[v]]; cost /= 2;
  long long best = cost; int *bestcol = malloc(4 * (size_t)n); memcpy(bestcol, col, 4 * (size_t)n);
  for (long long it = 0; it < iters; it++) {
    long long bd = (long long)1 << 60; int bv = -1, bc = -1, cnt = 0;
    for (int v = 0; v < n; v++) {
      long long g = gam[(size_t)v * k + col[v]]; if (g <= 0) continue;
      for (int c = 0; c < k; c++) { if (c == col[v]) continue;
        long long d = gam[(size_t)v * k + c] - g;
        if (tabu[(size_t)v * k + c] > it && !(cost + d < best)) continue;
        if (d < bd) { bd = d; bv = v; bc = c; cnt = 1; } else if (d == bd && rnd() % (++cnt) == 0) { bv = v; bc = c; } }
    }
    if (bv < 0) continue;
    int old = col[bv]; col[bv] = bc; cost += bd;
    for (int e = fp[bv]; e < fp[bv + 1]; e++) { int w = fn[e]; gam[(size_t)w * k + old] -= W; gam[(size_t)w * k + bc] += W; }
    for (int e = np_[bv]; e < np_[bv + 1]; e++) { int w = nn[e]; gam[(size_t)w * k + old] += 1; gam[(size_t)w * k + bc] -= 1; }
    tabu[(size_t)bv * k + old] = it + 7 + rnd() % 10;
    if (cost < best) { best = cost; memcpy(bestcol, col, 4 * (size_t)n); }
  }
  int fc = 0, nd = 0;
  for (int e = 0; e < mf; e++) fc += bestcol[FU[e]] == bestcol[FV[e]];
  for (int e = 0; e < mn; e++) nd += bestcol[NU[e]] != bestcol[NV[e]];
  fwrite(&fc, 4, 1, stdout); fwrite(&nd, 4, 1, stdout); fwrite(bestcol, 4, n, stdout); return 0;
}
