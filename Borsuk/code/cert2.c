// cert2.c: certify max_{x in K_i, y in K_j} |x - y| < target for convex K = balls ∩ halfspaces in R^4 (Borsuk P8).
// Support function bound by Lagrange duality (coordinate descent on the dual, long double), directions covered by an
// adaptive subdivision of the boundary of [-1,1]^4; |h(u) - h(u0)| <= Rmax |u - u0|.  Outward safety margin 1e-12.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef long double LD;
typedef struct { int nb, nh; LD c[16][4], r2[16], a[32][4], b[32]; } Body;
static Body *Bs;
static LD dual_ub(const Body *K, const LD *u, LD *lam, LD *mu) {
  // g = |u + 2 sum lam c - sum mu a|^2 / (4 L) - sum lam (|c|^2 - r2) + sum mu b
  for (int it = 0; it < 60; it++) {
    for (int i = 0; i < K->nb; i++) {           // 1-D minimisation over lam_i by golden section on log scale
      LD lo = -30, hi = 8;
      for (int t = 0; t < 40; t++) {
        LD m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618, g1, g2, keep = lam[i];
        lam[i] = expl(m1); { LD w[4] = {u[0],u[1],u[2],u[3]}, L = 0, s = 0;
          for (int k = 0; k < K->nb; k++) { L += lam[k]; for (int d = 0; d < 4; d++) w[d] += 2 * lam[k] * K->c[k][d];
            s -= lam[k] * (K->c[k][0]*K->c[k][0]+K->c[k][1]*K->c[k][1]+K->c[k][2]*K->c[k][2]+K->c[k][3]*K->c[k][3] - K->r2[k]); }
          for (int j = 0; j < K->nh; j++) { for (int d = 0; d < 4; d++) w[d] -= mu[j] * K->a[j][d]; s += mu[j] * K->b[j]; }
          g1 = (w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]) / (4 * L) + s; }
        lam[i] = expl(m2); { LD w[4] = {u[0],u[1],u[2],u[3]}, L = 0, s = 0;
          for (int k = 0; k < K->nb; k++) { L += lam[k]; for (int d = 0; d < 4; d++) w[d] += 2 * lam[k] * K->c[k][d];
            s -= lam[k] * (K->c[k][0]*K->c[k][0]+K->c[k][1]*K->c[k][1]+K->c[k][2]*K->c[k][2]+K->c[k][3]*K->c[k][3] - K->r2[k]); }
          for (int j = 0; j < K->nh; j++) { for (int d = 0; d < 4; d++) w[d] -= mu[j] * K->a[j][d]; s += mu[j] * K->b[j]; }
          g2 = (w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]) / (4 * L) + s; }
        lam[i] = keep; if (g1 < g2) hi = m2; else lo = m1;
      }
      lam[i] = expl((lo + hi) / 2);
    }
    for (int j = 0; j < K->nh; j++) {           // closed form: quadratic in mu_j, clipped at 0
      LD w[4] = {u[0],u[1],u[2],u[3]}, L = 0;
      for (int k = 0; k < K->nb; k++) { L += lam[k]; for (int d = 0; d < 4; d++) w[d] += 2 * lam[k] * K->c[k][d]; }
      for (int jj = 0; jj < K->nh; jj++) if (jj != j) for (int d = 0; d < 4; d++) w[d] -= mu[jj] * K->a[jj][d];
      LD aa = 0, wa = 0; for (int d = 0; d < 4; d++) { aa += K->a[j][d]*K->a[j][d]; wa += w[d]*K->a[j][d]; }
      LD m = (wa - 2 * L * K->b[j]) / aa; mu[j] = m > 0 ? m : 0;
    }
  }
  LD w[4] = {u[0],u[1],u[2],u[3]}, L = 0, s = 0;
  for (int k = 0; k < K->nb; k++) { L += lam[k]; for (int d = 0; d < 4; d++) w[d] += 2 * lam[k] * K->c[k][d];
    s -= lam[k] * (K->c[k][0]*K->c[k][0]+K->c[k][1]*K->c[k][1]+K->c[k][2]*K->c[k][2]+K->c[k][3]*K->c[k][3] - K->r2[k]); }
  for (int j = 0; j < K->nh; j++) { for (int d = 0; d < 4; d++) w[d] -= mu[j] * K->a[j][d]; s += mu[j] * K->b[j]; }
  return (w[0]*w[0]+w[1]*w[1]+w[2]*w[2]+w[3]*w[3]) / (4 * L) + s + 1e-12L;
}
static void dual_parts(const Body *K, const LD *lam, const LD *mu, LD *L, LD *w, LD *s) {
  *L = 0; *s = 0; for (int d = 0; d < 4; d++) w[d] = 0;
  for (int k = 0; k < K->nb; k++) { *L += lam[k]; for (int d = 0; d < 4; d++) w[d] += 2 * lam[k] * K->c[k][d];
    *s -= lam[k] * (K->c[k][0]*K->c[k][0]+K->c[k][1]*K->c[k][1]+K->c[k][2]*K->c[k][2]+K->c[k][3]*K->c[k][3] - K->r2[k]); }
  for (int j = 0; j < K->nh; j++) { for (int d = 0; d < 4; d++) w[d] -= mu[j] * K->a[j][d]; *s += mu[j] * K->b[j]; }
}
typedef struct { LD v[4], s; int ax; } Box;
int main(void) {
  int nK; if (scanf("%d", &nK) != 1) return 1; Bs = calloc(nK, sizeof(Body));
  for (int q = 0; q < nK; q++) { Body *K = &Bs[q]; if (scanf("%d %d", &K->nb, &K->nh) != 2) return 1;
    for (int i = 0; i < K->nb; i++) if (scanf("%Lf %Lf %Lf %Lf %Lf", &K->c[i][0], &K->c[i][1], &K->c[i][2], &K->c[i][3], &K->r2[i]) != 5) return 1;
    for (int j = 0; j < K->nh; j++) if (scanf("%Lf %Lf %Lf %Lf %Lf", &K->a[j][0], &K->a[j][1], &K->a[j][2], &K->a[j][3], &K->b[j]) != 5) return 1; }
  int nQ; if (scanf("%d", &nQ) != 1) return 1;
  for (int q = 0; q < nQ; q++) {
    int i, j; LD target, Rmax; long maxb; if (scanf("%d %d %Lf %Lf %ld", &i, &j, &target, &Rmax, &maxb) != 5) return 1;
    long cap = 1 << 22; Box *st = malloc(sizeof(Box) * cap); long top = 0, n = 0; LD worst = 0; int ok = 1;
    for (int ax = 0; ax < 4; ax++) for (int sg = -1; sg <= 1; sg += 2) for (int a = 0; a < 2; a++) for (int b = 0; b < 2; b++) for (int c = 0; c < 2; c++) {
      Box B; int id[3] = {a, b, c}, t = 0; for (int d = 0; d < 4; d++) B.v[d] = (d == ax) ? sg : (-0.5L + id[t++]); B.s = 0.5L; B.ax = ax; st[top++] = B; }
    LD lam1[16], mu1[32], lam2[16], mu2[32];
    while (top > 0) {
      Box B = st[--top]; n++; if (n > maxb) { ok = 0; break; }
      LD nv = sqrtl(B.v[0]*B.v[0]+B.v[1]*B.v[1]+B.v[2]*B.v[2]+B.v[3]*B.v[3]) * (1 - 1e-15L), mv[4];
      for (int d = 0; d < 4; d++) mv[d] = -B.v[d];
      for (int k = 0; k < 16; k++) lam1[k] = lam2[k] = 0.1L; for (int k = 0; k < 32; k++) mu1[k] = mu2[k] = 0;
      LD u0[4], mu0[4]; for (int d = 0; d < 4; d++) { u0[d] = B.v[d] / nv; mu0[d] = -u0[d]; }
      dual_ub(&Bs[i], u0, lam1, mu1); dual_ub(&Bs[j], mu0, lam2, mu2);
      // for all unit u: h_i(u) + h_j(-u) <= C + <u, z>, with C, z from the fixed multipliers (weak duality)
      LD L1, w1[4], s1, L2, w2[4], s2; dual_parts(&Bs[i], lam1, mu1, &L1, w1, &s1); dual_parts(&Bs[j], lam2, mu2, &L2, w2, &s2);
      LD C = (1 + w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2]+w1[3]*w1[3]) / (4 * L1) + s1 + (1 + w2[0]*w2[0]+w2[1]*w2[1]+w2[2]*w2[2]+w2[3]*w2[3]) / (4 * L2) + s2, z[4], zn = 0, zu = 0;
      for (int d = 0; d < 4; d++) { z[d] = w1[d] / (2 * L1) - w2[d] / (2 * L2); zn += z[d]*z[d]; zu += z[d]*u0[d]; }
      zn = sqrtl(zn); LD phi0 = zn > 0 ? acosl(fminl(1, fmaxl(-1, zu / zn))) : 0;
      LD chord = B.s * sqrtl(3.0L) * 2 / nv;                     // directions of the sub-cube lie within this chord of u0
      LD theta = 2 * asinl(fminl(1, chord / 2)) * (1 + 1e-12L);
      LD quad = C + zn * cosl(fmaxl(0, phi0 - theta)), lip = C + zu + 2 * Rmax * chord;
      LD bound = fminl(quad, lip) + 1e-12L;
      if (bound < target) { if (bound > worst) worst = bound; continue; }
      if (B.s < 1.0L / 65536) { ok = 0; worst = bound; break; }
      if (top + 8 >= cap) { ok = 0; break; }
      LD h2 = B.s / 2; for (int a = 0; a < 8; a++) { Box C = B; C.s = h2; int t = 0;
        for (int d = 0; d < 4; d++) if (d != B.ax) { C.v[d] = B.v[d] + ((a >> t) & 1 ? h2 : -h2); t++; } st[top++] = C; }
    }
    printf("%d %.12Lf %ld\n", ok, worst, n); fflush(stdout); free(st);
  }
  return 0;
}
