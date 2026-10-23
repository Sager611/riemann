/* verify_core.c -- floating-point core of verify.py (rigorous: see NOTE.md, section "Rounding").
 * Compile: cc -O3 -ffp-contract=off -o vcore verify_core.c;  run: vcore checks.bin pieces.bin out.bin
 * Input (float64): NC NP ETA QX QY QZ BS FLIM C12, then NC checks x 26: type, map1[12], map2[12], eps
 *   (map = R[:,0], R[:,1], R[:,2], T: world_m = R_m0 x + R_m1 y + R_m2 z + T_m),
 * Second input file: 16 blocks of NP float64: for corner c = 0..3: a_c, d_c, y_c, z_c of every piece
 * (corner c of the slice at s is (a_c + s d_c, y_c, z_c)).
 * Output: NP bytes, 1 iff every s in [0,1] passes every check.
 * Type 0 check passes at s if the slice image under map1 lies in one box of K shrunk by eps+ETA;
 * type 1: images under map1 AND map2 lie in the SAME box shrunk by eps+ETA.
 * Blocks of BS and of 8 consecutive checks are skipped when one representative map with margin
 * F = (dev + max eps) * C12 + ETA passes (dev bounds how far any map of the block moves a body point). */
#pragma STDC FP_CONTRACT OFF
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct { double lo, hi; } I;
typedef struct { double a[4], d[4], y[4], z[4]; } Pc;
static const double LO[3][3] = {{-INFINITY, 0, 0}, {0, 0, 0}, {0, 5, 0}};
static const double HI[3][3] = {{1, 1, 1}, {1, 6, 1}, {1, 6, INFINITY}};
static int nc, np_, bs; static double eta, qx, qy, qz, flim, c12;
static double *ck, *bk, *sb;

static double rd(FILE *f) { double v; if (fread(&v, 8, 1, f) != 1) { fprintf(stderr, "read error\n"); exit(2); } return v; }

/* al = ((R_m0 a + R_m1 y) + R_m2 z) + T_m and be = R_m0 d, for m = 0..2 and corners c = 0..3 */
static void lines(const double *M, const Pc *p, double *al, double *be) {
  for (int m = 0; m < 3; m++) for (int c = 0; c < 4; c++) {
    double t = M[m] * p->a[c];
    double t2 = M[3 + m] * p->y[c];
    t = t + t2;
    t2 = M[6 + m] * p->z[c];
    t = t + t2;
    al[4 * m + c] = t + M[9 + m];
    be[4 * m + c] = M[m] * p->d[c];
  }
}

/* shrink [lo,hi] to the s with  L_m + f_m <= al + be s <= H_m - f_m  for all m, c  (box b) */
static void cut(const double *al, const double *be, int b, const double *f, double *lo, double *hi) {
  for (int m = 0; m < 3 && *lo <= *hi; m++) {
    double L = LO[b][m] + f[m], H = HI[b][m] - f[m];
    for (int c = 0; c < 4 && *lo <= *hi; c++) {
      double A = al[4 * m + c], B = be[4 * m + c];
      if (H < INFINITY) {
        if (B > 0) { double s = (H - A) / B; if (s < *hi) *hi = s; }
        else if (B < 0) { double s = (H - A) / B; if (s > *lo) *lo = s; }
        else if (A > H) { *lo = 1; *hi = 0; }
      }
      if (L > -INFINITY) {
        if (B > 0) { double s = (L - A) / B; if (s > *lo) *lo = s; }
        else if (B < 0) { double s = (L - A) / B; if (s < *hi) *hi = s; }
        else if (A < L) { *lo = 1; *hi = 0; }
      }
    }
  }
}

static int unite(I *J, int n, I *U) {
  for (int i = 1; i < n; i++) for (int j = i; j > 0 && J[j].lo < J[j - 1].lo; j--) { I t = J[j]; J[j] = J[j - 1]; J[j - 1] = t; }
  int m = 0;
  for (int i = 0; i < n; i++) { if (J[i].lo > J[i].hi) continue; if (m && J[i].lo <= U[m - 1].hi) { if (J[i].hi > U[m - 1].hi) U[m - 1].hi = J[i].hi; } else U[m++] = J[i]; }
  return m;
}
static int inside(const I *S, int ns, const I *U, int nu) {
  for (int i = 0; i < ns; i++) { int ok = 0; for (int j = 0; j < nu; j++) if (U[j].lo <= S[i].lo && S[i].hi <= U[j].hi) { ok = 1; break; } if (!ok) return 0; }
  return 1;
}
static int meet(const I *S, int ns, const I *U, int nu, I *O) {
  int m = 0;
  for (int i = 0; i < ns; i++) for (int j = 0; j < nu; j++) {
    double a = S[i].lo > U[j].lo ? S[i].lo : U[j].lo, b = S[i].hi < U[j].hi ? S[i].hi : U[j].hi;
    if (a <= b && m < 64) { O[m].lo = a; O[m].hi = b; m++; }
  }
  for (int i = 1; i < m; i++) for (int j = i; j > 0 && O[j].lo < O[j - 1].lo; j--) { I t = O[j]; O[j] = O[j - 1]; O[j - 1] = t; }
  return m;
}

static void build(double *out, int size, int nblk) {
  for (int b = 0; b < nblk; b++) {
    double *R = out + 15 * b; const double *r0 = ck + 26 * (b * size) + 1;
    for (int i = 0; i < 12; i++) R[i] = r0[i];
    double dev[3] = {0, 0, 0}, me = 0;
    for (int k = b * size; k < nc && k < (b + 1) * size; k++) {
      const double *c = ck + 26 * k;
      for (int mp = 0; mp < (c[0] > 0.5 ? 2 : 1); mp++) {
        const double *M = c + 1 + 12 * mp;
        for (int m = 0; m < 3; m++) {
          double v = fabs(M[m] - R[m]) * qx;
          double t = fabs(M[3 + m] - R[3 + m]) * qy;
          v = v + t;
          t = fabs(M[6 + m] - R[6 + m]) * qz;
          v = v + t;
          v = v + fabs(M[9 + m] - R[9 + m]);
          if (v > dev[m]) dev[m] = v;
        }
      }
      if (c[25] > me) me = c[25];
    }
    for (int m = 0; m < 3; m++) {
      double F = dev[m] + me;
      F = F * c12;
      F = F + eta;
      if (!(F <= flim)) { fprintf(stderr, "block margin %g above FLIM\n", F); exit(4); }
      R[12 + m] = F;
    }
  }
}

int main(int argc, char **argv) {
  if (argc < 4) return 1;
  FILE *f = fopen(argv[1], "rb"); if (!f) return 1;
  nc = (int)rd(f); np_ = (int)rd(f); eta = rd(f); qx = rd(f); qy = rd(f); qz = rd(f); bs = (int)rd(f); flim = rd(f); c12 = rd(f);
  if (bs % 8 != 0 || bs <= 0) { fprintf(stderr, "BS must be a positive multiple of 8\n"); return 1; }
  ck = malloc(sizeof(double) * 26 * nc);
  for (long i = 0; i < 26L * nc; i++) { ck[i] = rd(f); if (!isfinite(ck[i])) { fprintf(stderr, "non-finite check data\n"); return 5; } }
  int nb = (nc + bs - 1) / bs, n8 = (nc + 7) / 8;
  bk = malloc(sizeof(double) * 15 * nb); sb = malloc(sizeof(double) * 15 * n8);
  build(bk, bs, nb); build(sb, 8, n8);
  /* chain[k] = 1 iff map1 of check k equals map2 of check k-1 (then its lines are reused) */
  char *chain = calloc(nc, 1);
  for (int k = 1; k < nc; k++) {
    const double *c0 = ck + 26 * (k - 1), *c1 = ck + 26 * k; int eq = c0[0] > 0.5 && c1[0] > 0.5;
    for (int i = 0; i < 12 && eq; i++) eq = (c1[1 + i] == c0[13 + i]);
    chain[k] = (char)eq;
  }
  FILE *g = fopen(argv[2], "rb"); if (!g) return 1;
  double *pd = malloc(sizeof(double) * 16 * (size_t)np_);
  if (fread(pd, sizeof(double), 16 * (size_t)np_, g) != 16 * (size_t)np_) { fprintf(stderr, "piece read error\n"); return 2; }
  fclose(g);
  for (size_t i = 0; i < 16 * (size_t)np_; i++) if (!isfinite(pd[i])) { fprintf(stderr, "non-finite piece data\n"); return 5; }
  FILE *o = fopen(argv[3], "wb");
  long long exact = 0; int nfull = 0;
  double buf[4][12]; double *al = buf[0], *be = buf[1], *al2 = buf[2], *be2 = buf[3];
  for (int p = 0; p < np_; p++) {
    Pc pc;
    for (int c = 0; c < 4; c++) {
      pc.a[c] = pd[(size_t)(4 * c + 0) * np_ + p]; pc.d[c] = pd[(size_t)(4 * c + 1) * np_ + p];
      pc.y[c] = pd[(size_t)(4 * c + 2) * np_ + p]; pc.z[c] = pd[(size_t)(4 * c + 3) * np_ + p];
    }
    I S[64], J[3], U[3], O[64]; int ns = 1, nu; S[0].lo = 0; S[0].hi = 1;
    for (int b = 0; b < nb && ns > 0; b++) {
      const double *R = bk + 15 * b;
      lines(R, &pc, al, be);
      for (int i = 0; i < 3; i++) { J[i].lo = S[0].lo; J[i].hi = S[ns - 1].hi; cut(al, be, i, R + 12, &J[i].lo, &J[i].hi); }
      nu = unite(J, 3, U);
      if (inside(S, ns, U, nu)) continue;
      for (int q8 = (b * bs) / 8; q8 * 8 < nc && q8 * 8 < (b + 1) * bs && ns > 0; q8++) {
        const double *R8 = sb + 15 * q8;
        lines(R8, &pc, al, be);
        for (int i = 0; i < 3; i++) { J[i].lo = S[0].lo; J[i].hi = S[ns - 1].hi; cut(al, be, i, R8 + 12, &J[i].lo, &J[i].hi); }
        nu = unite(J, 3, U);
        if (inside(S, ns, U, nu)) continue;
        int lastk = -2;
        for (int k = q8 * 8; k < nc && k < (q8 + 1) * 8 && ns > 0; k++) {
          const double *ch = ck + 26 * k;
          int two = ch[0] > 0.5;
          double fe0 = ch[25] + eta; double fe[3] = {fe0, fe0, fe0};
          exact++;
          if (k == lastk + 1 && chain[k]) { double *t = al; al = al2; al2 = t; t = be; be = be2; be2 = t; }
          else lines(ch + 1, &pc, al, be);
          if (two) lines(ch + 13, &pc, al2, be2);
          lastk = two ? k : -2;
          for (int i = 0; i < 3; i++) {
            J[i].lo = S[0].lo; J[i].hi = S[ns - 1].hi;
            cut(al, be, i, fe, &J[i].lo, &J[i].hi);
            if (two && J[i].lo <= J[i].hi) cut(al2, be2, i, fe, &J[i].lo, &J[i].hi);
          }
          nu = unite(J, 3, U);
          ns = meet(S, ns, U, nu, O);
          for (int i = 0; i < ns; i++) S[i] = O[i];
        }
      }
    }
    int full = (ns == 1 && S[0].lo <= 0.0 && S[0].hi >= 1.0);
    nfull += full; fputc(full, o);
  }
  fclose(o); fclose(f);
  fprintf(stderr, "pieces %d full %d checks %d blocks %d exact %lld\n", np_, nfull, nc, nb, exact);
  return 0;
}
