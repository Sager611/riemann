// carve2.c -- rigorous "s-carve" of sheared prisms (pieces) against rigid-motion
// checks, for the 3D snaky-corridor sofa (corridor K = P1 u P2 u P3, boxes).
// Single-threaded.  A piece is
//   Q = { (x, yc+u, zc+v) : |u|<=hy, |v|<=hz, x = (1-s)A(u,v) + s B(u,v), s in [0,1] }
//   A(u,v) = A0 + aA u + bA v,  B(u,v) = B0 + aB u + bB v.
// For each s the slice is a planar parallelogram; a check passes at s iff for
// some box i the slice image under the check's map(s) lies in P_i shrunk by f
// (per coordinate, exact AABB of the slice image).  Checks:
//   type 0: one map (point check of a pose); f = eta.
//   type 1: chord check: two maps (endpoint poses of a motion segment whose
//           trajectories stay within eps of the chord); both images in the SAME
//           box shrunk by eps+eta.
// world = d*x + c1*y + c2*z + t.  A binary tree over consecutive checks with a
// representative map and a deviation bound culls whole ranges.
// Output: surviving s-intervals per piece.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXI 256

typedef struct { double lo, hi; } Iv;
typedef struct { double d[3], c1[3], c2[3], t[3]; } Map;
typedef struct { Map m1, m2; double eps; int type; } Check;
typedef struct { Map rep; double F[3]; int lo, hi, left, right; } Node;
typedef struct { double yc, zc, hy, hz, A0, aA, bA, B0, aB, bB; int par; } Piece;
typedef struct { double au[3][4], bu[3][4]; int nl; double al[3][4], bl[3][4]; } Lines;

static const double PL[5][3] = {{-INFINITY,0,0},{0,0,0},{0,5,0},{0,0,0},{0,5,0}};
static const double PH[5][3] = {{1,1,1},{1,6,1},{1,6,INFINITY},{1,1,1},{1,6,1}};
/* transitions a->b through the intersection box ab, chord point at t* in TS */
static const int TRA[4][3] = {{0,1,3},{1,0,3},{1,2,4},{2,1,4}};
static int NT = 0; static double TS[16];

static int NC, NP, NN, OVFL = 0;
static double ETA, QX, QY, QZ;
static Check *CK; static Node *ND; static Piece *PC;
static long long NODEVIS = 0, LEAFVIS = 0;

static void mk_lines(const Map *M, const Piece *p, Lines *L) {
  for (int m = 0; m < 3; m++) {
    double K = M->d[m]*p->A0 + M->c1[m]*p->yc + M->c2[m]*p->zc + M->t[m];
    double G = M->d[m]*(p->B0 - p->A0);
    double P = M->c1[m] + M->d[m]*p->aA, R = M->c2[m] + M->d[m]*p->bA;
    if (p->par) {
      double e = fabs(P)*p->hy + fabs(R)*p->hz;
      L->au[m][0] = K + e; L->bu[m][0] = G; L->al[m][0] = K - e; L->bl[m][0] = G;
    } else {
      double Q = M->d[m]*(p->aB - p->aA), S = M->d[m]*(p->bB - p->bA);
      int k = 0;
      for (int s1 = -1; s1 <= 1; s1 += 2) for (int s2 = -1; s2 <= 1; s2 += 2) {
        double a = K + s1*P*p->hy + s2*R*p->hz, b = G + s1*Q*p->hy + s2*S*p->hz;
        L->au[m][k] = a; L->bu[m][k] = b; L->al[m][k] = a; L->bl[m][k] = b; k++;
      }
    }
  }
  L->nl = p->par ? 1 : 4;
}

/* intersect [*lo,*hi] with {s : lines of L inside box i shrunk by f} */
static void box_s(const Lines *L, int i, const double f[3], double *lo, double *hi) {
  double a0 = *lo, b0 = *hi;
  for (int m = 0; m < 3 && a0 <= b0; m++) {
    double H = PH[i][m] - f[m], Lo = PL[i][m] + f[m];
    int fh = isfinite(H), fl = isfinite(Lo);
    for (int k = 0; k < L->nl; k++) {
      if (fh) { double a = L->au[m][k], b = L->bu[m][k];
        if (b > 0) { double s = (H - a) / b; if (s < b0) b0 = s; }
        else if (b < 0) { double s = (H - a) / b; if (s > a0) a0 = s; }
        else if (!(a <= H)) { a0 = INFINITY; b0 = -INFINITY; } }
      if (fl) { double a = L->al[m][k], b = L->bl[m][k];
        if (b > 0) { double s = (Lo - a) / b; if (s > a0) a0 = s; }
        else if (b < 0) { double s = (Lo - a) / b; if (s < b0) b0 = s; }
        else if (!(a >= Lo)) { a0 = INFINITY; b0 = -INFINITY; } }
    }
  }
  *lo = a0; *hi = b0;
}

static int union3(Iv J[3], Iv *U) {
  Iv t[3]; int n = 0;
  for (int i = 0; i < 3; i++) if (J[i].lo <= J[i].hi) t[n++] = J[i];
  for (int i = 1; i < n; i++) for (int j = i; j > 0 && t[j].lo < t[j-1].lo; j--) { Iv s = t[j]; t[j] = t[j-1]; t[j-1] = s; }
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (m > 0 && t[i].lo <= U[m-1].hi) { if (t[i].hi > U[m-1].hi) U[m-1].hi = t[i].hi; }
    else U[m++] = t[i];
  }
  return m;
}

static int unionN(Iv *J, int n0, Iv *U) {
  Iv t[64]; int n = 0;
  for (int i = 0; i < n0; i++) if (J[i].lo <= J[i].hi) t[n++] = J[i];
  for (int i = 1; i < n; i++) for (int j = i; j > 0 && t[j].lo < t[j-1].lo; j--) { Iv s = t[j]; t[j] = t[j-1]; t[j-1] = s; }
  int m = 0;
  for (int i = 0; i < n; i++) {
    if (m > 0 && t[i].lo <= U[m-1].hi) { if (t[i].hi > U[m-1].hi) U[m-1].hi = t[i].hi; }
    else U[m++] = t[i];
  }
  return m;
}

static int intersect(const Iv *S, int ns, const Iv *U, int nu, Iv *O) {
  int i = 0, j = 0, m = 0;
  while (i < ns && j < nu) {
    double a = S[i].lo > U[j].lo ? S[i].lo : U[j].lo;
    double b = S[i].hi < U[j].hi ? S[i].hi : U[j].hi;
    if (a <= b) { if (m < MAXI) { O[m].lo = a; O[m].hi = b; } m++; }
    if (S[i].hi < U[j].hi) i++; else j++;
  }
  return m;
}

static int covered(const Iv *S, int ns, const Iv *U, int nu) {
  int j = 0;
  for (int i = 0; i < ns; i++) {
    while (j < nu && U[j].hi < S[i].hi) j++;
    if (j >= nu || U[j].lo > S[i].lo || U[j].hi < S[i].hi) return 0;
  }
  return 1;
}

static void apply_check(int k, const Piece *p, Iv *S, int *ns) {
  const Check *c = &CK[k];
  Lines L1, L2; mk_lines(&c->m1, p, &L1); if (c->type) mk_lines(&c->m2, p, &L2);
  double f[3]; for (int m = 0; m < 3; m++) f[m] = c->eps + ETA;
  Iv J[64], U[64], O[MAXI];
  int nj = 0;
  for (int i = 0; i < 3; i++) {
    double lo = S[0].lo, hi = S[*ns-1].hi;
    box_s(&L1, i, f, &lo, &hi);
    if (c->type && lo <= hi) box_s(&L2, i, f, &lo, &hi);
    J[nj].lo = lo; J[nj].hi = hi; nj++;
  }
  int nu = unionN(J, nj, U);
  if (c->type && NT > 0 && !covered(S, *ns, U, nu)) {
    for (int q = 0; q < NT; q++) {
      double ts = TS[q]; Map Mm; Lines Lm;
      for (int m = 0; m < 3; m++) {
        Mm.d[m] = (1 - ts) * c->m1.d[m] + ts * c->m2.d[m]; Mm.c1[m] = (1 - ts) * c->m1.c1[m] + ts * c->m2.c1[m];
        Mm.c2[m] = (1 - ts) * c->m1.c2[m] + ts * c->m2.c2[m]; Mm.t[m] = (1 - ts) * c->m1.t[m] + ts * c->m2.t[m];
      }
      mk_lines(&Mm, p, &Lm);
      for (int r = 0; r < 4; r++) {
        double lo = S[0].lo, hi = S[*ns-1].hi;
        box_s(&L1, TRA[r][0], f, &lo, &hi);
        if (lo <= hi) box_s(&L2, TRA[r][1], f, &lo, &hi);
        if (lo <= hi) box_s(&Lm, TRA[r][2], f, &lo, &hi);
        if (lo <= hi && nj < 64) { J[nj].lo = lo; J[nj].hi = hi; nj++; }
      }
    }
    nu = unionN(J, nj, U);
  }
  int m = intersect(S, *ns, U, nu, O);
  if (m > MAXI) { OVFL = 1; m = MAXI; }
  memcpy(S, O, sizeof(Iv)*m); *ns = m;
}

static void visit(int n, const Piece *p, Iv *S, int *ns) {
  if (*ns == 0) return;
  const Node *nd = &ND[n];
  Lines L; mk_lines(&nd->rep, p, &L);
  Iv J[3], U[3];
  for (int i = 0; i < 3; i++) { double lo = S[0].lo, hi = S[*ns-1].hi; box_s(&L, i, nd->F, &lo, &hi); J[i].lo = lo; J[i].hi = hi; }
  int nu = union3(J, U);
  if (covered(S, *ns, U, nu)) return;
  if (nd->left < 0) { LEAFVIS++; apply_check(nd->lo, p, S, ns); return; }
  NODEVIS++;
  visit(nd->left, p, S, ns);
  visit(nd->right, p, S, ns);
}

static void dev_acc(const Map *a, const Map *r, double dev[3]) {
  for (int m = 0; m < 3; m++) {
    double v = fabs(a->d[m]-r->d[m])*QX + fabs(a->c1[m]-r->c1[m])*QY + fabs(a->c2[m]-r->c2[m])*QZ + fabs(a->t[m]-r->t[m]);
    if (v > dev[m]) dev[m] = v;
  }
}

static int build(int lo, int hi) {
  int id = NN++;
  int mid = (lo + hi) / 2;
  ND[id].lo = lo; ND[id].hi = hi; ND[id].rep = CK[mid].m1;
  double dev[3] = {0,0,0}, me = 0;
  for (int k = lo; k < hi; k++) {
    dev_acc(&CK[k].m1, &ND[id].rep, dev);
    if (CK[k].type) dev_acc(&CK[k].m2, &ND[id].rep, dev);
    if (CK[k].eps > me) me = CK[k].eps;
  }
  for (int m = 0; m < 3; m++) ND[id].F[m] = (dev[m] + me) * (1 + 1e-12) + ETA + 1e-15;
  if (hi - lo == 1) { ND[id].left = ND[id].right = -1; }
  else { int L = build(lo, mid); int R = build(mid, hi); ND[id].left = L; ND[id].right = R; }
  return id;
}

static double rd(FILE *f) { double v; if (fread(&v, 8, 1, f) != 1) { fprintf(stderr, "read error\n"); exit(2); } return v; }
static void rdmap(FILE *f, Map *M) { for (int m = 0; m < 3; m++) M->d[m] = rd(f); for (int m = 0; m < 3; m++) M->c1[m] = rd(f); for (int m = 0; m < 3; m++) M->c2[m] = rd(f); for (int m = 0; m < 3; m++) M->t[m] = rd(f); }

int main(int argc, char **argv) {
  if (argc < 3) { fprintf(stderr, "usage: carve2 in.bin out.bin\n"); return 1; }
  FILE *f = fopen(argv[1], "rb"); if (!f) { perror("open"); return 1; }
  NC = (int)rd(f); NP = (int)rd(f); ETA = rd(f); QX = rd(f); QY = rd(f); QZ = rd(f);
  NT = (int)rd(f); if (NT < 0 || NT > 16) { fprintf(stderr, "bad NT\n"); return 1; }
  for (int q = 0; q < NT; q++) TS[q] = rd(f);
  CK = calloc(NC, sizeof(Check)); PC = calloc(NP, sizeof(Piece));
  for (int k = 0; k < NC; k++) { CK[k].type = (int)rd(f); rdmap(f, &CK[k].m1); rdmap(f, &CK[k].m2); CK[k].eps = rd(f); }
  for (int j = 0; j < NP; j++) {
    Piece *p = &PC[j];
    p->yc = rd(f); p->zc = rd(f); p->hy = rd(f); p->hz = rd(f);
    p->A0 = rd(f); p->aA = rd(f); p->bA = rd(f); p->B0 = rd(f); p->aB = rd(f); p->bB = rd(f);
    p->par = (p->aA == p->aB && p->bA == p->bB);
  }
  fclose(f);
  ND = calloc(2 * NC + 2, sizeof(Node)); NN = 0; build(0, NC);
  int *cnt = calloc(NP, sizeof(int));
  size_t cap = (size_t)NP + 1024, tot = 0; Iv *all = malloc(sizeof(Iv) * cap);
  Iv S[MAXI];
  for (int j = 0; j < NP; j++) {
    int ns = 1; S[0].lo = 0.0; S[0].hi = 1.0;
    visit(0, &PC[j], S, &ns);
    cnt[j] = ns;
    if (tot + ns > cap) { cap = 2 * cap + ns; all = realloc(all, sizeof(Iv) * cap); }
    memcpy(&all[tot], S, sizeof(Iv) * ns); tot += ns;
  }
  FILE *o = fopen(argv[2], "wb");
  int hdr[4] = {NP, (int)tot, 0, OVFL};
  fwrite(hdr, 4, 4, o); fwrite(cnt, sizeof(int), NP, o); fwrite(all, sizeof(Iv), tot, o);
  fclose(o);
  fprintf(stderr, "checks %d pieces %d nodes %d nodevisits %lld leafvisits %lld intervals %zu overflow %d\n", NC, NP, NN, NODEVIS, LEAFVIS, tot, OVFL);
  return OVFL ? 3 : 0;
}
