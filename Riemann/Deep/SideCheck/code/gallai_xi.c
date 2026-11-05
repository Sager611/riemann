// sidecheck: exact max order F(b1,b2,b3) of a Gallai 3-colouring (no rainbow triangle) of K_n with
// alpha_c <= b_c for each colour c (alpha_c = largest vertex set with no colour-c edge).
// F(5,5,5) = (largest Gallai colouring in which every K_6 sees all three colours).
// Recursion via the Gallai partition: X = V_1..V_t (t>=2), reduced K_t coloured with two colours p,q
// (third colour r absent between parts). alpha_r(X) = sum alpha_r(V_i); alpha_p(X) = max over sets of parts
// pairwise joined in colour q of sum alpha_p(V_i); alpha_q symmetric.  Parts with beta_r = 1 are r-cliques,
// F(.,.,1) = min of the other two budgets.
#include <stdio.h>
#include <string.h>
static int F[7][7][7]; static int done[7][7][7];
static int T, adj[6][6];           // adj[i][j] = 1 if parts i,j joined in colour p, 0 if colour q
static int bud[3], rr, pp, qq, best;
static int br[6], bp[6], bq[6];
static int Fv(int a, int b, int c);
static int ok_constraints(void) {
  // every set I of parts pairwise joined in colour q (adj=0): sum bp <= bud[pp]
  // every set J pairwise joined in colour p (adj=1): sum bq <= bud[qq]
  for (int m = 1; m < (1 << T); m++) {
    int indep = 1, cliq = 1, sp = 0, sq = 0;
    for (int i = 0; i < T; i++) if (m >> i & 1) {
      sp += bp[i]; sq += bq[i];
      for (int j = i + 1; j < T; j++) if (m >> j & 1) { if (adj[i][j]) indep = 0; else cliq = 0; }
    }
    if (indep && sp > bud[pp]) return 0;
    if (cliq && sq > bud[qq]) return 0;
  }
  return 1;
}
static int value(void) {
  int s = 0;
  for (int i = 0; i < T; i++) {
    int v[3]; v[rr] = br[i]; v[pp] = bp[i]; v[qq] = bq[i];
    s += Fv(v[0], v[1], v[2]);
  }
  return s;
}
static void rec(int i, int sumr) {
  if (i == T) { if (ok_constraints()) { int v = value(); if (v > best) best = v; } return; }
  for (int x = 1; sumr + x + (T - 1 - i) <= bud[rr]; x++) {
    br[i] = x;
    if (x == 1) { // r-clique part: only min(bp,bq) matters -> take bp = bq = m
      for (int m = 1; m <= 5; m++) { bp[i] = m; bq[i] = m; if (m > bud[pp] || m > bud[qq]) break; rec(i + 1, sumr + x); }
    } else {
      for (int a = 1; a <= bud[pp]; a++) for (int b = 1; b <= bud[qq]; b++) { bp[i] = a; bq[i] = b; rec(i + 1, sumr + x); }
    }
  }
}
static int Fv(int a, int b, int c) {
  if (a <= 0 || b <= 0 || c <= 0) return 0;
  if (done[a][b][c]) return F[a][b][c];
  int res = 1;                                   // single vertex
  int bb[3] = {a, b, c};
  for (int r = 0; r < 3; r++) {
    int p = (r + 1) % 3, q = (r + 2) % 3;
    for (int t = 2; t <= 5 && t <= bb[r]; t++) {
      int ne = t * (t - 1) / 2;
      for (int g = 0; g < (1 << ne); g++) {
        // save globals (recursion re-enters Fv)
        int sT = T, sadj[6][6], sb[3], sr = rr, sp = pp, sq = qq, sbest = best, sbr[6], sbp[6], sbq[6];
        memcpy(sadj, adj, sizeof adj); memcpy(sb, bud, sizeof bud); memcpy(sbr, br, sizeof br); memcpy(sbp, bp, sizeof bp); memcpy(sbq, bq, sizeof bq);
        T = t; int e = 0;
        for (int i = 0; i < t; i++) for (int j = i + 1; j < t; j++) { adj[i][j] = adj[j][i] = (g >> e) & 1; e++; }
        bud[0] = a; bud[1] = b; bud[2] = c; rr = r; pp = p; qq = q; best = 0;
        rec(0, 0);
        if (best > res) res = best;
        T = sT; memcpy(adj, sadj, sizeof adj); memcpy(bud, sb, sizeof bud); rr = sr; pp = sp; qq = sq; best = sbest;
        memcpy(br, sbr, sizeof br); memcpy(bp, sbp, sizeof bp); memcpy(bq, sbq, sizeof bq);
      }
    }
  }
  done[a][b][c] = 1; F[a][b][c] = res; return res;
}
int main(void) {
  // sanity: known Gallai values?  print F(k-1,k-1,k-1) for k = 3..6
  for (int k = 2; k <= 6; k++) printf("F(%d,%d,%d) = %d   (largest Gallai 3-colouring with every K_%d 3-coloured)\n", k-1, k-1, k-1, Fv(k-1, k-1, k-1), k);
  printf("F(5,5,1)=%d F(5,4,3)=%d F(4,4,4)=%d\n", Fv(5,5,1), Fv(5,4,3), Fv(4,4,4));
  printf("checks (Budden-DeJonge gr-1): F(2,3,3)=%d [4] F(2,4,4)=%d [5] F(2,5,5)=%d [6] F(3,3,3)=%d [4] F(4,3,3)=%d [5] F(5,3,3)=%d [6] F(4,4,3)=%d [6] F(5,4,3)=%d [7]\n", Fv(2,3,3), Fv(2,4,4), Fv(2,5,5), Fv(3,3,3), Fv(4,3,3), Fv(5,3,3), Fv(4,4,3), Fv(5,4,3));
  return 0;
}
