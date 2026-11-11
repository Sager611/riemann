// rank.c (explore28): numerical ranks of the off-diagonal blocks of the Cauchy-like matrix J'' (nu = 0), HODLR partition of the
// points sorted by y, after the exact Jacobi scaling used by slt.c. Each block M_IJ is Cauchy-like, M_ab = (P_a . Q_b)/(y_a - y_b)
// with P_a = (u_a, w_a), Q_b = (w_b, -u_b); adaptive cross approximation (rook pivoting) with exact generator (GKO) updates of the
// block's Schur complement; the pivot magnitudes |M_pq| after r steps estimate the (r+1)-th singular value scale.
// Usage: rank <gendump> <prec> <leaf> tol1 tol2 ...   (tolerances relative to max |J''_ab| = max |diag| after scaling)
#include <flint/arb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
static char buf[1 << 16];
static const char *tok(FILE *f) { if (fscanf(f, "%65535s", buf) != 1) { fprintf(stderr, "read error\n"); exit(2); } return buf; }
static slong n, prec; static arb_ptr y, u, w, d; static double *yd; static slong *ord;
static int cmpy(const void *a, const void *b) { double x = yd[*(const slong *) a], z = yd[*(const slong *) b]; return x < z ? -1 : x > z; }
static double getd(const arb_t x) { return arf_get_d(arb_midref(x), ARF_RND_NEAR); }
static int ntol; static double *tol; static double scale;
static long long storeH[64][32]; static int maxr[64][32]; static int nblk[64]; static slong bsz[64];
static void block(slong i0, slong m, slong j0, slong k, int lev) {   // rows ord[i0..i0+m), cols ord[j0..j0+k)
    arb_ptr P1 = _arb_vec_init(m), P2 = _arb_vec_init(m), Q1 = _arb_vec_init(k), Q2 = _arb_vec_init(k);
    arb_ptr yr = _arb_vec_init(m), yc = _arb_vec_init(k);
    char *ur = calloc(m, 1), *uc = calloc(k, 1);
    for (slong a = 0; a < m; a++) { slong s = ord[i0 + a]; arb_set(P1 + a, u + s); arb_set(P2 + a, w + s); arb_set(yr + a, y + s); }
    for (slong b = 0; b < k; b++) { slong s = ord[j0 + b]; arb_set(Q1 + b, w + s); arb_neg(Q2 + b, u + s); arb_set(yc + b, y + s); }
    arb_t t, e, piv, ip; arb_init(t); arb_init(e); arb_init(piv); arb_init(ip);
    #define ENT(res, a, b) do { arb_mul(res, P1 + (a), Q1 + (b), prec); arb_addmul(res, P2 + (a), Q2 + (b), prec); arb_sub(e, yr + (a), yc + (b), prec); arb_div(res, res, e, prec); } while (0)
    arb_ptr col = _arb_vec_init(m), row = _arb_vec_init(k);
    int r = 0, rk[32]; for (int i = 0; i < ntol; i++) rk[i] = -1;
    slong maxr_ = m < k ? m : k; double tmin = tol[ntol - 1] * scale;
    unsigned long long rs = 88172645463325252ULL;
    while (r < maxr_) {
        slong q = -1; for (slong b = 0; b < k; b++) if (!uc[b]) { q = b; break; }
        slong p = -1; double best = -1;
        for (int it = 0; it < 4; it++) {
            best = -1; for (slong a = 0; a < m; a++) if (!ur[a]) { ENT(col + a, a, q); double v = fabs(getd(col + a)); if (v > best) { best = v; p = a; } }
            double best2 = -1; slong q2 = q;
            for (slong b = 0; b < k; b++) if (!uc[b]) { ENT(row + b, p, b); double v = fabs(getd(row + b)); if (v > best2) { best2 = v; q2 = b; } }
            if (q2 == q) break; q = q2;
        }
        // fresh column and row of the pivot
        for (slong a = 0; a < m; a++) if (!ur[a]) ENT(col + a, a, q);
        for (slong b = 0; b < k; b++) if (!uc[b]) ENT(row + b, p, b);
        arb_set(piv, col + p); double pv = fabs(getd(piv));
        if (pv < tmin) {   // random residual probes before stopping
            double mx = 0; slong pa = -1, pb = -1;
            for (int s = 0; s < 64; s++) {
                rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; slong a = rs % m; rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; slong b = rs % k;
                if (ur[a] || uc[b]) continue; ENT(t, a, b); double v = fabs(getd(t)); if (v > mx) { mx = v; pa = a; pb = b; }
            }
            if (mx < tmin) break;
            p = pa; q = pb; for (slong a = 0; a < m; a++) if (!ur[a]) ENT(col + a, a, q); for (slong b = 0; b < k; b++) if (!uc[b]) ENT(row + b, p, b);
            arb_set(piv, col + p); pv = fabs(getd(piv));
        }
        for (int i = 0; i < ntol; i++) if (rk[i] < 0 && pv < tol[i] * scale) rk[i] = r;
        arb_inv(ip, piv, prec);
        arb_t g1, g2; arb_init(g1); arb_init(g2); arb_set(g1, P1 + p); arb_set(g2, P2 + p);
        for (slong a = 0; a < m; a++) if (!ur[a] && a != p) { arb_mul(t, col + a, ip, prec); arb_submul(P1 + a, t, g1, prec); arb_submul(P2 + a, t, g2, prec); }
        arb_set(g1, Q1 + q); arb_set(g2, Q2 + q);
        for (slong b = 0; b < k; b++) if (!uc[b] && b != q) { arb_mul(t, row + b, ip, prec); arb_submul(Q1 + b, t, g1, prec); arb_submul(Q2 + b, t, g2, prec); }
        arb_clear(g1); arb_clear(g2);
        ur[p] = 1; uc[q] = 1; r++;
    }
    for (int i = 0; i < ntol; i++) { int ri = rk[i] < 0 ? r : rk[i]; if (ri > maxr[lev][i]) maxr[lev][i] = ri; storeH[lev][i] += (long long) ri * (m + k); }
    nblk[lev]++; bsz[lev] = m;
    _arb_vec_clear(P1, m); _arb_vec_clear(P2, m); _arb_vec_clear(Q1, k); _arb_vec_clear(Q2, k); _arb_vec_clear(yr, m); _arb_vec_clear(yc, k);
    _arb_vec_clear(col, m); _arb_vec_clear(row, k); free(ur); free(uc); arb_clear(t); arb_clear(e); arb_clear(piv); arb_clear(ip);
}
static slong leafsz; static long long dense_leaf = 0; static int maxlev = 0;
static void rec(slong i0, slong m, int lev) {
    if (m <= leafsz) { dense_leaf += (long long) m * m; if (lev > maxlev) maxlev = lev; return; }
    slong h = m / 2; block(i0, h, i0 + h, m - h, lev); rec(i0, h, lev + 1); rec(i0 + h, m - h, lev + 1);
}
int main(int argc, char **argv) {
    FILE *fg = fopen(argv[1], "r"); prec = atol(argv[2]); leafsz = atol(argv[3]); ntol = argc - 4; tol = malloc(ntol * sizeof(double));
    for (int i = 0; i < ntol; i++) tol[i] = atof(argv[4 + i]);
    slong nU; if (fscanf(fg, "%ld %ld", &n, &nU) != 2) return 2;
    y = _arb_vec_init(n); u = _arb_vec_init(n); w = _arb_vec_init(n); d = _arb_vec_init(n); yd = malloc(n * sizeof(double)); ord = malloc(n * sizeof(slong));
    arb_t dl; arb_init(dl);
    for (slong a = 0; a < n; a++) {
        tok(fg); tok(fg);
        arb_ptr dst[5] = {y + a, u + a, w + a, d + a, dl};
        for (int q = 0; q < 5; q++) { char *s = (char *) tok(fg); for (char *c = s; *c; c++) if (*c == '_') *c = ' '; arb_set_str(dst[q], s, prec); }
        arb_add(d + a, d + a, dl, prec);
        double m = fabs(getd(d + a)); slong e = -(slong) lround(0.5 * log2(m));
        arb_mul_2exp_si(u + a, u + a, e); arb_mul_2exp_si(w + a, w + a, e); arb_mul_2exp_si(d + a, d + a, 2 * e);
        yd[a] = getd(y + a); ord[a] = a;
    }
    scale = 0; for (slong a = 0; a < n; a++) { double m = fabs(getd(d + a)); if (m > scale) scale = m; }
    qsort(ord, n, sizeof(slong), cmpy);
    rec(0, n, 0);
    printf("{\"n\": %ld, \"leaf\": %ld, \"scale\": %.3e, \"tols\": [", n, leafsz, scale); for (int i = 0; i < ntol; i++) printf("%s%.0e", i ? ", " : "", tol[i]); printf("], \"levels\": [\n");
    long long tot[32] = {0};
    for (int l = 0; l < 64 && nblk[l]; l++) {
        printf("  {\"level\": %d, \"blocks\": %d, \"size\": %ld, \"max_rank\": [", l, nblk[l], bsz[l]);
        for (int i = 0; i < ntol; i++) printf("%s%d", i ? ", " : "", maxr[l][i]); printf("], \"stored\": [");
        for (int i = 0; i < ntol; i++) { printf("%s%lld", i ? ", " : "", storeH[l][i]); tot[i] += storeH[l][i]; } printf("]},\n");
    }
    printf("  {\"dense_leaf\": %lld, \"hodlr_total\": [", dense_leaf); for (int i = 0; i < ntol; i++) printf("%s%lld", i ? ", " : "", tot[i] + dense_leaf);
    printf("], \"dense_lower\": %lld}]}\n", (long long) n * (n + 1) / 2);
    return 0;
}
