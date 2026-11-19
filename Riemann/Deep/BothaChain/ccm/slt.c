// slt.c (explore28): Chekmachev's Woodbury/Cholesky step as an N-body problem.
// J''(nu) = G + diag(delta(nu)), G_ab = (u_a w_b - w_a u_b)/(y_a - y_b) for a != b: the P0^{-1}-Gram of all columns is a
// Cauchy-like ("particle") matrix with points y_a and two charges (u_a, w_a); it is never stored.
// Inertia(J'') = (nU, nB)  <=>  Id - W - nu Id > 0  (Haynsworth; delta_B scaled by 1 - nu).
// Structured LDL^T (Schur/GKO generator recursion; the diagonal is tracked explicitly) in Arb balls:
//   pivot p: G_ap = (u_a w_p - w_a u_p)/(y_a - y_p), g = G_ap/d_p, u_a -= g u_p, w_a -= g w_p, d_a -= g G_ap  (a > p).
// O(n) memory, O(n^2) time; every pivot is an Arb ball, so the sign count is rigorous (Sylvester's law of inertia).
// Usage: slt <base> <prec> [-g gendump] [-l gendump-in] nu1 [nu2 ...]     reads base.idx, base.src (prep_slt.py)
// CCM (ccmfmm, defl sector): env SLT_BRD=<file> adds one dense U index g (the deflation column M g g^T of ../WeilSimple):
//   file: "n J''_gg delta_g" then n lines r_a = <g, P0^{-1} c_a> (particle order of base.idx). g is eliminated FIRST; its Schur
//   complement is Cauchy-like of displacement rank 4: G'_ab = [(u_a w_b - w_a u_b) + (u2_a w2_b - w2_a u2_b)]/(y_a - y_b),
//   u2 = r, w2 = y r/J''_gg (DERIVATION.md, Lemma P'); the recursion and Lemma B run on both generator pairs. Lines marked CCM.
#include <flint/arb.h>
#include <flint/fmpz.h>
#include <flint/fmpz_vec.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

static char buf[1 << 16];
static const char *tok(FILE *f) { if (fscanf(f, "%65535s", buf) != 1) { fprintf(stderr, "read error\n"); exit(2); } return buf; }
static void rd(arb_t x, FILE *f, slong prec) { const char *s = tok(f); if (arb_set_str(x, s, prec)) { fprintf(stderr, "bad number %s\n", s); exit(2); } }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }

static arb_ptr uu, ww, dd, yy; static int *bb; static double *hh; static int PIVOT = 1, NORM = 1, EQUIL = 1, POINT = 0;
static int BRD = 0; static arb_ptr u2, w2, rr;   // CCM: second generator pair and the border row of g
static void ws_init(slong n) { uu = _arb_vec_init(n); ww = _arb_vec_init(n); dd = _arb_vec_init(n); yy = _arb_vec_init(n); bb = malloc(n * sizeof(int)); hh = malloc(n * sizeof(double));
    u2 = _arb_vec_init(n); w2 = _arb_vec_init(n); rr = _arb_vec_init(n); }   // CCM
static void swp(slong i, slong j) { if (i == j) return; arb_swap(uu + i, uu + j); arb_swap(ww + i, ww + j); arb_swap(dd + i, dd + j); arb_swap(yy + i, yy + j); int t = bb[i]; bb[i] = bb[j]; bb[j] = t; double h = hh[i]; hh[i] = hh[j]; hh[j] = h;
    arb_swap(u2 + i, u2 + j); arb_swap(w2 + i, w2 + j); }   // CCM
static double radd(const arb_t x) { return mag_get_d(arb_radref(x)); }
static double absd(const arb_t x) { return fabs(arf_get_d(arb_midref(x), ARF_RND_NEAR)) * (1 + 1e-15) + 1e-300; }
// POINT mode (backward-error certificate, after Wiktorowicz's structured LDL): every step starts from exact points; the radii after
// the step bound the local error E_p (the computed Schur complement minus the exact one of the point matrix); with H_a = sum_b 1/|y_a - y_b|,
// ||E_p||_inf <= max_a H_a (ru_a (W + RW) + |u_a| RW + rw_a (U + RU) + |w_a| RU) + rd_a.  Then A_0 = J_s + E_tot with ||E_tot|| <= eta has the
// computed pivots exactly (Sylvester), and #pos = nU gives lambda_min(Id - W) >= nu - kappa eta (see REPORT).
static double errnorm(slong p0, slong n) {
    double U = 0, W = 0, RU = 0, RW = 0, mx = 0, U2 = 0, W2 = 0, RU2 = 0, RW2 = 0;
    for (slong a = p0; a < n; a++) { double x = absd(uu + a), z = absd(ww + a), rx = radd(uu + a), rz = radd(ww + a); if (x > U) U = x; if (z > W) W = z; if (rx > RU) RU = rx; if (rz > RW) RW = rz; }
    if (BRD) for (slong a = p0; a < n; a++) { double x = absd(u2 + a), z = absd(w2 + a), rx = radd(u2 + a), rz = radd(w2 + a); if (x > U2) U2 = x; if (z > W2) W2 = z; if (rx > RU2) RU2 = rx; if (rz > RW2) RW2 = rz; }   // CCM
    for (slong a = p0; a < n; a++) {
        double v = hh[a] * (radd(uu + a) * (W + RW) + absd(uu + a) * RW + radd(ww + a) * (U + RU) + absd(ww + a) * RU) + radd(dd + a);
        if (BRD) v += hh[a] * (radd(u2 + a) * (W2 + RW2) + absd(u2 + a) * RW2 + radd(w2 + a) * (U2 + RU2) + absd(w2 + a) * RU2);   // CCM: pair 2, same bound
        if (v > mx) mx = v;
    }
    return mx * (1 + 1e-12);
}
static void recenter(slong p0, slong n) { for (slong a = p0; a < n; a++) { arb_get_mid_arb(uu + a, uu + a); arb_get_mid_arb(ww + a, ww + a); arb_get_mid_arb(dd + a, dd + a);
    if (BRD) { arb_get_mid_arb(u2 + a, u2 + a); arb_get_mid_arb(w2 + a, w2 + a); } } }   // CCM
// generator balancing on the remaining indices [p, n): exact det-1 maps only (an exact dyadic shear w -= x u, then u *= 2^e, w *= 2^-e),
// so R = u w^T - w u^T is unchanged exactly; the map is chosen from double midpoints (rigour does not depend on the choice). Gu's
// observation: GKO-type recursions lose accuracy when the generator pair becomes ill-conditioned (u nearly parallel to w).
static void balance_pair(arb_ptr uu, arb_ptr ww, slong p, slong n, slong prec);   // CCM: the pair is a parameter
static void balance(slong p, slong n, slong prec) { balance_pair(uu, ww, p, n, prec); if (BRD) balance_pair(u2, w2, p, n, prec); }
static void balance_pair(arb_ptr uu, arb_ptr ww, slong p, slong n, slong prec) {
    double A = 0, B = 0, C = 0;
    for (slong a = p; a < n; a++) { double x = arf_get_d(arb_midref(uu + a), ARF_RND_NEAR), z = arf_get_d(arb_midref(ww + a), ARF_RND_NEAR); A += x * x; B += x * z; C += z * z; }
    if (!(A > 0) || !(C > 0) || !isfinite(A) || !isfinite(B) || !isfinite(C)) return;
    double x = B / A; arb_t t; arb_init(t);
    if (fabs(x) * sqrt(A / C) > 1e-3) {         // shear w -= x u (x exact double)
        arb_set_d(t, x);
        for (slong a = p; a < n; a++) arb_submul(ww + a, t, uu + a, prec);
        C = C - 2 * x * B + x * x * A; if (C < 0) C = 0;
    }
    if (C > 0) {
        int e = (int) lround(0.25 * log2(C / A));   // u *= 2^e, w *= 2^-e makes |u| ~ |w|
        if (e != 0) for (slong a = p; a < n; a++) { arb_mul_2exp_si(uu + a, uu + a, e); arb_mul_2exp_si(ww + a, ww + a, -e); }
    }
    arb_clear(t);
}
// returns 1 if inertia (nU, nB) (all U pivots > 0, all B pivots < 0 when PIVOT <= 1; #pos = nU and #neg = nB otherwise),
// 0 if the sign count is rigorously wrong, -1 if some pivot ball contains 0
static arb_ptr RB = NULL; static arb_t DGG, DELG;   // CCM: border r_a, J''_gg, delta_g = 1/M (read in main)
static int ldl(slong n, slong nU, arb_srcptr y, arb_srcptr u, arb_srcptr w, arb_srcptr dg, arb_srcptr del, const int *isB, const arb_t nu, slong prec, const char *label, int verbose) {
    double ts = now();
    static slong *eex = NULL; if (!eex) eex = malloc(n * sizeof(slong));
    arb_t om, inv, t1, G, g, piv, bmax; arb_init(om); arb_init(inv); arb_init(t1); arb_init(G); arb_init(g); arb_init(piv); arb_init(bmax);
    arb_sub_ui(om, nu, 1, prec); arb_neg(om, om);            // 1 - nu
    for (slong a = 0; a < n; a++) {
        arb_set(uu + a, u + a); arb_set(ww + a, w + a); arb_set(yy + a, y + a); bb[a] = isB[a];
        if (isB[a]) arb_mul(t1, del + a, om, prec); else arb_set(t1, del + a);
        arb_add(dd + a, dg + a, t1, prec);
        eex[a] = 0;
        if (EQUIL) {   // exact congruence by 2^e (Jacobi scaling): inertia and the Cauchy-like form are unchanged
            double m = fabs(arf_get_d(arb_midref(dd + a), ARF_RND_NEAR));
            if (m > 0 && isfinite(m)) { slong e = -(slong) lround(0.5 * log2(m)); eex[a] = e; arb_mul_2exp_si(uu + a, uu + a, e); arb_mul_2exp_si(ww + a, ww + a, e); arb_mul_2exp_si(dd + a, dd + a, 2 * e); }
        }
    }
    // CCM: border scaling e_g (exact 2^e congruence), pair 2 zero until g is eliminated
    slong eg = 0; arb_t dG; arb_init(dG);
    for (slong a = 0; a < n; a++) { arb_zero(u2 + a); arb_zero(w2 + a); }
    if (BRD) {
        double m = fabs(arf_get_d(arb_midref(DGG), ARF_RND_NEAR)); if (EQUIL && m > 0 && isfinite(m)) eg = -(slong) lround(0.5 * log2(m));
        arb_mul_2exp_si(dG, DGG, 2 * eg);
        for (slong a = 0; a < n; a++) arb_mul_2exp_si(rr + a, RB + a, eex[a] + eg);
    }
    double eta = 0, kappa = 0, sigma = 1e300, smax = 0, YF2 = 0;
    // |y_a - y_b| lower bounds: exact integer differences (fmpz_get_d rounds toward zero) when both points are integers, else double margin
    static fmpz *yz = NULL; static char *yint = NULL;
    if (POINT && !yz) { yz = _fmpz_vec_init(n); yint = malloc(n); for (slong a = 0; a < n; a++) yint[a] = arb_get_unique_fmpz(yz + a, y + a) && arb_is_exact(y + a); }
    #define DYLO(a, b, out) do { if (yint[a] && yint[b]) { fmpz_sub(ftmp, yz + (a), yz + (b)); out = fabs(fmpz_get_d(ftmp)); }         else { out = fabs(arf_get_d(arb_midref(y + (a)), ARF_RND_NEAR) - arf_get_d(arb_midref(y + (b)), ARF_RND_NEAR)) * (1 - 1e-12) - 1e-12; } } while (0)
    fmpz_t ftmp; fmpz_init(ftmp);
    if (POINT) {
        for (slong a = 0; a < n; a++) { double s = 0, dyy;
            for (slong b = 0; b < n; b++) if (b != a) { DYLO(a, b, dyy); s += 1.0 / dyy; }
            hh[a] = s * (1 + 1e-10); }
        // kappa: sigma <= lambda_min(S_s) (G_UU is a Gram matrix, so S_s >= diag(delta_U) scaled); smax = max_B s_a 2^{-2e}; ||Y_s||_F^2
        for (slong a = 0; a < n; a++) {
            double dl = arf_get_d(arb_midref(del + a), ARF_RND_NEAR), rdl = radd(del + a);
            if (!isB[a]) { double sg = (dl - rdl) * ldexp(1.0, 2 * eex[a]) * (1 - 1e-12); if (sg < sigma) sigma = sg; }
            else { double sa = 1.0 / ((fabs(dl) - rdl) * (1 - 1e-12)) * ldexp(1.0, -2 * eex[a]); if (sa > smax) smax = sa; }
        }
        for (slong a = 0; a < n; a++) if (!isB[a]) { double ua = absd(uu + a) + radd(uu + a), wa = absd(ww + a) + radd(ww + a), dyy;
            for (slong b = 0; b < n; b++) if (isB[b]) { DYLO(a, b, dyy); double v = (ua * (absd(ww + b) + radd(ww + b)) + wa * (absd(uu + b) + radd(uu + b))) / dyy; YF2 += v * v; } }
        if (BRD) {   // CCM: g is a U index: sigma <= delta_g 2^{2 e_g}; its row of Y_s is r_B
            double dl = arf_get_d(arb_midref(DELG), ARF_RND_NEAR), rdl = radd(DELG), sg = (dl - rdl) * ldexp(1.0, 2 * eg) * (1 - 1e-12); if (sg < sigma) sigma = sg;
            for (slong b = 0; b < n; b++) if (isB[b]) { double v = absd(rr + b) + radd(rr + b); YF2 += v * v; }
        }
        YF2 *= (1 + 1e-10);
        eta = errnorm(0, n);   // input radii
        if (BRD) {   // CCM: border input radii: row a gets rad(r_a), row g gets sum_a rad(r_a) + rad(J''_gg)
            double mr = 0, sr = 0; for (slong a = 0; a < n; a++) { double r = radd(rr + a); sr += r; if (r > mr) mr = r; }
            eta += mr; double rg = (sr + radd(dG)) * (1 + 1e-12); if (rg > eta) eta = rg;
            for (slong a = 0; a < n; a++) arb_get_mid_arb(rr + a, rr + a); arb_get_mid_arb(dG, dG);
        }
        recenter(0, n);
    }
    // block order: all U indices first, then B (PIVOT <= 1); PIVOT >= 1: largest |d| among the remaining indices of the current block
    slong cnt[2][3] = {{0,0,0},{0,0,0}};    // [block][pos, neg, undetermined]
    double maxrel = 0; slong firstbad = -1, bmaxi = -1; int have = 0;
    char gpiv[64] = "";
    if (BRD) {   // CCM: eliminate g first (a U pivot): d_a -= r_a^2/d_g, pair 2 = (r, y r/d_g)
        int s = arb_is_positive(dG) ? 0 : (arb_is_negative(dG) ? 1 : 2); cnt[0][s]++; if (s == 2) firstbad = -2;
        { char *sg = arb_get_str(dG, 8, 0); snprintf(gpiv, sizeof gpiv, "%s", sg); flint_free(sg); }
        arb_inv(inv, dG, prec);
        for (slong a = 0; a < n; a++) {
            arb_mul(t1, rr + a, inv, prec);
            arb_submul(dd + a, t1, rr + a, prec);
            arb_set(u2 + a, rr + a); arb_mul(w2 + a, t1, yy + a, prec);
        }
        if (POINT) { eta += errnorm(0, n); recenter(0, n); }
    }
    for (slong p = 0; p < n; p++) {
        if (PIVOT >= 1) {
            int needB = (p >= nU && PIVOT == 1); slong q = -1; double best = -1;
            for (slong a = p; a < n; a++) {
                if (PIVOT == 1 && bb[a] != needB) continue;
                double m = fabs(arf_get_d(arb_midref(dd + a), ARF_RND_NEAR)); if (m > best) { best = m; q = a; }
            }
            if (q >= 0) swp(p, q);
        }
        if (NORM && p % NORM == 0) balance(p, n, prec);
        if (POINT && NORM && p % NORM == 0) { eta += errnorm(p, n); recenter(p, n); }   // shear rounding: a perturbation F_p supported on [p, n), zero diagonal
        arb_set(piv, dd + p);
        int b = bb[p], s = arb_is_positive(piv) ? 0 : (arb_is_negative(piv) ? 1 : 2);
        cnt[b][s]++;
        if (s == 2 && firstbad < 0) firstbad = p;
        if (!arb_contains_zero(piv)) {
            double r = mag_get_d(arb_radref(piv)), m = fabs(arf_get_d(arb_midref(piv), ARF_RND_NEAR));
            if (m > 0 && r / m > maxrel) maxrel = r / m;
        }
        if (b && (!have || arf_cmp(arb_midref(piv), arb_midref(bmax)) > 0)) { arb_set(bmax, piv); bmaxi = p; have = 1; }
        arb_inv(inv, piv, prec);
        for (slong a = p + 1; a < n; a++) {
            arb_mul(t1, uu + a, ww + p, prec); arb_submul(t1, ww + a, uu + p, prec);
            if (BRD) { arb_addmul(t1, u2 + a, w2 + p, prec); arb_submul(t1, w2 + a, u2 + p, prec); }   // CCM
            arb_sub(G, yy + a, yy + p, prec); arb_div(G, t1, G, prec);
            arb_mul(g, G, inv, prec);
            arb_submul(uu + a, g, uu + p, prec); arb_submul(ww + a, g, ww + p, prec); arb_submul(dd + a, g, G, prec);
            if (BRD) { arb_submul(u2 + a, g, u2 + p, prec); arb_submul(w2 + a, g, w2 + p, prec); }   // CCM
        }
        if (POINT) { eta += errnorm(p + 1, n); recenter(p + 1, n); }
    }
    if (POINT) {   // lambda_min(Id - W) >= nu - kappa eta,  kappa = (1 + ||Y_s||^2/(sigma (sigma - eta))) smax   (with ||Y_s|| <= ||Y_s||_F + eta)
        double Yn = sqrt(YF2) + eta; kappa = (1 + Yn * Yn / (sigma * (sigma - eta))) * smax * (1 + 1e-10);
    }
    int ok; slong nUe = nU + BRD, ne = n + BRD;   // CCM: g counts as a U index
    if (PIVOT <= 1) ok = (cnt[0][0] == nUe && cnt[1][1] == ne - nUe) ? 1 : ((cnt[0][2] + cnt[1][2]) > 0 ? -1 : 0);
    if (POINT) ok = ((cnt[0][2] + cnt[1][2]) > 0) ? -1 : ((cnt[0][0] + cnt[1][0] == nUe) ? 1 : 0);
    else ok = ((cnt[0][2] + cnt[1][2]) > 0) ? -1 : ((cnt[0][0] + cnt[1][0] == nUe && cnt[0][1] + cnt[1][1] == ne - nUe) ? 1 : 0);
    if (verbose) {
        char *sb = have ? arb_get_str(bmax, 15, 0) : NULL;
        printf("{\"nu\": \"%s\", \"n\": %ld, \"nU\": %ld, \"U\": [%ld, %ld, %ld], \"B\": [%ld, %ld, %ld], \"inertia_ok\": %s, \"top_B_pivot\": \"%s\", \"top_B_index\": %ld, \"max_rel_rad\": %.3e, \"first_undetermined\": %ld, \"prec\": %ld, \"pivot\": %d, \"norm\": %d, \"point\": %d, \"eta\": %.3e, \"kappa\": %.3e, \"kappa_eta\": %.3e, \"sigma\": %.3e, \"border\": %d, \"g_pivot\": \"%s\", \"time\": %.1f}\n",
               label, ne, nUe, cnt[0][0], cnt[0][1], cnt[0][2], cnt[1][0], cnt[1][1], cnt[1][2], ok == 1 ? "true" : "false", sb ? sb : "", bmaxi, maxrel, firstbad, prec, PIVOT, NORM, POINT, eta, kappa, kappa * eta, sigma, BRD, gpiv, now() - ts);
        fflush(stdout);
        if (sb) flint_free(sb);
    }
    arb_clear(om); arb_clear(inv); arb_clear(t1); arb_clear(G); arb_clear(g); arb_clear(piv); arb_clear(bmax); arb_clear(dG);
    return ok;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: slt base prec [-g out] [-l in] nu...\n"); return 1; }
    const char *base = argv[1]; slong prec = atol(argv[2]); int ai = 3;
    const char *gout = NULL, *gin = NULL;
    if (getenv("SLT_PIVOT")) PIVOT = atoi(getenv("SLT_PIVOT"));
    if (getenv("SLT_NORM")) NORM = atoi(getenv("SLT_NORM"));
    if (getenv("SLT_EQUIL")) EQUIL = atoi(getenv("SLT_EQUIL"));
    if (getenv("SLT_POINT")) POINT = atoi(getenv("SLT_POINT"));
    while (ai < argc && argv[ai][0] == '-' && (argv[ai][1] == 'g' || argv[ai][1] == 'l')) { if (argv[ai][1] == 'g') gout = argv[ai+1]; else gin = argv[ai+1]; ai += 2; }
    char fn[4096]; double t0 = now();
    slong n, nU, nS;
    arb_ptr y, u, w, dg, del; int *isB; char *ty;
    if (!gin) {
        snprintf(fn, sizeof fn, "%s.idx", base); FILE *fi = fopen(fn, "r"); if (!fi) { perror(fn); return 1; }
        if (fscanf(fi, "%ld %ld %ld", &n, &nU, &nS) != 3) return 2;
        y = _arb_vec_init(n); u = _arb_vec_init(n); w = _arb_vec_init(n); dg = _arb_vec_init(n); del = _arb_vec_init(n);
        arb_ptr be = _arb_vec_init(n), kp = _arb_vec_init(n);
        isB = malloc(n * sizeof(int)); ty = malloc(n);
        for (slong a = 0; a < n; a++) {
            ty[a] = tok(fi)[0]; isB[a] = atoi(tok(fi));
            rd(y + a, fi, prec);
            if (ty[a] == 'N') { rd(be + a, fi, prec); rd(dg + a, fi, prec); rd(del + a, fi, prec); }
            else if (ty[a] == 'V') { rd(be + a, fi, prec); rd(kp + a, fi, prec); rd(dg + a, fi, prec); rd(del + a, fi, prec); }
            else if (ty[a] == 'P') { rd(w + a, fi, prec); rd(dg + a, fi, prec); rd(del + a, fi, prec); }
            else { fprintf(stderr, "bad type\n"); return 2; }
        }
        fclose(fi);
        snprintf(fn, sizeof fn, "%s.src", base); FILE *fs = fopen(fn, "r"); if (!fs) { perror(fn); return 1; }
        // Cauchy sums F(y) = sum_j gam_j/(y - x_j), F2(y) = sum_j gam_j/(y - x_j)^2. x_j = 4 j^2 and y = k^2 or 4 j^2 are exact
        // integers below 2^68, so each term is a fixed-point quotient: Gq_j = gam_j 2^P (rounded, error <= 1 + rad(gam_j) 2^P units),
        // q = tdiv(Gq_j, D), q2 = tdiv(q, D) with |D| = |y - x_j| >= 1: each truncation errs by < 1 unit, so
        // |F 2^P - sum q| <= sum_j (1 + rad_j 2^P + 1) and the same for F2. Error added as a ball radius: rigorous.
        slong P = prec + 64;
        fmpz *X = _fmpz_vec_init(nS), *Gq = _fmpz_vec_init(nS);
        double *rj = malloc(nS * sizeof(double));   // upper bounds of rad(gam_j)
        mag_t radsum, radmax, mt; mag_init(radsum); mag_init(radmax); mag_init(mt);
        arb_t gtmp; arb_init(gtmp);
        for (slong s = 0; s < nS; s++) {
            if (fmpz_set_str(X + s, tok(fs), 10)) { fprintf(stderr, "bad x\n"); return 2; }
            rd(gtmp, fs, prec);
            arf_get_fmpz_fixed_si(Gq + s, arb_midref(gtmp), -P);
            mag_add(radsum, radsum, arb_radref(gtmp)); if (mag_cmp(arb_radref(gtmp), radmax) > 0) mag_set(radmax, arb_radref(gtmp));
            { mag_t mm; mag_init(mm); mag_mul_2exp_si(mm, arb_radref(gtmp), P); rj[s] = mag_get_d(mm) * (1 + 1e-9); mag_clear(mm); }   // rad_j 2^P, upper bound
        }
        fclose(fs); arb_clear(gtmp);
        // radius: sum_j rad_j + 3 nS 2^-P  (units: rounding of Gq (<= 1), two truncations (< 1 each))
        arb_t Erad; arb_init(Erad); arb_set_ui(Erad, 3 * (ulong) nS + 3); arb_mul_2exp_si(Erad, Erad, -P);
        { arb_t t; arb_init(t); arf_set_mag(arb_midref(t), radsum); arb_add(Erad, Erad, t, prec); arb_clear(t); }
        double t1 = now();
        arb_t F, F2, t, e; arb_init(F); arb_init(F2); arb_init(t); arb_init(e);
        fmpz_t Y, D, q, q2, Fs, F2s; fmpz_init(Y); fmpz_init(D); fmpz_init(q); fmpz_init(q2); fmpz_init(Fs); fmpz_init(F2s);
        slong nt = 0;
        // chunking: SLT_RANGE="a0 a1" + SLT_FOUT=file computes F, F2 for indices in [a0, a1) only and writes them; SLT_FIN=f1,f2,... reads them
        slong r0 = 0, r1 = n; FILE *fout = NULL; char *fin = getenv("SLT_FIN");
        if (getenv("SLT_RANGE")) { sscanf(getenv("SLT_RANGE"), "%ld %ld", &r0, &r1); fout = fopen(getenv("SLT_FOUT"), "w"); }
        arb_ptr Fv = NULL, F2v = NULL; char *have = NULL;
        if (fin) {
            Fv = _arb_vec_init(n); F2v = _arb_vec_init(n); have = calloc(n, 1);
            char *list = strdup(fin);
            for (char *f = strtok(list, ","); f; f = strtok(NULL, ",")) {
                FILE *ff = fopen(f, "r"); if (!ff) { perror(f); return 1; }
                slong a; while (fscanf(ff, "%ld", &a) == 1) { char *s1 = (char *) tok(ff); for (char *c = s1; *c; c++) if (*c == '_') *c = ' '; arb_set_str(Fv + a, s1, prec);
                    s1 = (char *) tok(ff); for (char *c = s1; *c; c++) if (*c == '_') *c = ' '; arb_set_str(F2v + a, s1, prec); have[a] = 1; }
                fclose(ff);
            }
        }
        for (slong a = 0; a < n; a++) {
            if (ty[a] == 'P') { arb_zero(u + a); continue; }
            if (fin) { if (!have[a]) { fprintf(stderr, "missing F for %ld\n", a); return 2; } arb_set(F, Fv + a); arb_set(F2, F2v + a); }
            else if (a < r0 || a >= r1) continue;
            else if (ty[a] == 'N') {
                if (!arb_get_unique_fmpz(Y, y + a)) { fprintf(stderr, "y not integer\n"); return 2; }
                fmpz_zero(Fs); fmpz_zero(F2s);
                double es = 0;   // sum_j rad_j/|D_j| (|D_j| as a double lower bound)
                for (slong s = 0; s < nS; s++) {
                    fmpz_sub(D, Y, X + s);
                    fmpz_tdiv_q(q, Gq + s, D); fmpz_add(Fs, Fs, q);
                    fmpz_tdiv_q(q2, q, D); fmpz_add(F2s, F2s, q2);
                    es += rj[s] / (fabs(fmpz_get_d(D)) * (1 - 1e-12));      // units of 2^-P
                }
                // |F - Fs 2^-P| <= sum rad_j/|D_j| + (2 nS + 1) 2^-P ; |F2 - F2s 2^-P| <= sum rad_j/|D_j| + (3 nS + 1) 2^-P  (|D_j| >= 1)
                arb_set_d(e, es * (1 + 1e-9) + 3.0 * (double) nS + 3.0); arb_mul_2exp_si(e, e, -P);   // (sum rad_j 2^P/|D_j| + 3 nS + 3) 2^-P
                arb_set_fmpz(F, Fs); arb_mul_2exp_si(F, F, -P); arb_add_error(F, e);
                arb_set_fmpz(F2, F2s); arb_mul_2exp_si(F2, F2, -P); arb_add_error(F2, e);
            } else {   // pole (y_v not an integer): balls, gam_j in Gq_j 2^-P +- (2^-P + radmax)
                arb_zero(F); arb_zero(F2);
                arb_t gj, er; arb_init(gj); arb_init(er); arb_one(er); arb_mul_2exp_si(er, er, -P);
                { arb_t tt; arb_init(tt); arf_set_mag(arb_midref(tt), radmax); arb_add(er, er, tt, prec); arb_clear(tt); }
                for (slong s = 0; s < nS; s++) {
                    arb_set_fmpz(gj, Gq + s); arb_mul_2exp_si(gj, gj, -P); arb_add_error(gj, er);
                    arb_set_fmpz(t, X + s); arb_sub(t, y + a, t, prec); arb_inv(t, t, prec);
                    arb_mul(e, gj, t, prec); arb_add(F, F, e, prec); arb_addmul(F2, e, t, prec);
                }
                arb_clear(gj); arb_clear(er);
            }
            if (fout) {
                char *s1 = arb_get_str(F, prec * 0.30103 + 12, ARB_STR_MORE), *s2 = arb_get_str(F2, prec * 0.30103 + 12, ARB_STR_MORE);
                for (char *c = s1; *c; c++) if (*c == ' ') *c = '_'; for (char *c = s2; *c; c++) if (*c == ' ') *c = '_';
                fprintf(fout, "%ld %s %s\n", a, s1, s2); flint_free(s1); flint_free(s2); nt++;
                if (nt % 1000 == 0) { fprintf(stderr, "  %ld targets, %.0fs\n", nt, now() - t1); fflush(fout); }
                continue;
            }
            arb_set(u + a, be + a);
            arb_mul(w + a, be + a, F, prec);
            if (ty[a] == 'V') arb_add(w + a, w + a, kp + a, prec);
            arb_mul(t, be + a, be + a, prec); arb_addmul(dg + a, t, F2, prec);
            nt++;
            if (nt % 1000 == 0) { fprintf(stderr, "  %ld targets, %.0fs\n", nt, now() - t1); }
        }
        fprintf(stderr, "read %.1fs, Cauchy sums %ld x %ld in %.1fs\n", t1 - t0, nt, nS, now() - t1);
        if (fout) { fclose(fout); return 0; }
        _fmpz_vec_clear(X, nS); _fmpz_vec_clear(Gq, nS); free(rj); _arb_vec_clear(be, n); _arb_vec_clear(kp, n);
        arb_clear(F); arb_clear(F2); arb_clear(t); arb_clear(e); arb_clear(Erad);
        fmpz_clear(Y); fmpz_clear(D); fmpz_clear(q); fmpz_clear(q2); fmpz_clear(Fs); fmpz_clear(F2s);
        if (gout) {
            FILE *fg = fopen(gout, "w"); fprintf(fg, "%ld %ld\n", n, nU);
            for (slong a = 0; a < n; a++) {
                char *s1 = arb_get_str(y + a, prec * 0.30103 + 12, ARB_STR_MORE), *s2 = arb_get_str(u + a, prec * 0.30103 + 12, ARB_STR_MORE),
                     *s3 = arb_get_str(w + a, prec * 0.30103 + 12, ARB_STR_MORE), *s4 = arb_get_str(dg + a, prec * 0.30103 + 12, ARB_STR_MORE),
                     *s5 = arb_get_str(del + a, prec * 0.30103 + 12, ARB_STR_MORE);
                #define NS(s) for (char *q = s; *q; q++) if (*q == ' ') *q = '_';
                NS(s1) NS(s2) NS(s3) NS(s4) NS(s5)
                fprintf(fg, "%c %d %s %s %s %s %s\n", ty[a], isB[a], s1, s2, s3, s4, s5);
                flint_free(s1); flint_free(s2); flint_free(s3); flint_free(s4); flint_free(s5);
            }
            fclose(fg);
        }
    } else {
        FILE *fg = fopen(gin, "r"); if (!fg) { perror(gin); return 1; }
        if (fscanf(fg, "%ld %ld", &n, &nU) != 2) return 2;
        y = _arb_vec_init(n); u = _arb_vec_init(n); w = _arb_vec_init(n); dg = _arb_vec_init(n); del = _arb_vec_init(n);
        isB = malloc(n * sizeof(int)); ty = malloc(n);
        for (slong a = 0; a < n; a++) {
            ty[a] = tok(fg)[0]; isB[a] = atoi(tok(fg));
            arb_ptr dst[5] = {y + a, u + a, w + a, dg + a, del + a};
            for (int q = 0; q < 5; q++) { char *s = (char *) tok(fg); for (char *c = s; *c; c++) if (*c == '_') *c = ' '; if (arb_set_str(dst[q], s, prec)) { fprintf(stderr, "bad %s\n", s); return 2; } }
        }
        fclose(fg);
    }
    if (getenv("SLT_BRD") && getenv("SLT_BRD")[0]) {   // CCM: border file of the defl sector
        FILE *fb = fopen(getenv("SLT_BRD"), "r"); if (!fb) { perror(getenv("SLT_BRD")); return 1; }
        slong nb; if (fscanf(fb, "%ld", &nb) != 1 || nb != n) { fprintf(stderr, "border size mismatch\n"); return 2; }
        RB = _arb_vec_init(n); arb_init(DGG); arb_init(DELG);
        { char *s = (char *) tok(fb); for (char *c = s; *c; c++) if (*c == '_') *c = ' '; if (arb_set_str(DGG, s, prec)) { fprintf(stderr, "bad Jgg\n"); return 2; } }
        { char *s = (char *) tok(fb); for (char *c = s; *c; c++) if (*c == '_') *c = ' '; if (arb_set_str(DELG, s, prec)) { fprintf(stderr, "bad delta_g\n"); return 2; } }
        for (slong a = 0; a < n; a++) { char *s = (char *) tok(fb); for (char *c = s; *c; c++) if (*c == '_') *c = ' '; if (arb_set_str(RB + a, s, prec)) { fprintf(stderr, "bad r %s\n", s); return 2; } }
        fclose(fb); BRD = 1;
    }
    // structured LDL^T for each shift (or a bisection for lambda_min(Id - W): -b lo hi steps)
    ws_init(n);
    if (ai < argc && strcmp(argv[ai], "-b") == 0) {
        arb_t lo, hi, mid; arb_init(lo); arb_init(hi); arb_init(mid);
        arb_set_str(lo, argv[ai+1], prec); arb_set_str(hi, argv[ai+2], prec); int steps = atoi(argv[ai+3]);
        int okl = ldl(n, nU, y, u, w, dg, del, isB, lo, prec, argv[ai+1], 1), okh = ldl(n, nU, y, u, w, dg, del, isB, hi, prec, argv[ai+2], 1);
        if (!(okl == 1 && okh == 0)) { printf("{\"bisect\": \"bad bracket\", \"lo_ok\": %d, \"hi_ok\": %d}\n", okl, okh); return 3; }
        for (int it = 0; it < steps; it++) {
            arb_add(mid, lo, hi, prec); arb_mul_2exp_si(mid, mid, -1); arb_set_round(mid, mid, 64); arb_get_mid_arb(mid, mid);
            int ok = ldl(n, nU, y, u, w, dg, del, isB, mid, prec, "mid", 0);
            if (ok == 1) arb_set(lo, mid); else if (ok == 0) arb_set(hi, mid); else { printf("{\"bisect\": \"undetermined at step %d\"}\n", it); break; }
        }
        char *sl = arb_get_str(lo, 15, 0), *sh = arb_get_str(hi, 15, 0);
        printf("{\"lambda_min_bracket\": [\"%s\", \"%s\"], \"meaning\": \"inertia ok at lo (lambda_min(Id-W) > lo), fails at hi (lambda_min < hi)\"}\n", sl, sh);
        flint_free(sl); flint_free(sh);
    } else {
        arb_t nu; arb_init(nu);
        for (; ai < argc; ai++) {
            if (arb_set_str(nu, argv[ai], prec)) { fprintf(stderr, "bad nu\n"); return 2; }
            ldl(n, nU, y, u, w, dg, del, isB, nu, prec, argv[ai], 1);
        }
    }
    fprintf(stderr, "total %.1fs\n", now() - t0);
    return 0;
}
