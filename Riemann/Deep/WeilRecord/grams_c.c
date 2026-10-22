/* grams_c.c -- C port of bc.py grams(): rigorous raw Gram sums for A0 = E*C_pE and B0 = E*C_p^2E.
   Same mathematics and the same ball-arithmetic rigour as the Python version (Arb, now part of FLINT):
   on each piece of [0,1] between the breakpoints |1 - t_n| the integrands are polynomials of degree <= 2N-2,
   so (N+1)-point Gauss is exact; nodes, weights and Legendre values are Arb balls.  Legendre recurrence
   in ball arithmetic at bprec = prec + ceil(1.28 N) + 64 bits (wrapping growth <= (1+sqrt 2)^j), in the
   factorial-scaled form R_j = j! P_j: R_{j+1} = (2j+1) y R_j - j^2 R_{j-1} (integer coefficients, no divisions),
   P_j = R_j * (1/j!) with 1/j! a ball at prec + 64; shifted sums accumulated at prec + 64, stored at prec.
   Output: bc.py binary format (JSON header + length-prefixed integers), keys A0raw_p, B0raw_p (unscaled,
   unsymmetrised partial sums over the pieces ip with ip % K == k), N, prec, part -> finish with 'bc.py merge'.
   Build: cc -O2 -o grams_c grams_c.c -I/opt/homebrew/include -L/opt/homebrew/lib -lflint -lgmp -lmpfr
   Usage: grams_c TAG N PREC k/K OUTFILE   TAG = e3 | <integer mu> | L<p>o<q> (half-width p/q)
   2026-10-21 (weilrec4): env GRAMS_PIECES="i,j,..." selects an explicit set of pieces (load balancing across
   parallel processes) instead of ip % K == k; the header records "tag", "npieces" and the processed "pieces",
   so a merger can check that the parts cover every piece exactly once (cert4.py does).                    */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <flint/flint.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/arb_mat.h>
#include <flint/arb_hypgeom.h>
#include <flint/acb.h>

/* MODE "bound" (optional 6th argument): per piece, n-point Gauss with n chosen from a rigorous error bound instead of
   the exact (N+1)-point rule.  Bernstein-Walsh: |P_j(z)| <= rho(z)^j, rho(z) = a + sqrt(a^2-1), a(z) = (|z-1|+|z+1|)/2
   (P_j is bounded by 1 on [-1,1]).  a is convex, so its max over a (filled) ellipse is at most its max over the
   vertices of a circumscribed K-gon.  Gauss error for g analytic in the Bernstein ellipse E_r of the piece with
   |g| <= M there (Trefethen ATAP Thm 19.3, scaled to half-length h): h (64/15) M r^{-2n}/(r^2 - 1).  Integrands
   P_i Q_k (A0) and Q_i Q_k (B0), Q_k = sum_a c_a P_k(. + t_a): M_A = M_P M_Q, M_B = M_Q^2 with
   M_P = max rho^{N-1}, M_Q = sum_a c_a max rho(. + t_a)^{N-1}.  The bounds are added to the radius of every raw entry. */
static void logrho_max(arb_t out, const arb_t c, const arb_t hh, const arb_t r, slong K, slong prec) {
    arb_t A, B, sec, th, cs, sn, a, s, x, y, best; acb_t z, zz;
    arb_init(A); arb_init(B); arb_init(sec); arb_init(th); arb_init(cs); arb_init(sn); arb_init(a); arb_init(s);
    arb_init(x); arb_init(y); arb_init(best); acb_init(z); acb_init(zz);
    arb_inv(s, r, prec); arb_add(A, r, s, prec); arb_mul(A, A, hh, prec); arb_mul_2exp_si(A, A, -1);
    arb_sub(B, r, s, prec); arb_mul(B, B, hh, prec); arb_mul_2exp_si(B, B, -1);
    arb_const_pi(th, prec); arb_div_ui(th, th, K, prec); arb_cos(sec, th, prec); arb_inv(sec, sec, prec);
    arb_mul(A, A, sec, prec); arb_mul(B, B, sec, prec);
    arb_zero(best);
    for (slong k = 0; k < K; k++) {
        arb_const_pi(th, prec); arb_mul_ui(th, th, 2*k, prec); arb_div_ui(th, th, K, prec); arb_sin_cos(sn, cs, th, prec);
        arb_mul(x, A, cs, prec); arb_add(x, x, c, prec); arb_mul(y, B, sn, prec);
        acb_set_arb_arb(z, x, y);
        acb_sub_ui(zz, z, 1, prec); acb_abs(a, zz, prec);
        acb_add_ui(zz, z, 1, prec); acb_abs(s, zz, prec);
        arb_add(a, a, s, prec); arb_mul_2exp_si(a, a, -1);
        arb_sqr(s, a, prec); arb_sub_ui(s, s, 1, prec); arb_sqrtpos(s, s, prec); arb_add(a, a, s, prec);
        arb_log(a, a, prec);
        arb_max(best, best, a, prec);
    }
    arf_t ub; arf_init(ub); arb_get_ubound_arf(ub, best, prec); arb_set_arf(out, ub); arf_clear(ub);
    arb_clear(A); arb_clear(B); arb_clear(sec); arb_clear(th); arb_clear(cs); arb_clear(sn); arb_clear(a); arb_clear(s);
    arb_clear(x); arb_clear(y); arb_clear(best); acb_clear(z); acb_clear(zz);
}
/* returns n (n = N+1 means: exact rule, eA = eB = 0) */
static slong plan_piece(mag_t eA, mag_t eB, double *rused, const arb_t md, const arb_t hh, arb_srcptr at, arb_srcptr ac, int na, slong N, slong prec) {
    const double rs[] = {1.25, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 12, 16, 24, 32};
    slong pp = 128, K = 64, bestn = N + 1;
    arb_t r, lp, lq, MP, MQ, MA, MB, fac, tt, cc, e1, e2; arb_init(r); arb_init(lp); arb_init(lq); arb_init(MP); arb_init(MQ); arb_init(MA); arb_init(MB);
    arb_init(fac); arb_init(tt); arb_init(cc); arb_init(e1); arb_init(e2);
    mag_zero(eA); mag_zero(eB); *rused = 0;
    for (int ir = 0; ir < (int)(sizeof rs/sizeof rs[0]); ir++) {
        arb_set_d(r, rs[ir]);
        logrho_max(lp, md, hh, r, K, pp); arb_mul_ui(lp, lp, N - 1, pp); arb_exp(MP, lp, pp);
        arb_zero(MQ);
        for (int a = 0; a < na; a++) {
            arb_add(cc, md, at + a, pp); logrho_max(lq, cc, hh, r, K, pp); arb_mul_ui(lq, lq, N - 1, pp); arb_exp(tt, lq, pp);
            arb_addmul(MQ, ac + a, tt, pp);
        }
        arb_mul(MA, MP, MQ, pp); arb_sqr(MB, MQ, pp);
        arb_sqr(fac, r, pp); arb_sub_ui(fac, fac, 1, pp); arb_inv(fac, fac, pp); arb_mul(fac, fac, hh, pp);
        arb_mul_ui(fac, fac, 64, pp); arb_div_ui(fac, fac, 15, pp);
        arb_max(tt, MA, MB, pp); arb_mul(tt, tt, fac, pp); arb_log(tt, tt, pp);
        double lr = log(rs[ir]);
        double need = (arf_get_d(arb_midref(tt), ARF_RND_UP) + prec*0.6931471805599453)/(2*lr);
        slong n = (slong) ceil(need) + 1; if (n < 2) n = 2;
        if (n >= bestn) continue;
        arb_set_d(cc, rs[ir]); arb_pow_ui(cc, cc, 2*n, pp); arb_inv(cc, cc, pp); arb_mul(cc, cc, fac, pp);
        arb_mul(e1, cc, MA, pp); arb_mul(e2, cc, MB, pp);
        mag_t m1, m2; mag_init(m1); mag_init(m2); arb_get_mag(m1, e1); arb_get_mag(m2, e2);
        if (mag_cmp_2exp_si(m1, -prec + 2) > 0 || mag_cmp_2exp_si(m2, -prec + 2) > 0) { mag_clear(m1); mag_clear(m2); continue; }
        bestn = n; mag_set(eA, m1); mag_set(eB, m2); *rused = rs[ir]; mag_clear(m1); mag_clear(m2);
    }
    arb_clear(r); arb_clear(lp); arb_clear(lq); arb_clear(MP); arb_clear(MQ); arb_clear(MA); arb_clear(MB);
    arb_clear(fac); arb_clear(tt); arb_clear(cc); arb_clear(e1); arb_clear(e2);
    return bestn;
}

static double t0;
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }
static void die(const char *m) { fprintf(stderr, "FATAL: %s\n", m); fflush(stdout); exit(2); }

static long vm_int(long n) {             /* prime p if n = p^k (k >= 1), else 0 */
    for (long p = 2; p <= n; p++) if (n % p == 0) { long m = n; while (m % p == 0) m /= p; return m == 1 ? p : 0; }
    return 0;
}

/* tag: 0 = integer mu, 1 = e3 (L = 3/2), 2 = L = P/Qd */
static int kind; static long MU, LP, LQ;
static void L_at(arb_t L, slong prec) {
    if (kind == 0) { arb_log_ui(L, MU, prec); arb_mul_2exp_si(L, L, -1); }
    else if (kind == 1) { arb_set_ui(L, 3); arb_mul_2exp_si(L, L, -1); }
    else { fmpq_t q; fmpq_init(q); fmpq_set_si(q, LP, LQ); arb_set_fmpq(L, q, prec); fmpq_clear(q); }
}
/* exact rational 2b/a when n^a = mu^b (integer mu only), as in bc.py Setup.tau */
static int tau_exact(arb_t t, long n, slong prec) {
    if (kind != 0) return 0;
    fmpz_t x, y; fmpz_init(x); fmpz_init(y); int found = 0;
    for (long a = 1; a <= 64 && !found; a++) {
        fmpz_set_ui(x, n); fmpz_pow_ui(x, x, a);
        for (long b = 1; b <= 2*a; b++) {
            fmpz_set_ui(y, MU); fmpz_pow_ui(y, y, b);
            if (fmpz_equal(x, y)) { fmpq_t q; fmpq_init(q); fmpq_set_si(q, 2*b, a); arb_set_fmpq(t, q, prec); fmpq_clear(q); found = 1; break; }
            if (fmpz_cmp(y, x) > 0) break;
        }
    }
    fmpz_clear(x); fmpz_clear(y); return found;
}

/* ---- bc.py binary writer ---- */
static void wint_fmpz(FILE *f, const fmpz_t x) {
    signed char s = fmpz_sgn(x) < 0 ? -1 : 1;
    mpz_t z; mpz_init(z); fmpz_get_mpz(z, x); mpz_abs(z, z);
    size_t cnt = (mpz_sizeinbase(z, 2) + 7)/8; if (mpz_sgn(z) == 0) cnt = 1;
    unsigned char *buf = calloc(cnt, 1); size_t got = 0;
    if (mpz_sgn(z) != 0) mpz_export(buf, &got, -1, 1, 0, 0, z);
    unsigned int len = (unsigned int) cnt;
    fwrite(&s, 1, 1, f); fwrite(&len, 4, 1, f); fwrite(buf, 1, cnt, f);   /* '<bI' little-endian (arm64/x86 are LE) */
    free(buf); mpz_clear(z);
}
static void write_arb(FILE *f, const arb_t x) {
    fmpz_t m, e; fmpz_init(m); fmpz_init(e); arf_t r; arf_init(r);
    if (!arb_is_finite(x)) die("non-finite entry");
    arf_get_fmpz_2exp(m, e, arb_midref(x));
    wint_fmpz(f, m); wint_fmpz(f, e);
    arf_set_mag(r, arb_radref(x));
    if (arf_is_zero(r)) { fmpz_zero(m); fmpz_zero(e); } else arf_get_fmpz_2exp(m, e, r);
    wint_fmpz(f, m); wint_fmpz(f, e);
    fmpz_clear(m); fmpz_clear(e); arf_clear(r);
}

int main(int argc, char **argv) {
    if (argc < 6) die("usage: grams_c TAG N PREC k/K OUTFILE");
    const char *tag = argv[1]; slong N = atol(argv[2]); slong prec = atol(argv[3]);
    long pk, pK; if (sscanf(argv[4], "%ld/%ld", &pk, &pK) != 2 || pK < 1 || pk < 0 || pk >= pK) die("bad part");
    const char *out = argv[5];
    int bound_mode = (argc > 6 && strcmp(argv[6], "bound") == 0);
    if (strcmp(tag, "e3") == 0) kind = 1;
    else if (tag[0] == 'L') { kind = 2; if (sscanf(tag + 1, "%ldo%ld", &LP, &LQ) != 2 || LP <= 0 || LQ <= 0) die("bad L tag"); }
    else { kind = 0; MU = atol(tag); if (MU < 2) die("bad mu"); }
    slong bprec = prec + (slong) ceil(1.28*N) + 64, accprec = prec + 64;
    t0 = now();
    arb_t L, one, mone, tmp, tmp2; arb_init(L); arb_init(one); arb_init(mone); arb_init(tmp); arb_init(tmp2);
    L_at(L, bprec); arb_one(one); arb_set_si(mone, -1);
    /* prime powers n < mu = e^{2L} */
    long ns[256], ps[256]; int nn = 0;
    { arb_t twoL, ln; arb_init(twoL); arb_init(ln); arb_mul_2exp_si(twoL, L, 1);
      long top = (kind == 0) ? MU : (long) ceil(exp(2.0*(kind == 1 ? 1.5 : (double)LP/LQ))) + 1;
      for (long n = 2; n <= top; n++) {
          long p = vm_int(n); if (!p) continue;
          int ok = 0;
          if (kind == 0) ok = n < MU;
          else { arb_log_ui(ln, n, bprec); if (arb_lt(ln, twoL)) ok = 1; else if (arb_gt(ln, twoL)) ok = 0; else die("cannot decide n < mu"); }
          if (ok) { ns[nn] = n; ps[nn] = p; nn++; }
      }
      arb_clear(twoL); arb_clear(ln); }
    arb_ptr tau = _arb_vec_init(nn), cf = _arb_vec_init(nn), cfl = _arb_vec_init(nn);
    for (int i = 0; i < nn; i++) {
        if (!tau_exact(tau + i, ns[i], bprec)) { arb_log_ui(tau + i, ns[i], bprec); arb_div(tau + i, tau + i, L, bprec); }
        arb_log_ui(cf + i, ps[i], bprec); arb_sqrt_ui(tmp, ns[i], bprec); arb_div(cf + i, cf + i, tmp, bprec);
        arb_set_round(cfl + i, cf + i, accprec);
    }
    /* breakpoints 0, 1, |1 - t_n| (all in [0,1] since t_n < 2), sorted by midpoint */
    int nb = 0; arb_ptr bps = _arb_vec_init(nn + 2);
    arb_zero(bps + nb++); arb_one(bps + nb++);
    for (int i = 0; i < nn; i++) {
        arb_sub(tmp, one, tau + i, bprec); arb_abs(bps + nb, tmp);
        arb_get_lbound_arf(arb_midref(tmp2), bps + nb, 64);  /* filter as bc.py: lower >= -1e-30, upper <= 1 + 1e-30 */
        double lo = arf_get_d(arb_midref(tmp2), ARF_RND_DOWN);
        arb_get_ubound_arf(arb_midref(tmp2), bps + nb, 64);
        double hi = arf_get_d(arb_midref(tmp2), ARF_RND_UP);
        if (lo >= -1e-30 && hi <= 1 + 1e-30) nb++;
    }
    for (int i = 1; i < nb; i++) for (int j = i; j > 0 && arf_cmp(arb_midref(bps + j), arb_midref(bps + j - 1)) < 0; j--) arb_swap(bps + j, bps + j - 1);
    /* piece selection: GRAMS_PIECES (explicit list) or ip % K == k */
    int sel[300]; int nsel = 0; char plist[2048]; plist[0] = 0;
    if (nb - 1 > 300) die("too many pieces");
    for (int ip = 0; ip + 1 < nb; ip++) sel[ip] = (ip % pK == pk);
    if (getenv("GRAMS_PIECES")) {
        for (int ip = 0; ip + 1 < nb; ip++) sel[ip] = 0;
        char buf[2048]; strncpy(buf, getenv("GRAMS_PIECES"), sizeof buf - 1); buf[sizeof buf - 1] = 0;
        for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
            int ip = atoi(tok); if (ip < 0 || ip + 1 >= nb) die("GRAMS_PIECES: piece index out of range");
            if (sel[ip]) die("GRAMS_PIECES: duplicate piece"); sel[ip] = 1;
        }
    }
    for (int ip = 0; ip + 1 < nb; ip++) if (sel[ip]) {
        char t[16]; snprintf(t, sizeof t, "%s%d", nsel ? ", " : "", ip); strncat(plist, t, sizeof plist - strlen(plist) - 1); nsel++; }
    if (nsel == 0) die("no pieces selected");
    /* Gauss-Legendre rules, computed on demand (exact rule: N+1 nodes) */
    slong ng = N + 1, cache_n[64]; arb_ptr cache_x[64], cache_w[64]; int ncache = 0;
    printf("grams_c: mu=%s N=%ld prec=%ld bprec=%ld pieces=%d gauss(%ld) shifts n=", tag, (long)N, (long)prec, (long)bprec, nb - 1, (long)ng);
    for (int i = 0; i < nn; i++) printf("%ld%s", ns[i], i + 1 < nn ? "," : "");
    printf("  part %ld/%ld  %.0fs\n", pk, pK, now() - t0); fflush(stdout);
    slong dims[2] = {(N + 1)/2, N/2};
    int act[2] = {1, 1};                 /* GRAMS_PARITY=p: accumulate parity p only (half the products, half the file) */
    if (getenv("GRAMS_PARITY")) { int gp = atoi(getenv("GRAMS_PARITY")); if (gp != 0 && gp != 1) die("GRAMS_PARITY must be 0 or 1"); act[1 - gp] = 0; }
    arb_mat_t A0[2], B0[2];
    for (int p = 0; p < 2; p++) { slong dd = act[p] ? dims[p] : 0; arb_mat_init(A0[p], dd, dd); arb_mat_init(B0[p], dd, dd); }
    const slong chunk = 256;
    arb_t md, hh, u, w, y, p0, p1, t; arb_init(md); arb_init(hh); arb_init(u); arb_init(w); arb_init(y); arb_init(p0); arb_init(p1); arb_init(t);
    arb_ptr cv = _arb_vec_init(N), at = _arb_vec_init(2*nn), ac = _arb_vec_init(2*nn), ifac = _arb_vec_init(N);
    arb_one(ifac + 0);
    for (slong j = 1; j < N; j++) arb_div_ui(ifac + j, ifac + j - 1, j, accprec);     /* 1/j! */
    double nsteps = 0;
    long maxchunks = getenv("GRAMS_MAXCHUNKS") ? atol(getenv("GRAMS_MAXCHUNKS")) : -1, nchunks = 0;   /* profiling only */
    for (int ip = 0; ip + 1 < nb; ip++) {
        if (!sel[ip]) continue;
        arb_add(md, bps + ip, bps + ip + 1, bprec); arb_mul_2exp_si(md, md, -1);
        arb_sub(hh, bps + ip + 1, bps + ip, bprec); arb_mul_2exp_si(hh, hh, -1);
        if (arb_is_zero(hh)) { printf("  piece %d [%.6f,%.6f] has exact length 0: skipped (contribution exactly 0)\n", ip, arf_get_d(arb_midref(bps+ip), ARF_RND_NEAR), arf_get_d(arb_midref(bps+ip+1), ARF_RND_NEAR)); continue; }
        if (!arb_is_positive(hh)) die("piece of undetermined length (breakpoints not separated)");
        int na = 0;
        for (int i = 0; i < nn; i++) {
            arb_add(tmp, md, tau + i, bprec);
            if (arb_lt(tmp, one)) { arb_set(at + na, tau + i); arb_set(ac + na, cfl + i); na++; }
            else if (!arb_gt(tmp, one)) die("ambiguous + shift");
            arb_sub(tmp, md, tau + i, bprec);
            if (arb_gt(tmp, mone)) { arb_neg(at + na, tau + i); arb_set(ac + na, cfl + i); na++; }
            else if (!arb_lt(tmp, mone)) die("ambiguous - shift");
        }
        mag_t eA, eB; mag_init(eA); mag_init(eB); double rused = 0;
        slong ngp = bound_mode ? plan_piece(eA, eB, &rused, md, hh, at, ac, na, N, prec) : ng;
        if (getenv("GRAMS_PLANONLY")) {     /* planning only: node count per piece, no quadrature */
            printf("  plan piece %d [%.6f,%.6f] active shifts %d  gauss %ld (rho %.2f, err A %.2e, B %.2e)\n", ip, arf_get_d(arb_midref(bps+ip), ARF_RND_NEAR), arf_get_d(arb_midref(bps+ip+1), ARF_RND_NEAR), na, (long)ngp, rused, mag_get_d(eA), mag_get_d(eB)); fflush(stdout);
            mag_clear(eA); mag_clear(eB); continue;
        }
        int ci = -1;
        for (int z = 0; z < ncache; z++) if (cache_n[z] == ngp) ci = z;
        if (ci < 0) {
            if (ncache >= 64) die("gauss cache full");
            ci = ncache++; cache_n[ci] = ngp; cache_x[ci] = _arb_vec_init(ngp); cache_w[ci] = _arb_vec_init(ngp);
            for (slong k = 0; k < ngp; k++) arb_hypgeom_legendre_p_ui_root(cache_x[ci] + k, cache_w[ci] + k, ngp, k, bprec);
        }
        arb_ptr gx = cache_x[ci], gw = cache_w[ci];
        printf("  piece %d [%.6f,%.6f] active shifts %d  gauss %ld%s (rho %.2f, err A %.2e, B %.2e)  %.0fs\n", ip, arf_get_d(arb_midref(bps+ip), ARF_RND_NEAR), arf_get_d(arb_midref(bps+ip+1), ARF_RND_NEAR), na, (long)ngp, ngp == ng ? " exact" : "", rused, mag_get_d(eA), mag_get_d(eB), now() - t0); fflush(stdout);
        for (slong c0 = 0; c0 < ngp; c0 += chunk) {
            slong m = FLINT_MIN(chunk, ngp - c0);
            arb_mat_t Pm[2], Cm[2], WT[2], T;
            for (int p = 0; p < 2; p++) { slong dd = act[p] ? dims[p] : 0, mm = act[p] ? m : 0; arb_mat_init(Pm[p], dd, mm); arb_mat_init(Cm[p], dd, mm); arb_mat_init(WT[p], mm, dd); }
            double tl0 = now();
            for (slong q = 0; q < m; q++) {
                arb_mul(u, hh, gx + c0 + q, bprec); arb_add(u, u, md, bprec);
                arb_mul(w, hh, gw + c0 + q, bprec);
                /* P_j(u) */
                arb_one(p0); arb_set(p1, u);
                for (slong j = 0; j < N; j++) {
                    arb_srcptr v = (j == 0) ? p0 : p1;
                    if (j >= 2) {
                        ulong jj = j - 1;               /* R_j from R_{j-1} = p1, R_{j-2} = p0 */
                        arb_mul(t, u, p1, bprec); arb_mul_ui(t, t, 2*jj + 1, bprec);
                        arb_mul_ui(tmp, p0, jj*jj, bprec); arb_sub(t, t, tmp, bprec);
                        arb_swap(p0, p1); arb_swap(p1, t); v = p1;
                    }
                    if (!act[j % 2]) continue;
                    arb_set_round(tmp, v, accprec); arb_mul(tmp, tmp, ifac + j, accprec);
                    arb_set_round(arb_mat_entry(Pm[j % 2], j/2, q), tmp, prec);
                    if (!arb_is_finite(arb_mat_entry(Pm[j % 2], j/2, q))) die("non-finite Legendre value");
                }
                /* (C_p P_j)(u) = sum over active shifts c * P_j(u + t) */
                for (slong j = 0; j < N; j++) arb_zero(cv + j);
                for (int a = 0; a < na; a++) {
                    arb_add(y, u, at + a, bprec);
                    arb_one(p0); arb_set(p1, y);
                    for (slong j = 0; j < N; j++) {
                        arb_srcptr v = (j == 0) ? p0 : p1;
                        if (j >= 2) {
                            ulong jj = j - 1;
                            arb_mul(t, y, p1, bprec); arb_mul_ui(t, t, 2*jj + 1, bprec);
                            arb_mul_ui(tmp2, p0, jj*jj, bprec); arb_sub(t, t, tmp2, bprec);
                            arb_swap(p0, p1); arb_swap(p1, t); v = p1;
                        }
                        arb_set_round(tmp, v, accprec);
                        arb_addmul(cv + j, ac + a, tmp, accprec);
                    }
                }
                nsteps += (double) N*(1 + na);
                for (slong j = 0; j < N; j++) {
                    if (!act[j % 2]) continue;
                    arb_ptr cj = arb_mat_entry(Cm[j % 2], j/2, q);
                    arb_mul(cv + j, cv + j, ifac + j, accprec);          /* sum_a c_a R_j(y_a) / j! */
                    arb_set_round(cj, cv + j, prec);
                    if (!arb_is_finite(cj)) die("non-finite shifted Legendre sum");
                    arb_mul(arb_mat_entry(WT[j % 2], q, j/2), cj, w, prec);
                }
            }
            double tl1 = now();
            for (int p = 0; p < 2; p++) {
                if (act[p]) {
                    arb_mat_init(T, dims[p], dims[p]);
                    arb_mat_mul(T, Pm[p], WT[p], prec); arb_mat_add(A0[p], A0[p], T, prec);
                    arb_mat_mul(T, Cm[p], WT[p], prec); arb_mat_add(B0[p], B0[p], T, prec);
                    arb_mat_clear(T);
                }
                arb_mat_clear(Pm[p]); arb_mat_clear(Cm[p]); arb_mat_clear(WT[p]);
            }
            printf("    nodes %ld/%ld  %.0fs  (legendre steps %.3e; chunk: values %.1fs, products %.1fs)\n", (long)(c0 + m), (long)ngp, now() - t0, nsteps, tl1 - tl0, now() - tl1); fflush(stdout);
            if (maxchunks >= 0 && ++nchunks >= maxchunks) { printf("PROFILE STOP (GRAMS_MAXCHUNKS): output not written\n"); return 3; }
        }
        if (!mag_is_zero(eA) || !mag_is_zero(eB))           /* quadrature error of this piece, every raw entry */
            for (int p = 0; p < 2; p++) if (act[p]) for (slong i = 0; i < dims[p]; i++) for (slong j = 0; j < dims[p]; j++) {
                arb_add_error_mag(arb_mat_entry(A0[p], i, j), eA); arb_add_error_mag(arb_mat_entry(B0[p], i, j), eB); }
        mag_clear(eA); mag_clear(eB);
    }
    if (getenv("GRAMS_PLANONLY")) { printf("PLAN ONLY: output not written  %.0fs\n", now() - t0); return 4; }
    /* write */
    FILE *f = fopen(out, "wb"); if (!f) die("cannot open output");
    char head[4096], mh[512] = "";
    for (int k = 0; k < 4; k++) { int pp = k % 2; if (!act[pp]) continue;      /* key order A0raw_0, A0raw_1, B0raw_0, B0raw_1 */
        char t[96]; snprintf(t, sizeof t, "\"%s_%d\": [\"M\", %ld, %ld], ", k < 2 ? "A0raw" : "B0raw", pp, (long)dims[pp], (long)dims[pp]); strncat(mh, t, sizeof mh - strlen(mh) - 1); }
    snprintf(head, sizeof head, "{%s\"N\": [\"o\", %ld], \"prec\": [\"o\", %ld], \"part\": [\"o\", [%ld, %ld]], \"tag\": [\"o\", \"%s\"], \"npieces\": [\"o\", %d], \"pieces\": [\"o\", [%s]], \"mode\": [\"o\", \"%s\"], \"parities\": [\"o\", [%s]]}",
             mh, (long)N, (long)prec, pk, pK, tag, nb - 1, plist, bound_mode ? "bound" : "exact", act[0] && act[1] ? "0, 1" : (act[0] ? "0" : "1"));
    unsigned long long hl = strlen(head); fwrite(&hl, 8, 1, f); fwrite(head, 1, hl, f);
    arb_mat_struct *mats[4] = {A0[0], A0[1], B0[0], B0[1]};
    for (int k = 0; k < 4; k++)
        if (act[k % 2]) for (slong i = 0; i < arb_mat_nrows(mats[k]); i++)
            for (slong j = 0; j < arb_mat_ncols(mats[k]); j++) write_arb(f, arb_mat_entry(mats[k], i, j));
    if (fflush(f) != 0 || ferror(f)) { fclose(f); remove(out); die("write error (disk full?): output removed"); }   /* 2026-10-21: a full disk silently truncated a part */
    if (fclose(f) != 0) { remove(out); die("write failed: output removed"); }
    double ra = 0, rb = 0;
    for (int p = 0; p < 2; p++) if (act[p]) for (slong i = 0; i < dims[p]; i++) for (slong j = 0; j < dims[p]; j++) {
        double x = mag_get_d(arb_radref(arb_mat_entry(A0[p], i, j))); if (x > ra) ra = x;
        x = mag_get_d(arb_radref(arb_mat_entry(B0[p], i, j))); if (x > rb) rb = x; }
    printf("grams_c done: %s  max raw radius A0 %.2e B0 %.2e  %.0fs  flint %s\n", out, ra, rb, now() - t0, FLINT_VERSION);
    return 0;
}
