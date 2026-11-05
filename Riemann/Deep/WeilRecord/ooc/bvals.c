/* bvals.c -- spherical Bessel values at the band nodes, for the out-of-core Weil certificate (weilooc, 2026-11-05).
   Same values as cert2.sbessel (start values j_{N+1}, j_N from arb's J_{n+1/2} at raised precision, downward
   recurrence j_{n-1} = (2n+1)/x j_n - j_{n+1} in ball arithmetic at BPREC = F + x_max + 128 bits), at x = L t_q.
   Output, per parity p and node q (node-major, fixed point F/T, offset-binary uint16 digits):
       Bhat[q][i] = sqrt((2L/pi) (2j+1) w_q |g_q|) * sgn_j * j_j(L t_q),  j = 2i + p,  sgn_j = (-1)^(floor(j/2)),
   with w_q, g_q the MIDPOINTS of the weight and of a(t_q) - gamma (their radii are reported for the error budget),
   so K = sum_{g>0} Bhat^T Bhat - sum_{g<0} Bhat^T Bhat.  Also the first column of the band Gram,
       b0[i] = (2L/pi) sqrt((2j+1)(2p+1)) sum_q w_q s_j(t_q) s_p(t_q)   (ball arithmetic, weights as balls),
   written as "mid_man mid_exp rad" lines (exact dyadic midpoint, decimal radius upper bound).
   Usage: bvals TAG N F T BPREC NODEFILE OUTDIR   (TAG: L<p>o<q> or integer mu)                                       */
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
#include <flint/arb_hypgeom.h>

static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }
static void die(const char *m) { fprintf(stderr, "FATAL: %s\n", m); fflush(stdout); exit(2); }
static int kind; static long MU, LP, LQ;
static void L_at(arb_t L, slong prec) {
    if (kind == 0) { arb_log_ui(L, MU, prec); arb_mul_2exp_si(L, L, -1); }
    else { fmpq_t q; fmpq_init(q); fmpq_set_si(q, LP, LQ); arb_set_fmpq(L, q, prec); fmpq_clear(q); }
}
static slong F_, T_, ND;
static mpz_t offT;
static void put_digits(unsigned short *dst, const fmpz_t Xf) {
    mpz_t X, Y; mpz_init(X); mpz_init(Y); fmpz_get_mpz(X, Xf); mpz_add(Y, X, offT);
    if (mpz_sgn(Y) < 0 || mpz_sizeinbase(Y, 2) > (size_t) T_) die("value out of fixed-point range");
    memset(dst, 0, ND*sizeof(unsigned short));
    size_t cnt = 0; mpz_export(dst, &cnt, -1, 2, -1, 0, Y);
    if (cnt > (size_t) ND) die("export overflow");
    mpz_clear(X); mpz_clear(Y);
}
/* j_n(x) = sqrt(pi/(2x)) J_{n+1/2}(x) with relative accuracy >= bprec - 40 bits (precision raised as needed) */
static void bessel_start(arb_t v, slong n, const arb_t t, slong bprec) {
    slong p = bprec;
    arb_t x, nu, L, c; arb_init(x); arb_init(nu); arb_init(L); arb_init(c);
    while (1) {
        L_at(L, p); arb_mul(x, L, t, p);
        arb_set_ui(nu, 2*n + 1); arb_mul_2exp_si(nu, nu, -1);
        arb_hypgeom_bessel_j(v, nu, x, p);
        arb_const_pi(c, p); arb_div(c, c, x, p); arb_mul_2exp_si(c, c, -1); arb_sqrt(c, c, p); arb_mul(v, v, c, p);
        if (!arb_is_zero(v) && arb_rel_accuracy_bits(v) >= bprec - 40) break;
        if (p > 60000) die("bessel start value not accurate");
        p = (slong)(p*1.5);
    }
    arb_set_round(v, v, bprec);
    arb_clear(x); arb_clear(nu); arb_clear(L); arb_clear(c);
}
int main(int argc, char **argv) {
    if (argc < 8) die("usage: bvals TAG N F T BPREC NODEFILE OUTDIR");
    const char *tag = argv[1]; slong N = atol(argv[2]); F_ = atol(argv[3]); T_ = atol(argv[4]); slong bprec = atol(argv[5]);
    const char *nf = argv[6], *outdir = argv[7];
    ND = T_/16; if (F_ % 16 || T_ % 16) die("F, T must be multiples of 16");
    if (tag[0] == 'L') { kind = 2; if (sscanf(tag + 1, "%ldo%ld", &LP, &LQ) != 2) die("bad tag"); } else { kind = 0; MU = atol(tag); }
    mpz_init(offT); mpz_set_ui(offT, 1); mpz_mul_2exp(offT, offT, T_ - 1);
    slong prec = F_ + 64;                     /* precision of scale factors and products (balls) */
    double t0 = now();
    FILE *fn = fopen(nf, "r"); if (!fn) die("cannot open node file");
    /* read nodes */
    slong cap = 1024, nq = 0; arb_ptr ts = _arb_vec_init(cap), ws = _arb_vec_init(cap), gs = _arb_vec_init(cap);
    char s1[4096], s2[256], s3[4096], s4[256], s5[64], s6[4096], s7[256], s8[64];
    while (fscanf(fn, "%4095s %255s %4095s %255s %63s %4095s %255s %63s", s1, s2, s3, s4, s5, s6, s7, s8) == 8) {
        if (nq == cap) {
            arb_ptr a = _arb_vec_init(2*cap), b = _arb_vec_init(2*cap), c = _arb_vec_init(2*cap);
            for (slong i = 0; i < nq; i++) { arb_swap(a + i, ts + i); arb_swap(b + i, ws + i); arb_swap(c + i, gs + i); }
            _arb_vec_clear(ts, cap); _arb_vec_clear(ws, cap); _arb_vec_clear(gs, cap); ts = a; ws = b; gs = c; cap *= 2;
        }
        fmpz_t m, e; fmpz_init(m); fmpz_init(e);
        fmpz_set_str(m, s1, 10); fmpz_set_str(e, s2, 10); arb_set_fmpz_2exp(ts + nq, m, e);
        fmpz_set_str(m, s3, 10); fmpz_set_str(e, s4, 10); arb_set_fmpz_2exp(ws + nq, m, e);
        { arb_t r; arb_init(r); arb_set_d(r, strtod(s5, NULL)); arb_get_mag(arb_radref(ws + nq), r); arb_clear(r); }   /* radius as given (upper bound) */
        fmpz_set_str(m, s6, 10); fmpz_set_str(e, s7, 10); arb_set_fmpz_2exp(gs + nq, m, e);
        { arb_t r; arb_init(r); arb_set_d(r, strtod(s8, NULL)); arb_get_mag(arb_radref(gs + nq), r); arb_clear(r); }
        fmpz_clear(m); fmpz_clear(e);
        nq++;
    }
    fclose(fn);
    printf("bvals: tag=%s N=%ld F=%ld T=%ld bprec=%ld nodes=%ld\n", tag, (long)N, (long)F_, (long)T_, (long)bprec, (long)nq); fflush(stdout);
    slong dims[2] = {(N + 1)/2, N/2};
    arb_t L, Lh, c2Lpi, x, ix, tmp, sv, sig, wm, gm; arb_init(L); arb_init(Lh); arb_init(c2Lpi); arb_init(x); arb_init(ix); arb_init(tmp); arb_init(sv); arb_init(sig); arb_init(wm); arb_init(gm);
    L_at(L, prec); L_at(Lh, bprec);
    arb_const_pi(tmp, prec); arb_mul_2exp_si(c2Lpi, L, 1); arb_div(c2Lpi, c2Lpi, tmp, prec);       /* 2L/pi */
    arb_ptr js = _arb_vec_init(N + 2);
    arb_ptr col[2]; col[0] = _arb_vec_init(dims[0]); col[1] = _arb_vec_init(dims[1]);
    unsigned short *row[2]; row[0] = malloc(sizeof(unsigned short)*ND*dims[0]); row[1] = malloc(sizeof(unsigned short)*ND*dims[1]);
    FILE *fo[2]; char path[1024];
    for (int p = 0; p < 2; p++) { snprintf(path, sizeof path, "%s/bv_p%d.bin", outdir, p); fo[p] = fopen(path, "wb"); if (!fo[p]) die("cannot open output"); }
    mag_t maxerr[2], maxrad_s, wrad, grad; mag_init(maxerr[0]); mag_init(maxerr[1]); mag_init(maxrad_s); mag_init(wrad); mag_init(grad);
    mag_t half; mag_init(half); mag_set_ui_2exp_si(half, 1, -F_ - 1);
    fmpz_t Xi; fmpz_init(Xi); arf_t tf; arf_init(tf);
    arb_t radcap; arb_init(radcap);
    for (slong q = 0; q < nq; q++) {
        /* x = L t at bprec (t exact) */
        arb_mul(x, Lh, ts + q, bprec);
        bessel_start(js + N + 1, N + 1, ts + q, bprec);
        bessel_start(js + N, N, ts + q, bprec);
        arb_inv(ix, x, bprec);
        for (slong n = N; n >= 1; n--) {
            arb_mul_ui(tmp, ix, 2*n + 1, bprec); arb_mul(tmp, tmp, js + n, bprec); arb_sub(js + n - 1, tmp, js + n + 1, bprec);
        }
        for (slong n = 0; n < N; n++) {
            if (!arb_is_finite(js + n)) die("non-finite Bessel value");
            mag_max(maxrad_s, maxrad_s, arb_radref(js + n));
        }
        mag_max(wrad, wrad, arb_radref(ws + q)); mag_max(grad, grad, arb_radref(gs + q));
        /* sigma_q = sqrt((2L/pi) w_mid |g_mid|) (ball);  scaled values */
        arb_set_arf(wm, arb_midref(ws + q)); arb_set_arf(gm, arb_midref(gs + q)); arb_abs(gm, gm);
        arb_mul(sig, wm, gm, prec); arb_mul(sig, sig, c2Lpi, prec); arb_sqrtpos(sig, sig, prec);
        for (slong j = 0; j < N; j++) {
            int p = j % 2; slong i = j/2;
            arb_sqrt_ui(tmp, 2*j + 1, prec); arb_mul(sv, sig, tmp, prec);
            arb_set_round(tmp, js + j, prec);
            if ((j/2) % 2) arb_neg(tmp, tmp);
            arb_mul(sv, sv, tmp, prec);
            arf_mul_2exp_si(tf, arb_midref(sv), F_); arf_get_fmpz(Xi, tf, ARF_RND_NEAR);
            put_digits(row[p] + ND*i, Xi);
            mag_t e; mag_init(e); mag_add(e, arb_radref(sv), half); mag_max(maxerr[p], maxerr[p], e); mag_clear(e);
            /* b0 column: w_q s_j s_p (balls; weight as ball) */
            arb_set_round(tmp, js + j, prec); if ((j/2) % 2) arb_neg(tmp, tmp);
            arb_t s0; arb_init(s0); arb_set_round(s0, js + p, prec); if ((p/2) % 2) arb_neg(s0, s0);
            arb_mul(s0, s0, tmp, prec); arb_addmul(col[p] + i, s0, ws + q, prec); arb_clear(s0);
        }
        for (int p = 0; p < 2; p++) if (fwrite(row[p], sizeof(unsigned short), ND*dims[p], fo[p]) != (size_t)(ND*dims[p])) die("write error");
        if (q % 500 == 499) { printf("  nodes %ld/%ld  %.0fs\n", (long)q + 1, (long)nq, now() - t0); fflush(stdout); }
    }
    for (int p = 0; p < 2; p++) { if (fflush(fo[p]) || ferror(fo[p])) die("write error (disk full?)"); fclose(fo[p]); }
    for (int p = 0; p < 2; p++) {
        snprintf(path, sizeof path, "%s/bv_b0_p%d.txt", outdir, p);
        FILE *fb = fopen(path, "w"); if (!fb) die("cannot write b0");
        for (slong i = 0; i < dims[p]; i++) {
            slong j = 2*i + p;
            arb_mul_ui(tmp, c2Lpi, 1, prec); arb_sqrt_ui(sv, (2*j + 1)*(2*p + 1), prec); arb_mul(tmp, tmp, sv, prec);
            arb_mul(tmp, tmp, col[p] + i, prec);
            fmpz_t m, e; fmpz_init(m); fmpz_init(e); arf_get_fmpz_2exp(m, e, arb_midref(tmp));
            char *ms = fmpz_get_str(NULL, 10, m), *es = fmpz_get_str(NULL, 10, e);
            fprintf(fb, "%s %s %.17e\n", ms, es, mag_get_d(arb_radref(tmp))*(1 + 1e-14) + 1e-300);
            flint_free(ms); flint_free(es); fmpz_clear(m); fmpz_clear(e);
        }
        if (fflush(fb) || ferror(fb)) die("b0 write error"); fclose(fb);
    }
    snprintf(path, sizeof path, "%s/bv_meta.json", outdir);
    FILE *fm = fopen(path, "w"); if (!fm) die("cannot write meta");
    fprintf(fm, "{\"N\": %ld, \"F\": %ld, \"T\": %ld, \"bprec\": %ld, \"nodes\": %ld, \"maxerr0\": \"%.17e\", \"maxerr1\": \"%.17e\", \"maxrad_bessel\": \"%.17e\", \"wrad\": \"%.17e\", \"grad\": \"%.17e\", \"time\": %.1f}\n",
            (long)N, (long)F_, (long)T_, (long)bprec, (long)nq, mag_get_d(maxerr[0])*(1 + 1e-14), mag_get_d(maxerr[1])*(1 + 1e-14), mag_get_d(maxrad_s)*(1 + 1e-14),
            mag_get_d(wrad)*(1 + 1e-14), mag_get_d(grad)*(1 + 1e-14), now() - t0);
    if (fflush(fm) || ferror(fm)) die("meta write error"); fclose(fm);
    printf("bvals done %.0fs  maxerr %.2e %.2e  max Bessel radius %.2e\n", now() - t0, mag_get_d(maxerr[0]), mag_get_d(maxerr[1]), mag_get_d(maxrad_s));
    return 0;
}
