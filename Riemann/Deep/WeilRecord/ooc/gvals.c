/* gvals.c -- Legendre values at the Gram quadrature nodes, for the out-of-core Weil certificate (weilooc, 2026-11-05).
   Replaces the value part of grams_c.c (same pieces, shifts, coefficients and bound-mode node plan, code copied), but
   (1) evaluates P_j by the "rotation" recurrence (ooc/DERIVATION.md section 1) in fixed point with FV fractional bits:
         a_{j+1} = j/(j+1) (c a_j + s b_j),  b_{j+1} = -s a_j + c b_j,  P_{j+1} = a_{j+1} + c P_j,
         a_0 = 1, b_0 = 0, P_0 = 1, c = u, s = sqrt(1 - u^2);  (a_j, b_j) = D_j = P_j - e^{i theta} P_{j-1}.
       The exact transfer map is diag(j/(j+1), 1) R(theta) with norm exactly 1, so rounding errors do not grow:
       |P_j(computed) - P_j(u)| <= 2.1 (j+1)^2 2^-FV (no ball arithmetic, no 1.28 N extra bits);
   (2) writes the scaled values  Phat[q][i] = sqrt(w_q (2j+1)) P_j(u_q),  Chat[q][i] = sqrt(w_q (2j+1)) (C_p P_j)(u_q)
       (j = 2i + parity) as fixed-point digits (F fractional bits, T total, offset binary, uint16 little endian) to
       OUTDIR/gv_p<parity>_pc<piece>.bin, shape (nodes of the piece, d_parity, T/16), and the error data to
       OUTDIR/gv_meta_pc<piece>.json.  A0 = Phat^T Chat and B0 = Chat^T Chat are then formed by fxmm (exact products).
   Build: cc -O2 -o gvals gvals.c -I/opt/homebrew/include -L/opt/homebrew/lib -lflint -lgmp -lmpfr
   Usage: gvals TAG N F T FV OUTDIR [piece list "i,j,.." or "all"]                                                 */
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
#include <flint/acb.h>

static double t0;
static double now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9*ts.tv_nsec; }
static void die(const char *m) { fprintf(stderr, "FATAL: %s\n", m); fflush(stdout); exit(2); }

/* ---------- copied from grams_c.c: bound-mode plan (Bernstein-Walsh + ATAP Thm 19.3) ---------- */
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
static long vm_int(long n) { for (long p = 2; p <= n; p++) if (n % p == 0) { long m = n; while (m % p == 0) m /= p; return m == 1 ? p : 0; } return 0; }
static int kind; static long MU, LP, LQ;
static void L_at(arb_t L, slong prec) {
    if (kind == 0) { arb_log_ui(L, MU, prec); arb_mul_2exp_si(L, L, -1); }
    else if (kind == 1) { arb_set_ui(L, 3); arb_mul_2exp_si(L, L, -1); }
    else { fmpq_t q; fmpq_init(q); fmpq_set_si(q, LP, LQ); arb_set_fmpq(L, q, prec); fmpq_clear(q); }
}
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

/* ---------- fixed point helpers (mpz, FV fractional bits) ---------- */
static slong FV;
static mpz_t half;                        /* 2^(FV-1) */
static void rshift_round(mpz_t r, const mpz_t x, slong k) {   /* r = round_half_up(x / 2^k) */
    mpz_t h; mpz_init_set_ui(h, 1); mpz_mul_2exp(h, h, k - 1); mpz_add(r, x, h); mpz_fdiv_q_2exp(r, r, k); mpz_clear(h);
}
/* exact dyadic midpoint of a ball, rounded to FV bits; returns |result - x| bound added into err (mag) */
static void arb_to_fix(mpz_t r, const arb_t x, mag_t err) {
    fmpz_t m, e; fmpz_init(m); fmpz_init(e);
    arf_t t; arf_init(t);
    arf_mul_2exp_si(t, arb_midref(x), FV);
    arf_get_fmpz(m, t, ARF_RND_NEAR);        /* |m - mid 2^FV| <= 1/2 */
    fmpz_get_mpz(r, m);
    mag_t h; mag_init(h); mag_set_ui_2exp_si(h, 1, -FV - 1);
    mag_add(err, arb_radref(x), h);
    mag_clear(h); arf_clear(t); fmpz_clear(m); fmpz_clear(e);
}

/* Legendre values P_0..P_{N-1} at u (fixed point, |u| <= 1), rotation recurrence; s = sqrt(1-u^2) fixed point */
static void legendre_fix(mpz_t *out, slong N, const mpz_t u, const mpz_t s, mpz_t *tmp) {
    mpz_ptr a = tmp[0], b = tmp[1], t1 = tmp[2], t2 = tmp[3], na = tmp[4], nb = tmp[5];
    mpz_set_ui(a, 1); mpz_mul_2exp(a, a, FV); mpz_set_ui(b, 0);
    mpz_set(out[0], a);
    for (slong j = 0; j + 1 < N; j++) {
        /* na = round( round((c a + s b)/2^FV) * j / (j+1) ) */
        mpz_mul(t1, u, a); mpz_addmul(t1, s, b);
        mpz_add(t1, t1, half); mpz_fdiv_q_2exp(t1, t1, FV);
        mpz_mul_ui(t1, t1, (unsigned long) j);
        /* round(t1/(j+1)) = floor((2 t1 + (j+1)) / (2 (j+1))) */
        mpz_mul_2exp(t1, t1, 1); mpz_add_ui(t1, t1, (unsigned long)(j + 1)); mpz_fdiv_q_ui(na, t1, (unsigned long)(2*(j + 1)));
        /* nb = round((c b - s a)/2^FV) */
        mpz_mul(t2, u, b); mpz_submul(t2, s, a);
        mpz_add(t2, t2, half); mpz_fdiv_q_2exp(nb, t2, FV);
        /* P_{j+1} = na + round(c P_j / 2^FV) */
        mpz_mul(t1, u, out[j]); mpz_add(t1, t1, half); mpz_fdiv_q_2exp(t1, t1, FV);
        mpz_add(out[j + 1], na, t1);
        mpz_swap(a, na); mpz_swap(b, nb);
    }
}

static slong F_, T_, ND;
static mpz_t offT;    /* 2^(T-1) */
static void put_digits(unsigned short *dst, const mpz_t X) {   /* offset binary, ND uint16 little endian */
    mpz_t Y; mpz_init(Y); mpz_add(Y, X, offT);
    if (mpz_sgn(Y) < 0 || mpz_sizeinbase(Y, 2) > (size_t) T_) die("value out of fixed-point range");
    memset(dst, 0, ND*sizeof(unsigned short));
    size_t cnt = 0; mpz_export(dst, &cnt, -1, 2, -1, 0, Y);    /* least significant word first, little endian words */
    if (cnt > (size_t) ND) die("export overflow");
    mpz_clear(Y);
}

int main(int argc, char **argv) {
    if (argc < 7) die("usage: gvals TAG N F T FV OUTDIR [pieces]");
    const char *tag = argv[1]; slong N = atol(argv[2]); F_ = atol(argv[3]); T_ = atol(argv[4]); FV = atol(argv[5]);
    const char *outdir = argv[6]; const char *plist = argc > 7 ? argv[7] : "all";
    if (F_ % 16 || T_ % 16 || FV < F_ + 32 || N < 4) die("bad F/T/FV/N");
    ND = T_/16;
    if (strcmp(tag, "e3") == 0) kind = 1;
    else if (tag[0] == 'L') { kind = 2; if (sscanf(tag + 1, "%ldo%ld", &LP, &LQ) != 2 || LP <= 0 || LQ <= 0) die("bad L tag"); }
    else { kind = 0; MU = atol(tag); if (MU < 2) die("bad mu"); }
    slong prec = FV + 128;              /* ball precision for setup quantities */
    t0 = now();
    mpz_init(half); mpz_set_ui(half, 1); mpz_mul_2exp(half, half, FV - 1);
    mpz_init(offT); mpz_set_ui(offT, 1); mpz_mul_2exp(offT, offT, T_ - 1);
    arb_t L, one, mone, tmp, tmp2; arb_init(L); arb_init(one); arb_init(mone); arb_init(tmp); arb_init(tmp2);
    L_at(L, prec); arb_one(one); arb_set_si(mone, -1);
    long ns[256], ps[256]; int nn = 0;
    { arb_t twoL, ln; arb_init(twoL); arb_init(ln); arb_mul_2exp_si(twoL, L, 1);
      long top = (kind == 0) ? MU : (long) ceil(exp(2.0*(kind == 1 ? 1.5 : (double)LP/LQ))) + 1;
      for (long n = 2; n <= top; n++) {
          long p = vm_int(n); if (!p) continue;
          int ok = 0;
          if (kind == 0) ok = n < MU;
          else { arb_log_ui(ln, n, prec); if (arb_lt(ln, twoL)) ok = 1; else if (arb_gt(ln, twoL)) ok = 0; else die("cannot decide n < mu"); }
          if (ok) { ns[nn] = n; ps[nn] = p; nn++; }
      }
      arb_clear(twoL); arb_clear(ln); }
    arb_ptr tau = _arb_vec_init(nn), cf = _arb_vec_init(nn);
    for (int i = 0; i < nn; i++) {
        if (!tau_exact(tau + i, ns[i], prec)) { arb_log_ui(tau + i, ns[i], prec); arb_div(tau + i, tau + i, L, prec); }
        arb_log_ui(cf + i, ps[i], prec); arb_sqrt_ui(tmp, ns[i], prec); arb_div(cf + i, cf + i, tmp, prec);
    }
    int nb = 0; arb_ptr bps = _arb_vec_init(nn + 2);
    arb_zero(bps + nb++); arb_one(bps + nb++);
    for (int i = 0; i < nn; i++) {
        arb_sub(tmp, one, tau + i, prec); arb_abs(bps + nb, tmp);
        arb_get_lbound_arf(arb_midref(tmp2), bps + nb, 64);
        double lo = arf_get_d(arb_midref(tmp2), ARF_RND_DOWN);
        arb_get_ubound_arf(arb_midref(tmp2), bps + nb, 64);
        double hi = arf_get_d(arb_midref(tmp2), ARF_RND_UP);
        if (lo >= -1e-30 && hi <= 1 + 1e-30) nb++;
    }
    for (int i = 1; i < nb; i++) for (int j = i; j > 0 && arf_cmp(arb_midref(bps + j), arb_midref(bps + j - 1)) < 0; j--) arb_swap(bps + j, bps + j - 1);
    int sel[300]; if (nb - 1 > 300) die("too many pieces");
    for (int ip = 0; ip + 1 < nb; ip++) sel[ip] = strcmp(plist, "all") == 0;
    if (strcmp(plist, "all") != 0) {
        char buf[2048]; strncpy(buf, plist, sizeof buf - 1); buf[sizeof buf - 1] = 0;
        for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) { int ip = atoi(tok); if (ip < 0 || ip + 1 >= nb) die("piece out of range"); sel[ip] = 1; }
    }
    printf("gvals: tag=%s N=%ld F=%ld T=%ld FV=%ld pieces=%d shifts n=", tag, (long)N, (long)F_, (long)T_, (long)FV, nb - 1);
    for (int i = 0; i < nn; i++) printf("%ld%s", ns[i], i + 1 < nn ? "," : "");
    printf("\n"); fflush(stdout);
    /* fixed-point coefficients and sqrt(2j+1) */
    mpz_t *cfx = malloc(sizeof(mpz_t)*nn); mag_t cferr; mag_init(cferr);     /* max |cfx - cf| */
    for (int i = 0; i < nn; i++) { mag_t e; mag_init(e); mpz_init(cfx[i]); arb_to_fix(cfx[i], cf + i, e); mag_max(cferr, cferr, e); mag_clear(e); }
    arb_ptr sq = _arb_vec_init(N);
    for (slong j = 0; j < N; j++) arb_sqrt_ui(sq + j, 2*j + 1, prec);
    slong dims[2] = {(N + 1)/2, N/2};
    mpz_t *Pu = malloc(sizeof(mpz_t)*N), *Pt = malloc(sizeof(mpz_t)*N), *Cs = malloc(sizeof(mpz_t)*N), tmpz[8];
    for (slong j = 0; j < N; j++) { mpz_init(Pu[j]); mpz_init(Pt[j]); mpz_init(Cs[j]); }
    for (int k = 0; k < 8; k++) mpz_init(tmpz[k]);
    mpz_t ufix, sfix, yfix, X, sig; mpz_init(ufix); mpz_init(sfix); mpz_init(yfix); mpz_init(X); mpz_init(sig);
    mpz_t taufix[512]; mag_t tauerr; mag_init(tauerr);
    for (int i = 0; i < nn; i++) { mag_t e; mag_init(e); mpz_init(taufix[i]); arb_to_fix(taufix[i], tau + i, e); mag_max(tauerr, tauerr, e); mag_clear(e); }
    unsigned short *rowP[2], *rowC[2];
    for (int p = 0; p < 2; p++) { rowP[p] = malloc(sizeof(unsigned short)*ND*dims[p]); rowC[p] = malloc(sizeof(unsigned short)*ND*dims[p]); }
    arb_t md, hh, u, w, x, sw, sv, csum; arb_init(md); arb_init(hh); arb_init(u); arb_init(w); arb_init(x); arb_init(sw); arb_init(sv); arb_init(csum);
    mpz_t one_fix; mpz_init(one_fix); mpz_set_ui(one_fix, 1); mpz_mul_2exp(one_fix, one_fix, FV);
    for (int ip = 0; ip + 1 < nb; ip++) {
        if (!sel[ip]) continue;
        double tp = now();
        arb_add(md, bps + ip, bps + ip + 1, prec); arb_mul_2exp_si(md, md, -1);
        arb_sub(hh, bps + ip + 1, bps + ip, prec); arb_mul_2exp_si(hh, hh, -1);
        char path[1024];
        snprintf(path, sizeof path, "%s/gv_meta_pc%d.json", outdir, ip);
        if (arb_is_zero(hh)) {
            FILE *fm = fopen(path, "w"); if (!fm) die("cannot write meta");
            fprintf(fm, "{\"piece\": %d, \"nodes\": 0, \"empty\": true}\n", ip); fclose(fm);
            printf("  piece %d: exact length 0, skipped\n", ip); continue;
        }
        if (!arb_is_positive(hh)) die("piece of undetermined length");
        /* active shifts (as grams_c): +t if md + t < 1, -t if md - t > -1 */
        int na = 0; arb_ptr at = _arb_vec_init(2*nn), ac = _arb_vec_init(2*nn); int aidx[512]; int asgn[512];
        for (int i = 0; i < nn; i++) {
            arb_add(tmp, md, tau + i, prec);
            if (arb_lt(tmp, one)) { arb_set(at + na, tau + i); arb_set(ac + na, cf + i); aidx[na] = i; asgn[na] = 1; na++; }
            else if (!arb_gt(tmp, one)) die("ambiguous + shift");
            arb_sub(tmp, md, tau + i, prec);
            if (arb_gt(tmp, mone)) { arb_neg(at + na, tau + i); arb_set(ac + na, cf + i); aidx[na] = i; asgn[na] = -1; na++; }
            else if (!arb_lt(tmp, mone)) die("ambiguous - shift");
        }
        mag_t eA, eB; mag_init(eA); mag_init(eB); double rused = 0;
        slong ng = plan_piece(eA, eB, &rused, md, hh, at, ac, na, N, F_);
        /* c_sum = sum of active |c_a| (upper bound) */
        arb_zero(csum); for (int a = 0; a < na; a++) arb_add(csum, csum, ac + a, prec);
        /* Gauss-Legendre rule at prec */
        arb_ptr gx = _arb_vec_init(ng), gw = _arb_vec_init(ng);
        for (slong k = 0; k < ng; k++) arb_hypgeom_legendre_p_ui_root(gx + k, gw + k, ng, k, prec);
        FILE *fo[2][2];
        for (int p = 0; p < 2; p++) for (int z = 0; z < 2; z++) {
            snprintf(path, sizeof path, "%s/gv_p%d_pc%d_%s.bin", outdir, p, ip, z ? "C" : "P");
            fo[p][z] = fopen(path, "wb"); if (!fo[p][z]) die("cannot open output");
        }
        mag_t du, dy, dsg, dsc, maxw; mag_init(du); mag_init(dy); mag_init(dsg); mag_init(dsc); mag_init(maxw);
        arb_t sumw, sumsw; arb_init(sumw); arb_init(sumsw); arb_zero(sumw); arb_zero(sumsw);
        for (slong q = 0; q < ng; q++) {
            arb_mul(u, hh, gx + q, prec); arb_add(u, u, md, prec);
            arb_mul(w, hh, gw + q, prec); arb_add(sumw, sumw, w, prec);
            mag_t e; mag_init(e);
            arb_to_fix(ufix, u, e); mag_max(du, du, e);
            if (mpz_cmpabs(ufix, one_fix) > 0) die("node outside [-1,1]");
            /* s = sqrt(1 - u^2) at the fixed-point node (exact dyadic) */
            { fmpz_t zz; fmpz_init(zz); fmpz_set_mpz(zz, ufix); arb_set_fmpz(x, zz); arb_mul_2exp_si(x, x, -FV); fmpz_clear(zz); }
            arb_sqr(tmp, x, prec); arb_sub_ui(tmp, tmp, 1, prec); arb_neg(tmp, tmp); arb_sqrtpos(tmp, tmp, prec);
            mag_zero(e); arb_to_fix(sfix, tmp, e); mag_max(dsc, dsc, e);
            legendre_fix(Pu, N, ufix, sfix, tmpz);
            for (slong j = 0; j < N; j++) mpz_set_ui(Cs[j], 0);
            for (int a = 0; a < na; a++) {
                /* y = u +- tau, rounded to FV bits (exact integer add of fixed-point values) */
                if (asgn[a] > 0) mpz_add(yfix, ufix, taufix[aidx[a]]); else mpz_sub(yfix, ufix, taufix[aidx[a]]);
                if (mpz_cmpabs(yfix, one_fix) > 0) die("shifted node outside [-1,1]");
                { fmpz_t zz; fmpz_init(zz); fmpz_set_mpz(zz, yfix); arb_set_fmpz(x, zz); arb_mul_2exp_si(x, x, -FV); fmpz_clear(zz); }
                arb_sqr(tmp, x, prec); arb_sub_ui(tmp, tmp, 1, prec); arb_neg(tmp, tmp); arb_sqrtpos(tmp, tmp, prec);
                mag_zero(e); arb_to_fix(sfix, tmp, e); mag_max(dsc, dsc, e);
                legendre_fix(Pt, N, yfix, sfix, tmpz);
                for (slong j = 0; j < N; j++) mpz_addmul(Cs[j], cfx[aidx[a]], Pt[j]);       /* scale 2^-2FV */
            }
            for (slong j = 0; j < N; j++) { mpz_add(Cs[j], Cs[j], half); mpz_fdiv_q_2exp(Cs[j], Cs[j], FV); }  /* round to FV */
            /* scaled outputs: round_F( sigma_j * value ), sigma_j = sqrt(w (2j+1)) */
            arb_sqrtpos(sw, w, prec); arb_add(sumsw, sumsw, sw, prec);
            mag_t mw; mag_init(mw); arb_get_mag(mw, w); mag_max(maxw, maxw, mw); mag_clear(mw);
            for (slong j = 0; j < N; j++) {
                int p = j % 2; slong i = j/2;
                arb_mul(sv, sw, sq + j, prec);
                mag_zero(e); arb_to_fix(sig, sv, e); mag_max(dsg, dsg, e);
                mpz_mul(X, sig, Pu[j]); rshift_round(X, X, 2*FV - F_); put_digits(rowP[p] + ND*i, X);
                mpz_mul(X, sig, Cs[j]); rshift_round(X, X, 2*FV - F_); put_digits(rowC[p] + ND*i, X);
            }
            for (int p = 0; p < 2; p++) {
                if (fwrite(rowP[p], sizeof(unsigned short), ND*dims[p], fo[p][0]) != (size_t)(ND*dims[p])) die("write error");
                if (fwrite(rowC[p], sizeof(unsigned short), ND*dims[p], fo[p][1]) != (size_t)(ND*dims[p])) die("write error");
            }
            mag_clear(e);
        }
        for (int p = 0; p < 2; p++) for (int z = 0; z < 2; z++) { if (fflush(fo[p][z]) || ferror(fo[p][z])) die("write error (disk full?)"); fclose(fo[p][z]); }
        /* shifted-node displacement: |y_fix - y_true| <= du + tauerr + 0 (fixed-point add is exact) */
        mag_add(dy, du, tauerr);
        snprintf(path, sizeof path, "%s/gv_meta_pc%d.json", outdir, ip);
        FILE *fm = fopen(path, "w"); if (!fm) die("cannot write meta");
        fprintf(fm, "{\"piece\": %d, \"nodes\": %ld, \"N\": %ld, \"F\": %ld, \"T\": %ld, \"FV\": %ld, \"active_shifts\": %d, \"rho\": %.3f, ",
                ip, (long)ng, (long)N, (long)F_, (long)T_, (long)FV, na, rused);
        fprintf(fm, "\"eA\": \"%.17e\", \"eB\": \"%.17e\", \"du\": \"%.17e\", \"dy\": \"%.17e\", \"ds\": \"%.17e\", \"dsig\": \"%.17e\", \"cferr\": \"%.17e\", ",
                mag_get_d(eA)*(1 + 1e-15), mag_get_d(eB)*(1 + 1e-15), mag_get_d(du)*(1 + 1e-15), mag_get_d(dy)*(1 + 1e-15), mag_get_d(dsc)*(1 + 1e-15), mag_get_d(dsg)*(1 + 1e-15), mag_get_d(cferr)*(1 + 1e-15));
        arf_t ub; arf_init(ub);
        arb_get_ubound_arf(ub, csum, 64); double csu = arf_get_d(ub, ARF_RND_UP);
        arb_get_ubound_arf(ub, sumw, 64); double swu = arf_get_d(ub, ARF_RND_UP);
        arb_get_ubound_arf(ub, sumsw, 64); double sswu = arf_get_d(ub, ARF_RND_UP);
        arf_clear(ub);
        fprintf(fm, "\"csum\": \"%.17e\", \"sumw\": \"%.17e\", \"sumsqrtw\": \"%.17e\", \"maxw\": \"%.17e\", \"lo\": %.9f, \"hi\": %.9f, \"time\": %.1f}\n",
                csu*(1 + 1e-15), swu*(1 + 1e-15), sswu*(1 + 1e-15), mag_get_d(maxw)*(1 + 1e-15),
                arf_get_d(arb_midref(bps + ip), ARF_RND_NEAR), arf_get_d(arb_midref(bps + ip + 1), ARF_RND_NEAR), now() - tp);
        if (fflush(fm) || ferror(fm)) die("meta write error"); fclose(fm);
        printf("  piece %d [%.6f,%.6f] shifts %d nodes %ld (rho %.2f, eA %.2e, eB %.2e)  %.1fs (total %.0fs)\n", ip,
               arf_get_d(arb_midref(bps + ip), ARF_RND_NEAR), arf_get_d(arb_midref(bps + ip + 1), ARF_RND_NEAR), na, (long)ng, rused,
               mag_get_d(eA), mag_get_d(eB), now() - tp, now() - t0);
        fflush(stdout);
        _arb_vec_clear(gx, ng); _arb_vec_clear(gw, ng); _arb_vec_clear(at, 2*nn); _arb_vec_clear(ac, 2*nn);
        mag_clear(eA); mag_clear(eB); mag_clear(du); mag_clear(dy); mag_clear(dsg); mag_clear(dsc); mag_clear(maxw); arb_clear(sumw); arb_clear(sumsw);
    }
    printf("gvals done %.0fs\n", now() - t0);
    return 0;
}
