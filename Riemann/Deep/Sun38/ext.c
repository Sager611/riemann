/* sun38 ext.c: Sun 1309.1679 TeX Conj3.10 (printed 3.8) for PRIME q in [lo, hi], every a0 in F_q.
   For each (q, a0) find a "P2" certificate: odd t and d with gcd(d, n) = 2 (n = q-1) such that
   a0 + g^t and a0 + g^(t+d) are primitive (g = least primitive root). Then the circular permutation
   x_{2m} = g^(m d), x_{2m+1} = g^(t - m d) (m = 0..n/2-1) is a Hamilton cycle of G_{q,a0}.
   The search uses exp/log tables; every certificate found is then re-checked independently:
   g^t by modular exponentiation and primitivity of z by z^(n/r) != 1 for every prime r | n.
   Output per q: "q n_cert n_fail [a0 ...]" (fail = no P2 certificate among the first 64 odd t in T).
   Usage: ext LO HI > out.txt          (single thread, memory O(HI)) */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

static int isprime(long m) {
    if (m < 2) return 0;
    for (long d = 2; d * d <= m; d++) if (m % d == 0) return 0;
    return 1;
}
static int factor(long m, long *ps) {
    int k = 0;
    for (long d = 2; d * d <= m; d++)
        if (m % d == 0) { ps[k++] = d; while (m % d == 0) m /= d; }
    if (m > 1) ps[k++] = m;
    return k;
}
static uint64_t mpow(uint64_t b, uint64_t e, uint64_t m) {
    uint64_t r = 1 % m; b %= m;
    while (e) { if (e & 1) r = r * b % m; b = b * b % m; e >>= 1; }
    return r;
}
static long gcdl(long a, long b) { if (a < 0) a = -a; while (b) { long t = a % b; a = b; b = t; } return a; }
static int isprim(uint64_t z, long q, const long *ps, int np) {
    if (z % q == 0) return 0;
    for (int i = 0; i < np; i++) if (mpow(z, (q - 1) / ps[i], q) == 1) return 0;
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: ext LO HI\n"); return 2; }
    long lo = atol(argv[1]), hi = atol(argv[2]);
    int *expt = malloc(sizeof(int) * (hi + 1)), *logt = malloc(sizeof(int) * (hi + 1));
    unsigned char *unit = malloc(hi + 1);
    long ps[32]; int odds[80];
    long long tot = 0, totfail = 0, totbad = 0;
    for (long q = lo < 3 ? 3 : lo; q <= hi; q++) {
        if (!isprime(q)) continue;
        long n = q - 1;
        int np = factor(n, ps);
        long g = 2;
        while (!isprim(g, q, ps, np)) g++;
        long e = 1;
        for (long i = 0; i < n; i++) { expt[i] = (int)e; logt[e] = (int)i; e = e * g % q; }
        for (long s = 0; s < n; s++) {
            int u = 1;
            for (int i = 0; i < np; i++) if (s % ps[i] == 0) { u = 0; break; }
            unit[s] = (unsigned char)u;
        }
        long ncert = 0, nfail = 0;
        long fails[64]; int nf = 0;
        for (long a0 = 0; a0 < q; a0++) {
            int no = 0, found = 0; long T = -1, D = -1;
            for (long s = 1; s < n && !found; s += 2) {
                long z = a0 + expt[s]; if (z >= q) z -= q;
                if (z == 0 || !unit[logt[z]]) continue;
                for (int i = 0; i < no; i++)
                    if (gcdl(s - odds[i], n) == 2) { T = odds[i]; D = s - odds[i]; found = 1; break; }
                if (!found) { odds[no++] = (int)s; if (no >= 64) break; }
            }
            if (!found) { nfail++; if (nf < 64) fails[nf++] = a0; continue; }
            /* independent re-check */
            int ok = (T % 2 == 1) && (D > 0) && (D < n) && gcdl(D, n) == 2;
            uint64_t z1 = (a0 + mpow(g, T, q)) % q, z2 = (a0 + mpow(g, (T + D) % n, q)) % q;
            ok = ok && isprim(z1, q, ps, np) && isprim(z2, q, ps, np);
            if (!ok) { totbad++; fprintf(stderr, "BAD q=%ld a0=%ld t=%ld d=%ld\n", q, a0, T, D); }
            else ncert++;
        }
        tot += q; totfail += nfail;
        printf("%ld %ld %ld", q, ncert, nfail);
        for (int i = 0; i < nf; i++) printf(" %ld", fails[i]);
        printf("\n");
        fflush(stdout);
    }
    fprintf(stderr, "range %ld..%ld: pairs %lld, no-certificate %lld, bad %lld\n", lo, hi, tot, totfail, totbad);
    return totbad ? 1 : 0;
}
