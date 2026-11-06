/* sievecheck.c: independent C re-implementation of the split-prime sieve (Corollary B) over the FULL range of cert.py:
   every n in [1001, nmax] with Delta = n^2 - 4N < 10^18 (both N = +1, -1); prints n surviving "(Delta/p) != 1 for all primes p with 4p^4 < Delta".
   Quadratic-residue tables by squaring (no Euler criterion), 128-bit comparisons.  cc -O2 -o sievecheck sievecheck.c && ./sievecheck */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
typedef unsigned __int128 u128;
#define PMAX 23000
int main(void) {
    static char comp[PMAX + 1];
    int primes[4000], np = 0;
    for (int i = 2; i <= PMAX; i++) { if (!comp[i]) { primes[np++] = i; for (int j = 2 * i; j <= PMAX; j += i) comp[j] = 1; } }
    const uint64_t DS = 1000000000000000000ULL;   /* 10^18 */
    for (int N = 1; N >= -1; N -= 2) {
        /* split tables: sp[k][r] = 1 iff ((r^2 - 4N) mod p) is a nonzero square mod p, p = primes[k] odd */
        char **sp = malloc(sizeof(char *) * np);
        for (int k = 1; k < np; k++) {
            int p = primes[k];
            char *qr = calloc(p, 1);
            for (long x = 1; x < p; x++) qr[(x * x) % p] = 1;
            sp[k] = calloc(p, 1);
            for (long r = 0; r < p; r++) { long v = ((r * r - 4L * N) % p + p) % p; sp[k][r] = (v != 0 && qr[v]); }
            free(qr);
        }
        uint64_t n = 1001, count = 0, surv = 0;
        for (;; n++) {
            u128 D = (u128)n * n - (u128)(4 * (long)N);   /* N = -1 adds 4 */
            if (N == -1) D = (u128)n * n + 4;
            if (D >= (u128)DS && (u128)(n - 1) * (n - 1) + 4 >= (u128)DS + 8) break;
            if (D >= (u128)DS + 8) break;   /* covers every Delta_n < 10^18 (superset by a few n) */
            count++;
            int ok = 1;
            if ((uint64_t)(D % 8) == 1) ok = 0;                 /* p = 2 splits iff Delta = 1 mod 8 */
            for (int k = 1; ok && k < np; k++) {
                u128 p = primes[k];
                if (4 * p * p * p * p >= D) break;
                if (sp[k][n % primes[k]]) ok = 0;
            }
            if (ok) { surv++; printf("N=%+d n=%llu\n", N, (unsigned long long)n); }
        }
        fprintf(stderr, "N=%+d: tested %llu values of n in [1001, %llu], survivors %llu\n", N, (unsigned long long)count, (unsigned long long)(n - 1), (unsigned long long)surv);
        for (int k = 1; k < np; k++) free(sp[k]);
        free(sp);
    }
    return 0;
}
