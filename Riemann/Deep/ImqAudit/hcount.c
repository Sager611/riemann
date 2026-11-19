/* hcount.c (ImqAudit): exact h(D), D = -d < -4 fundamental, by counting reduced forms (a, b, c): b^2 - 4ac = D,
 * |b| <= a <= c, b >= 0 if |b| = a or a = c, gcd(a, b, c) = 1 (unconditional; every class has exactly one).
 * b runs over b = d mod 2, 0 <= b <= sqrt(d/3); N_b = (b^2 + d)/4 is factored by sieving the roots of x^2 = -d mod p
 * (Euler criterion, Tonelli-Shanks), and the divisors a of N_b with b <= a <= N_b/a are counted.
 * Also M = #(reduced forms with a <= X = floor(sqrt(d/4))) and M2 = sum_{a <= X} r(a) by a multiplicative
 * depth-first enumeration (Legendre symbols by Euler's criterion; written independently of isieve2.c).
 * Usage: hcount < file_of_d    -> lines "d fund h M M2"; with MONLY=1 in the environment only M2 ("d fund - - M2")     Build: cc -O2 -o hcount hcount.c -lm */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef unsigned long long u64; typedef __uint128_t u128;
static u64 isqrt(u64 x) { u64 r = (u64)sqrtl((long double)x); while (r && r * r > x) r--; while ((r + 1) * (r + 1) <= x) r++; return r; }
static u64 gcd(u64 a, u64 b) { while (b) { u64 t = a % b; a = b; b = t; } return a; }
static u64 mulmod(u64 a, u64 b, u64 m) { return (u64)((u128)a * b % m); }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = mulmod(r, a, m); a = mulmod(a, a, m); e >>= 1; } return r; }
static int legendre(u64 a, u64 p) { a %= p; if (!a) return 0; return powmod(a, (p - 1) / 2, p) == 1 ? 1 : -1; }
static u64 sqrtmod(u64 n, u64 p) { /* p odd prime, n a nonzero QR */
  n %= p; if (p % 4 == 3) return powmod(n, (p + 1) / 4, p);
  u64 q = p - 1; int s = 0; while (!(q & 1)) { q >>= 1; s++; }
  u64 z = 2; while (legendre(z, p) != -1) z++;
  int m = s; u64 c = powmod(z, q, p), t = powmod(n, q, p), r = powmod(n, (q + 1) / 2, p);
  while (t != 1) { int i = 0; u64 tt = t; while (tt != 1) { tt = mulmod(tt, tt, p); i++; }
    u64 b = c; for (int j = 0; j < m - i - 1; j++) b = mulmod(b, b, p); m = i; c = mulmod(b, b, p); t = mulmod(t, c, p); r = mulmod(r, b, p); }
  if (mulmod(r, r, p) != n) { fprintf(stderr, "sqrtmod fail\n"); exit(3); }
  return r; }
static u64 *pr; static int npr;
static int fundamental(u64 d) { int r = (int)(d & 15); u64 u;
  if (r == 3 || r == 7 || r == 11 || r == 15) u = d; else if (r == 4 || r == 8) u = d / 4; else return 0;
  for (int i = 0; i < npr && pr[i] * pr[i] <= u; i++) if (u % (pr[i] * pr[i]) == 0) return 0;
  return 1; }
#define MF 16
static u64 *rem, *fp; static unsigned char *fe, *nf;
static u64 X_, B_; static u64 cntH, cntM; static u64 curb, curN;
static u64 dvp[MF]; static int dve[MF], dvn;
static void divs(int k, u64 a) {
  if (k == dvn) { if (a < curb || a > curN / a) return; u64 c = curN / a;
    if (gcd(gcd(a, curb), c) != 1) return;
    u64 w = (curb == 0) ? 1 : 1 + ((curb < a && a < c) ? 1 : 0);
    cntH += w; if (a <= X_) cntM += w; return; }
  u64 x = a; for (int e = 0; e <= dve[k]; e++) { divs(k + 1, x); if (e < dve[k]) { if (x > curN / dvp[k]) break; x *= dvp[k]; } } }
static signed char *chi; static int nchi;
static u64 dfs(u64 n, int idx, u64 w) { u64 tot = w;
  for (int i = idx; i < nchi && n * pr[i] <= X_; i++) { int s = chi[i]; if (s < 0) continue;
    if (s == 0) { tot += dfs(n * pr[i], i + 1, w); continue; }
    u64 m = n * pr[i]; for (;;) { tot += dfs(m, i + 1, 2 * w); if (m > X_ / pr[i]) break; m *= pr[i]; } }
  return tot; }
int main(void) {
  u64 PMAX = 1100000; char *comp = calloc(PMAX + 1, 1); pr = malloc(sizeof(u64) * 100000); npr = 0;
  for (u64 i = 2; i <= PMAX; i++) { if (comp[i]) continue; pr[npr++] = i; for (u64 j = i * i; j <= PMAX; j += i) comp[j] = 1; }
  free(comp); chi = malloc(npr);
  size_t cap = 600000; rem = malloc(8 * cap); fp = malloc(8 * cap * MF); fe = malloc(cap * MF); nf = malloc(cap);
  u64 d; int monly = getenv("MONLY") != 0;
  while (scanf("%llu", &d) == 1) {
    if (d < 5 || d > 1000000000000ULL) { printf("%llu skip\n", d); continue; }
    int fu = fundamental(d); if (!fu) { printf("%llu 0 - - -\n", d); fflush(stdout); continue; }
    X_ = isqrt(d / 4); B_ = isqrt(d / 3); u64 b0 = d & 1;
    if (monly) goto mpart;
    size_t nb = (B_ >= b0) ? (size_t)((B_ - b0) / 2 + 1) : 0; if (nb > cap) { fprintf(stderr, "cap\n"); return 2; }
    u64 maxN = 0;
    for (size_t i = 0; i < nb; i++) { u64 b = b0 + 2 * i; if ((b * b + d) % 4) { fprintf(stderr, "parity\n"); return 2; }
      rem[i] = (b * b + d) / 4; nf[i] = 0; if (rem[i] > maxN) maxN = rem[i]; }
    u64 P = isqrt(maxN) + 1; if (P > PMAX) { fprintf(stderr, "PMAX\n"); return 2; }
    for (size_t i = 0; i < nb; i++) { int e = 0; while (!(rem[i] & 1)) { rem[i] >>= 1; e++; } if (e) { fp[i * MF] = 2; fe[i * MF] = (unsigned char)e; nf[i] = 1; } }
    for (int k = 1; k < npr && pr[k] <= P; k++) { u64 p = pr[k]; u64 n = (p - d % p) % p; u64 roots[2]; int nr = 0;
      if (n == 0) { roots[0] = 0; nr = 1; } else if (legendre(n, p) == 1) { u64 r = sqrtmod(n, p); roots[0] = r; roots[1] = p - r; nr = 2; } else continue;
      u64 inv2 = (p + 1) / 2;
      for (int t = 0; t < nr; t++) { u64 i0 = mulmod((roots[t] + p - b0 % p) % p, inv2, p);
        for (u64 i = i0; i < nb; i += p) { int e = 0; while (rem[i] % p == 0) { rem[i] /= p; e++; }
          if (!e) { fprintf(stderr, "sieve miss d=%llu p=%llu\n", d, p); return 2; }
          if (nf[i] >= MF) { fprintf(stderr, "MF\n"); return 2; } fp[i * MF + nf[i]] = p; fe[i * MF + nf[i]] = (unsigned char)e; nf[i]++; } } }
    cntH = cntM = 0;
    for (size_t i = 0; i < nb; i++) { curb = b0 + 2 * i; curN = (curb * curb + d) / 4; dvn = nf[i];
      for (int k = 0; k < dvn; k++) { dvp[k] = fp[i * MF + k]; dve[k] = fe[i * MF + k]; }
      if (rem[i] > 1) { if (dvn >= MF) { fprintf(stderr, "MF\n"); return 2; } dvp[dvn] = rem[i]; dve[dvn] = 1; dvn++; }
      { u64 chk = 1; for (int k = 0; k < dvn; k++) for (int e = 0; e < dve[k]; e++) chk *= dvp[k]; if (chk != curN) { fprintf(stderr, "factor mismatch\n"); return 2; } }
      divs(0, 1); }
  mpart:
    nchi = 0; while (nchi < npr && pr[nchi] <= X_) { u64 p = pr[nchi]; int s;
      if (p == 2) s = (d & 1) ? (((d & 7) == 7) ? 1 : -1) : 0; else s = legendre((p - d % p) % p, p);
      chi[nchi++] = (signed char)s; }
    u64 M2 = dfs(1, 0, 1);
    if (monly) printf("%llu 1 - - %llu\n", d, M2); else printf("%llu 1 %llu %llu %llu\n", d, cntH, cntM, M2);
  }
  return 0; }
