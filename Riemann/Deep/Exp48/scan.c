/* exp48 scan.c: unconditional scan of ALL negative discriminants D (fundamental or not), 3 <= |D| <= X,
 * for form class groups Cl(D) of exponent dividing 8. No L-function input, no GRH: exact form arithmetic.
 *   1. For split primes q (q prime, q !| D, (D/q) = 1) in increasing order:
 *        if 4 q^8 < |D|: the class of the prime form (q,b,c) has order > 8 (EKN Lemma 1; holds for orders,
 *        q !| conductor), so exp(Cl(D)) does not divide 8 -> reject;
 *        else compute f_q^8 by three squarings; if not principal -> reject (explicit element of order !| 8).
 *   2. Survivors (all tested split primes have f^8 = 1): enumerate every reduced primitive form, h = their
 *      number (exact), test f^2, f^4, f^8 for each; E | 8 iff all f^8 = 1; structure from the counts.
 * Output (stdout): D fund h E n2 n4  for every D with E | 8.   Stats on stderr.
 * Build: cc -O2 -o scan scan.c -lm      Run: ./scan X [Xlo]
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef __int128 i128;
typedef int64_t i64;
typedef struct { i64 a, b, c; } form;
static i64 DD; /* current discriminant (negative) */

static i64 gcdl(i64 x, i64 y) { if (x < 0) x = -x; if (y < 0) y = -y; while (y) { i64 t = x % y; x = y; y = t; } return x; }
/* extended gcd: returns g = gcd(x,y) >= 0 and u,v with u x + v y = g */
static i64 egcd(i64 x, i64 y, i64 *u, i64 *v) {
  i64 u0 = 1, v0 = 0, u1 = 0, v1 = 1;
  while (y) { i64 q = x / y, t;
    t = x - q * y; x = y; y = t;
    t = u0 - q * u1; u0 = u1; u1 = t;
    t = v0 - q * v1; v0 = v1; v1 = t; }
  if (x < 0) { x = -x; u0 = -u0; v0 = -v0; }
  *u = u0; *v = v0; return x;
}
static i64 cfromab(i64 a, i64 b) { i128 t = (i128)b * b - DD; return (i64)(t / (4 * (i128)a)); }
static void normalize(form *f) { /* b into (-a, a] */
  i64 a = f->a, b = f->b, m = 2 * a;
  i64 r = b % m; if (r < 0) r += m; if (r > a) r -= m;
  f->b = r; f->c = cfromab(a, r);
}
static void reduce(form *f) {
  normalize(f);
  while (f->a > f->c || (f->a == f->c && f->b < 0)) {
    i64 t = f->a; f->a = f->c; f->c = t; f->b = -f->b; normalize(f);
  }
}
/* Cohen, A Course in Comp. Alg. NT, Alg. 5.4.7 (composition of positive definite forms) */
static form compose(form f1, form f2) {
  if (f1.a > f2.a) { form t = f1; f1 = f2; f2 = t; }
  i64 a1 = f1.a, b1 = f1.b, a2 = f2.a, b2 = f2.b, c2 = f2.c;
  i64 s = (b1 + b2) / 2, n = b2 - s;
  i64 u, v, d, y1, d1, x2, y2;
  if (a2 % a1 == 0) { y1 = 0; d = a1; }
  else { d = egcd(a2, a1, &u, &v); y1 = u; }
  if (s % d == 0) { y2 = -1; x2 = 0; d1 = d; }
  else { d1 = egcd(s, d, &x2, &y2); y2 = -y2; }
  i64 v1 = a1 / d1, v2 = a2 / d1;
  i128 r = ((i128)y1 * y2 % v1 * (n % v1) - (i128)x2 * (c2 % v1)) % v1; if (r < 0) r += v1;
  form g; g.a = v1 * v2; g.b = (i64)(b2 + 2 * (i128)v2 * r);
  g.c = cfromab(g.a, g.b);
  reduce(&g); return g;
}
static int isprinc(form f) { return f.a == 1; }
static form pow8(form f, int *o2, int *o4) { /* returns f^8; flags f^2 = 1, f^4 = 1 */
  form g = compose(f, f); *o2 = isprinc(g);
  g = compose(g, g); *o4 = isprinc(g);
  g = compose(g, g); return g;
}
#define NP 200
static int PR[NP];
static signed char *SQ[NP];   /* SQ[i][r] = square root of r mod p (or -1) */
static int kron(i64 D, int p, int *root) { /* Kronecker (D/p), D<0 ; root mod p if split, odd p */
  if (p == 2) { i64 m = ((D % 8) + 8) % 8; if (m % 2 == 0) return 0; return (m == 1 || m == 7) ? 1 : -1; }
  i64 r = D % p; if (r < 0) r += p;
  if (r == 0) return 0;
  return 0; /* unused */
}
static i64 issqfree_odd_part(i64 m) { /* returns 1 if m squarefree */
  for (i64 p = 2; p * p <= m; p++) if (m % (p * p) == 0) return 0;
  return 1;
}
static int isfund(i64 D) { i64 n = -D;
  if (n % 4 == 3) return issqfree_odd_part(n);            /* D = 1 mod 4 */
  if (n % 4 == 0) { i64 m = n / 4; if (m % 4 == 1 || m % 4 == 2) return issqfree_odd_part(m); } /* D/4 = 2,3 mod 4 */
  return 0;
}
int main(int argc, char **argv) {
  i64 X = atoll(argv[1]); i64 Xlo = argc > 2 ? atoll(argv[2]) : 3;
  /* primes */
  int cnt = 0; for (int p = 2; cnt < NP; p++) { int ok = 1; for (int q = 2; q * q <= p; q++) if (p % q == 0) { ok = 0; break; } if (ok) PR[cnt++] = p; }
  int *sqt[NP];
  for (int i = 1; i < NP; i++) { int p = PR[i]; sqt[i] = malloc(sizeof(int) * p); for (int r = 0; r < p; r++) sqt[i][r] = -1;
    for (int t = 0; t < p; t++) { int r = (int)((i64)t * t % p); if (sqt[i][r] < 0) sqt[i][r] = t; } }
  i64 nrejL1 = 0, nrejpow = 0, nfull = 0, nfullfail = 0, nout = 0;
  for (i64 n = Xlo; n <= X; n++) {
    if (n % 4 == 1 || n % 4 == 2) continue;   /* D = -n must be 0,1 mod 4 */
    i64 D = -n; DD = D;
    int tested = 0, rejected = 0;
    for (int i = 0; i < NP && tested < 6; i++) {
      int p = PR[i]; form f; int k;
      if (p == 2) { i64 m = n % 8; /* D mod 8 = (8 - m) mod 8 */ i64 Dm = (8 - m) % 8; if (Dm % 2 == 0) continue; if (Dm != 1) continue; f.a = 2; f.b = 1; }
      else { i64 r = (D % p + p) % p; if (r == 0) continue; int t = sqt[i][r]; if (t < 0) continue;
             i64 b = t; if (((b - D) & 1) != 0) b += p; f.a = p; f.b = b; }
      /* split prime p */
      i128 p8 = (i128)p * p; p8 = p8 * p8; p8 = p8 * p8;
      if (4 * p8 < (i128)n) { rejected = 1; nrejL1++; break; }
      f.c = cfromab(f.a, f.b); reduce(&f);
      int o2, o4; form g = pow8(f, &o2, &o4);
      if (!isprinc(g)) { rejected = 1; nrejpow++; break; }
      tested++;
      (void)k;
    }
    if (rejected) continue;
    /* full check */
    nfull++;
    i64 h = 0, n2 = 0, n4 = 0, n8 = 0; int bad = 0;
    i64 amax = (i64)sqrt((double)n / 3.0) + 2;
    for (i64 a = 1; a <= amax && !bad; a++) {
      for (i64 b = (n & 1); b <= a; b += 2) {
        i128 t = (i128)b * b + n; if (t % (4 * a) != 0) continue;
        i64 c = (i64)(t / (4 * a)); if (c < a) continue;
        if (gcdl(gcdl(a, b), c) != 1) continue;
        int mult = (b == 0 || b == a || c == a) ? 1 : 2;
        form f = { a, b, c }; int o2, o4; form g = pow8(f, &o2, &o4);
        h += mult;
        if (!isprinc(g)) { bad = 1; break; }
        n8 += mult; if (o4) n4 += mult; if (o2) n2 += mult;
      }
    }
    if (bad) { nfullfail++; continue; }
    int E = (h == 1) ? 1 : (n2 == h ? 2 : (n4 == h ? 4 : 8));
    printf("%lld %d %lld %d %lld %lld\n", (long long)D, isfund(D), (long long)h, E, (long long)n2, (long long)n4);
    nout++;
  }
  fprintf(stderr, "X=%lld Xlo=%lld rejL1=%lld rejpow=%lld full=%lld fullfail=%lld out=%lld\n", (long long)X, (long long)Xlo, (long long)nrejL1, (long long)nrejpow, (long long)nfull, (long long)nfullfail, (long long)nout);
  return 0;
}
