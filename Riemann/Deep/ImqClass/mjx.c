/* mjx.c: extract (d, h) with h <= NMAX from one 2^28-slab of the Mosunov-Jacobson class group tables
 * (LMFDB files cl{3,7}mod8.K.gz, cl{4,8}mod16.K.gz), with checks, and update a running state:
 * per-h counts and largest d, and SHA-256 of the sorted list "d h\n" for each threshold N in THR.
 * Checks per file: (i) every d read is the next fundamental discriminant of its residue class in the
 * slab (so the file lists exactly the fundamental d of its class: completeness and positions);
 * (ii) the invariants multiply to h and form a divisibility chain; (iii) gzip stream complete (CRC).
 * Usage: mjx STATE K DIR      (state created at K = 0)
 *        mjx --report STATE OUTPREFIX   (writes OUTPREFIX_perh.txt and OUTPREFIX_sha.txt)
 * Build: cc -O2 -o mjx mjx.c -lz
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <zlib.h>
#include <CommonCrypto/CommonDigest.h>

typedef unsigned long long u64;
#define NMAX 10000
#define SLAB (1ULL << 28)
static const int THR[] = {1000, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000};
#define NTHR ((int)(sizeof(THR) / sizeof(THR[0])))
static const int CR[4] = {3, 7, 4, 8}, CMOD[4] = {8, 8, 16, 16};
static const char *CNAME[4] = {"cl3mod8", "cl7mod8", "cl4mod16", "cl8mod16"};

typedef struct {
  u64 magic, nmax, next_k;
  u64 cnt[NMAX + 1], maxd[NMAX + 1];
  u64 thr_cnt[NTHR], thr_maxd[NTHR];
  u64 entries[4], total_entries;
  CC_SHA256_CTX ctx[NTHR];
} state_t;

static state_t S;
static uint8_t *bad; /* bit i set: odd p^2 divides lo + i */
static u64 lo;

static void sieve_squares(u64 k) {
  lo = k * SLAB;
  u64 hi = lo + SLAB; /* exclusive */
  size_t nb = SLAB / 8;
  bad = calloc(nb, 1);
  u64 lim = 1;
  while (lim * lim < hi) lim++;
  char *comp = calloc(lim + 1, 1);
  for (u64 p = 3; p <= lim; p += 2) {
    if (comp[p]) continue;
    for (u64 j = p * p; j <= lim; j += 2 * p) comp[j] = 1;
    u64 q = p * p;
    if (q >= hi) break;
    u64 first = ((lo + q - 1) / q) * q;
    for (u64 x = first; x < hi; x += q) {
      u64 i = x - lo;
      bad[i >> 3] |= (uint8_t)(1u << (i & 7));
    }
  }
  free(comp);
}
static inline int is_fund(u64 d) { /* d in the slab, residue class already valid */
  u64 i = d - lo;
  if (d == 0) return 0;
  return !(bad[i >> 3] & (1u << (i & 7)));
}

static u64 dump_lo = 0, dump_hi = 0; static FILE *dumpf = 0; /* optional: all (d, h) with dump_lo <= d < dump_hi (file dump.cfg: "lo hi path") */
static u64 *arr[4];
static size_t alen[4], acap[4];

static void die(const char *m, const char *f, u64 x) {
  fprintf(stderr, "FAIL %s %s %llu\n", m, f, x);
  exit(2);
}

static void process_file(int c, u64 k, const char *path) {
  gzFile g = gzopen(path, "rb");
  if (!g) die("open", path, 0);
  gzbuffer(g, 1 << 20);
  const u64 m = CMOD[c], r = CR[c];
  u64 d = k * SLAB + r;           /* d_0 */
  u64 e = k * SLAB + r;           /* next expected fundamental d (candidate) */
  while (e < lo + SLAB && !is_fund(e)) e += m;
  size_t BUFSZ = 1 << 24;
  char *buf = malloc(BUFSZ + 1);
  size_t have = 0;
  u64 n = 0, first = 0, nle[NTHR] = {0};
  int eof = 0;
  for (;;) {
    if (!eof) {
      int got = gzread(g, buf + have, (unsigned)(BUFSZ - have));
      if (got < 0) { int en; die(gzerror(g, &en), path, n); }
      if (got == 0) eof = 1;
      have += (size_t)got;
    }
    if (have == 0 && eof) break;
    /* process complete lines */
    size_t pos = 0;
    for (;;) {
      char *nl = memchr(buf + pos, '\n', have - pos);
      if (!nl) break;
      char *p = buf + pos, *end = nl;
      u64 vals[64]; int nv = 0;
      while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t' || *p == '\r')) p++;
        if (p >= end) break;
        u64 v = 0;
        if (*p < '0' || *p > '9') die("nondigit", path, n);
        while (p < end && *p >= '0' && *p <= '9') { v = v * 10 + (u64)(*p - '0'); p++; }
        if (nv >= 64) die("too many fields", path, n);
        vals[nv++] = v;
      }
      if (nv < 2) die("short line", path, n);
      d += m * vals[0];
      u64 h = vals[1];
      /* invariants */
      u64 prod = 1;
      for (int i = 2; i < nv; i++) {
        prod *= vals[i];
        if (i > 2 && vals[i - 1] % vals[i] != 0) die("chain", path, d);
      }
      if (nv == 2) { if (h != 1) die("no invariants but h>1", path, d); }
      else if (prod != h) die("product", path, d);
      /* completeness / position */
      if (d != e) die("unexpected d (expected next fundamental)", path, d);
      e += m;
      while (e < lo + SLAB && !is_fund(e)) e += m;
      if (n == 0) first = d;
      n++;
      if (dumpf && d >= dump_lo && d < dump_hi) fprintf(dumpf, "%llu %llu\n", d, h);
      if (h <= NMAX) {
        if (alen[c] == acap[c]) { acap[c] = acap[c] ? 2 * acap[c] : (1 << 20); arr[c] = realloc(arr[c], acap[c] * 8); if (!arr[c]) die("oom", path, 0); }
        arr[c][alen[c]++] = (d << 14) | h;
        for (int j = 0; j < NTHR; j++) if (h <= (u64)THR[j]) nle[j]++;
      }
      pos = (size_t)(nl - buf) + 1;
    }
    memmove(buf, buf + pos, have - pos);
    have -= pos;
    if (eof) { if (have) die("trailing partial line", path, n); break; }
  }
  int en; const char *em = gzerror(g, &en);
  if (en != Z_OK && en != Z_STREAM_END) die(em, path, n);
  gzclose(g);
  free(buf);
  if (e < lo + SLAB) die("missing fundamental d at end of file", path, e);
  if (n == 0) die("empty", path, 0);
  S.entries[c] += n; S.total_entries += n;
  printf("%s.%llu: entries %llu, first %llu, last %llu, h<=2000: %llu, h<=%d: %llu\n", CNAME[c], k, n, first, d, nle[1], NMAX, nle[NTHR - 1]);
}

#define SHABUF (1 << 20)
static char *sbuf[NTHR];
static size_t slen[NTHR];
static void merge_update(void) {
  size_t idx[4] = {0, 0, 0, 0};
  char line[64];
  u64 prev = 0;
  for (int j = 0; j < NTHR; j++) { if (!sbuf[j]) sbuf[j] = malloc(SHABUF); slen[j] = 0; }
  for (;;) {
    int best = -1; u64 bv = 0;
    for (int c = 0; c < 4; c++) if (idx[c] < alen[c] && (best < 0 || arr[c][idx[c]] < bv)) { best = c; bv = arr[c][idx[c]]; }
    if (best < 0) break;
    idx[best]++;
    u64 d = bv >> 14, h = bv & 16383;
    if (d <= prev && prev) die("merge order", "", d);
    prev = d;
    S.cnt[h]++; S.maxd[h] = d;
    int L = snprintf(line, sizeof line, "%llu %llu\n", d, h);
    for (int j = 0; j < NTHR; j++) if (h <= (u64)THR[j]) {
      if (slen[j] + (size_t)L > SHABUF) { CC_SHA256_Update(&S.ctx[j], sbuf[j], (CC_LONG)slen[j]); slen[j] = 0; }
      memcpy(sbuf[j] + slen[j], line, (size_t)L); slen[j] += (size_t)L;
      S.thr_cnt[j]++; S.thr_maxd[j] = d;
    }
  }
  for (int j = 0; j < NTHR; j++) if (slen[j]) { CC_SHA256_Update(&S.ctx[j], sbuf[j], (CC_LONG)slen[j]); slen[j] = 0; }
}

static void load_state(const char *path, u64 k) {
  FILE *f = fopen(path, "rb");
  if (!f) {
    if (k != 0) die("no state but K != 0", path, k);
    memset(&S, 0, sizeof S);
    S.magic = 0x6d6a7831; S.nmax = NMAX; S.next_k = 0;
    for (int j = 0; j < NTHR; j++) CC_SHA256_Init(&S.ctx[j]);
    return;
  }
  if (fread(&S, sizeof S, 1, f) != 1) die("state read", path, 0);
  fclose(f);
  if (S.magic != 0x6d6a7831 || S.nmax != NMAX) die("state magic", path, 0);
  if (S.next_k != k) die("state next_k != K", path, S.next_k);
}
static void save_state(const char *path) {
  char tmp[4096]; snprintf(tmp, sizeof tmp, "%s.tmp", path);
  FILE *f = fopen(tmp, "wb");
  if (!f || fwrite(&S, sizeof S, 1, f) != 1 || fclose(f)) die("state write", tmp, 0);
  if (rename(tmp, path)) die("rename", path, 0);
}

int main(int argc, char **argv) {
  if (argc == 4 && !strcmp(argv[1], "--report")) {
    FILE *f = fopen(argv[2], "rb");
    if (!f || fread(&S, sizeof S, 1, f) != 1) die("state", argv[2], 0);
    fclose(f);
    char p1[4096], p2[4096];
    snprintf(p1, sizeof p1, "%s_perh.txt", argv[3]); snprintf(p2, sizeof p2, "%s_sha.txt", argv[3]);
    FILE *a = fopen(p1, "w"), *b = fopen(p2, "w");
    fprintf(a, "# h  number_of_fields  largest_d   (fields with d < %llu*2^28 from the Mosunov-Jacobson tables)\n", S.next_k);
    for (int h = 1; h <= NMAX; h++) fprintf(a, "%d %llu %llu\n", h, S.cnt[h], S.maxd[h]);
    fprintf(b, "# slabs k < %llu (d < %llu), total entries %llu (cl3mod8 %llu, cl7mod8 %llu, cl4mod16 %llu, cl8mod16 %llu)\n", S.next_k, S.next_k * SLAB, S.total_entries, S.entries[0], S.entries[1], S.entries[2], S.entries[3]);
    fprintf(b, "# N  fields_with_h<=N  largest_d  sha256(sorted lines \"d h\\n\")\n");
    for (int j = 0; j < NTHR; j++) {
      CC_SHA256_CTX c = S.ctx[j]; unsigned char md[32];
      CC_SHA256_Final(md, &c);
      fprintf(b, "%d %llu %llu ", THR[j], S.thr_cnt[j], S.thr_maxd[j]);
      for (int i = 0; i < 32; i++) fprintf(b, "%02x", md[i]);
      fprintf(b, "\n");
    }
    fclose(a); fclose(b);
    return 0;
  }
  if (argc != 4) { fprintf(stderr, "usage: mjx STATE K DIR | mjx --report STATE OUTPREFIX\n"); return 1; }
  u64 k = strtoull(argv[2], 0, 10);
  load_state(argv[1], k);
  { FILE *cf = fopen("dump.cfg", "r"); char dp[4096];
    if (cf) { if (fscanf(cf, "%llu %llu %4095s", &dump_lo, &dump_hi, dp) == 3 && dump_hi > k * SLAB && dump_lo < (k + 1) * SLAB) dumpf = fopen(dp, "a"); fclose(cf); } }
  sieve_squares(k);
  for (int c = 0; c < 4; c++) {
    char path[4096]; snprintf(path, sizeof path, "%s/%s.%llu.gz", argv[3], CNAME[c], k);
    process_file(c, k, path);
  }
  if (dumpf) fclose(dumpf);
  merge_update();
  S.next_k = k + 1;
  save_state(argv[1]);
  printf("slab %llu done: total entries %llu; h<=2000: %llu (largest %llu); h<=%d: %llu (largest %llu)\n", k, S.total_entries, S.thr_cnt[1], S.thr_maxd[1], NMAX, S.thr_cnt[NTHR - 1], S.thr_maxd[NTHR - 1]);
  return 0;
}
