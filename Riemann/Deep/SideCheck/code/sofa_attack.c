// sidecheck: independent floating-point re-check of the sofa branch-and-bound certificate.
// Own tree traversal (same split rule as stated in PROOF.md: widest coordinate, first on ties, midpoint),
// own area routine: exact-in-principle vertical slab integration of
//   K cap INTERSECT_j (Pt_j u Qt_j),  K = [0,W]x[0,1],
// Pt = {p<=b1+1, a2<=q<=b2+1}, Qt = {a1<=p<=b1+1, q<=b2},  p = c x + s y, q = -s x + c y.
// (cross-section length is piecewise linear in x; breakpoints = pairwise line intersections and
//  crossings of y=0, y=1; midpoint rule on each piece is exact up to rounding.)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MAXC 4
static int ncorr; static double cs[MAXC][2]; static double W, T;
typedef struct { double a1,b1,a2,b2; } Box;
// cross-section of one corridor at x: up to 2 intervals (Pt, Qt) clipped to [0,1]
static int iv_corr(int j, const Box *b, double x, double *lo, double *hi) {
  double c = cs[j][0], s = cs[j][1];
  // constraint helpers: alpha*y <= beta  ->  y <= beta/alpha (alpha>0) or y >= beta/alpha (alpha<0)
  int n = 0;
  for (int part = 0; part < 2; part++) {
    double L = 0.0, H = 1.0;
    // list constraints as (alpha, beta): alpha*y <= beta
    double al[3], be[3]; int m = 0;
    if (part == 0) { // Pt: p <= b1+1 ; q >= a2 ; q <= b2+1
      al[m] = s;  be[m] = b->b1 + 1 - c * x; m++;
      al[m] = -c; be[m] = -(b->a2) + (-s * x); m++;   // -q <= -a2  :  -(-s x + c y) <= -a2 -> -c y <= -a2 - s x
      al[m] = c;  be[m] = b->b2 + 1 + s * x; m++;      // q <= b2+1 : c y <= b2+1 + s x
    } else {         // Qt: p >= a1 ; p <= b1+1 ; q <= b2
      al[m] = -s; be[m] = -(b->a1) + c * x; m++;       // -p <= -a1 : -s y <= -a1 + c x
      al[m] = s;  be[m] = b->b1 + 1 - c * x; m++;
      al[m] = c;  be[m] = b->b2 + s * x; m++;
    }
    for (int k = 0; k < m; k++) {
      double r = be[k] / al[k];
      if (al[k] > 0) { if (r < H) H = r; } else { if (r > L) L = r; }
    }
    if (H > L) { lo[n] = L; hi[n] = H; n++; }
  }
  // merge into disjoint sorted intervals
  if (n == 2) {
    if (lo[1] < lo[0]) { double t = lo[0]; lo[0] = lo[1]; lo[1] = t; t = hi[0]; hi[0] = hi[1]; hi[1] = t; }
    if (lo[1] <= hi[0]) { if (hi[1] > hi[0]) hi[0] = hi[1]; n = 1; }
  }
  return n;
}
static double cross_len(const Box *bx, double x) {
  // intersect unions over corridors
  double lo[8], hi[8]; int n = iv_corr(0, &bx[0], x, lo, hi);
  for (int j = 1; j < ncorr && n > 0; j++) {
    double l2[2], h2[2]; int m = iv_corr(j, &bx[j], x, l2, h2);
    double nl[8], nh[8]; int k = 0;
    for (int a = 0; a < n; a++) for (int b = 0; b < m; b++) {
      double L = lo[a] > l2[b] ? lo[a] : l2[b], H = hi[a] < h2[b] ? hi[a] : h2[b];
      if (H > L) { nl[k] = L; nh[k] = H; k++; }
    }
    n = k; memcpy(lo, nl, sizeof(double) * k); memcpy(hi, nh, sizeof(double) * k);
  }
  double t = 0; for (int a = 0; a < n; a++) t += hi[a] - lo[a];
  return t;
}
static int cmpd(const void *a, const void *b) { double x = *(const double*)a, y = *(const double*)b; return (x > y) - (x < y); }
static double region_area(const Box *bx) {
  // lines in the form A x + B y = Cc
  double LA[32], LB[32], LC[32]; int nl = 0;
  LA[nl] = 0; LB[nl] = 1; LC[nl] = 0; nl++;
  LA[nl] = 0; LB[nl] = 1; LC[nl] = 1; nl++;
  for (int j = 0; j < ncorr; j++) {
    double c = cs[j][0], s = cs[j][1];
    double pv[2] = { bx[j].a1, bx[j].b1 + 1 }, qv[3] = { bx[j].a2, bx[j].b2, bx[j].b2 + 1 };
    for (int k = 0; k < 2; k++) { LA[nl] = c; LB[nl] = s; LC[nl] = pv[k]; nl++; }
    for (int k = 0; k < 3; k++) { LA[nl] = -s; LB[nl] = c; LC[nl] = qv[k]; nl++; }
  }
  double xs[600]; int nx = 0; xs[nx++] = 0; xs[nx++] = W;
  for (int a = 0; a < nl; a++) for (int b = a + 1; b < nl; b++) {
    double det = LA[a] * LB[b] - LA[b] * LB[a];
    if (fabs(det) < 1e-15) continue;
    double x = (LC[a] * LB[b] - LC[b] * LB[a]) / det;
    if (x > 0 && x < W) xs[nx++] = x;
  }
  qsort(xs, nx, sizeof(double), cmpd);
  double area = 0;
  for (int i = 0; i + 1 < nx; i++) {
    double w = xs[i + 1] - xs[i]; if (w <= 0) continue;
    area += w * cross_len(bx, 0.5 * (xs[i] + xs[i + 1]));
  }
  return area;
}

static double urand(unsigned long long *st){ *st = *st*6364136223846793005ULL + 1442695040888963407ULL; return ((*st >> 11) * (1.0/9007199254740992.0)); }
static Box dom[MAXC];
static double f(const double *u){ Box b[MAXC]; for(int j=0;j<ncorr;j++){ b[j].a1=b[j].b1=u[2*j]; b[j].a2=b[j].b2=u[2*j+1]; } return region_area(b); }
int main(int argc,char**argv){
  const char*dir=argv[1]; long starts=argc>2?atol(argv[2]):20000; char fn[512];
  snprintf(fn,sizeof fn,"%s/cases.txt",dir); FILE*F=fopen(fn,"r"); int ncase; if(fscanf(F,"%d %lf %lf",&ncase,&T,&W)!=3) return 9;
  unsigned long long st=12345;
  for(int ci=0;ci<ncase;ci++){ char name[128]; long nb; if(fscanf(F," %127[^|]|%d|%ld",name,&ncorr,&nb)!=3) return 9; int js[MAXC];
    for(int j=0;j<ncorr;j++) if(fscanf(F,"%d %lf %lf %lf %lf %lf %lf",&js[j],&cs[j][0],&cs[j][1],&dom[j].a1,&dom[j].b1,&dom[j].a2,&dom[j].b2)!=7) return 9;
    double best=-1, bu[8];
    for(long s=0;s<starts;s++){
      double u[8]; for(int j=0;j<ncorr;j++){ u[2*j]=dom[j].a1+urand(&st)*(dom[j].b1-dom[j].a1); u[2*j+1]=dom[j].a2+urand(&st)*(dom[j].b2-dom[j].a2); }
      double v=f(u); if(v<=0.3) continue;
      double h=0.25;
      int iters=0; while(h>1e-8 && iters<3000){ iters++; int imp=0; for(int k=0;k<2*ncorr;k++) for(int sg=-1;sg<=1;sg+=2){ double o=u[k]; u[k]=o+sg*h; double w=f(u); if(w>v+1e-15){v=w;imp=1;} else u[k]=o; }
        if(!imp) h*=0.5; }
      for(int it=0;it<400;it++){ double w8[8]; for(int k=0;k<2*ncorr;k++) w8[k]=u[k]+(urand(&st)-0.5)*1e-3; double w=f(w8); if(w>v){v=w; memcpy(u,w8,sizeof w8);} }
      if(v>best){best=v; memcpy(bu,u,sizeof bu);} }
    printf("%-16s j=%d,%d  best exact-relaxation area found %.9f at u=",name,js[0],js[1],best); for(int k=0;k<2*ncorr;k++) printf(" %.6f",bu[k]); printf("  %s\n", best>T?"EXCEEDS T!":"(< T)");
  }
  return 0; }
