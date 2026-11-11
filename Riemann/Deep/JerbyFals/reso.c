// reso.c -- cheap large-|zeta| finder: maximise F(t) = sum_{p<=P} p^(-1/2) cos(t log p) + (1/2) sum_{p<=7} p^(-1) cos(2 t log p)
// on a grid t = T + j h, j < H/h (phases by rotation, re-anchored in double-double every 4096 steps), keep the best grid
// value per bin of length SEP, print the K best bins as Gram windows "id n0 cnt" (cnt Gram points centred on the peak).
// usage: reso T H h P SEP K CNT TAG
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct { double hi, lo; } dd;
static inline dd mk(double a, double b){ dd r; r.hi=a; r.lo=b; return r; }
static inline dd two_sum(double a, double b){ double s=a+b, bb=s-a; return mk(s,(a-(s-bb))+(b-bb)); }
static inline dd qts(double a, double b){ double s=a+b; return mk(s,b-(s-a)); }
static inline dd two_prod(double a, double b){ double p=a*b; return mk(p,fma(a,b,-p)); }
static inline dd dd_add(dd a, dd b){ dd s=two_sum(a.hi,b.hi), t=two_sum(a.lo,b.lo); s.lo+=t.hi; s=qts(s.hi,s.lo); s.lo+=t.lo; return qts(s.hi,s.lo); }
static inline dd dd_mul(dd a, dd b){ dd p=two_prod(a.hi,b.hi); p.lo+=a.hi*b.lo+a.lo*b.hi; return qts(p.hi,p.lo); }
static inline dd dd_muld(dd a, double b){ dd p=two_prod(a.hi,b); p.lo+=a.lo*b; return qts(p.hi,p.lo); }
static inline dd dd_divd(dd a, double b){ double q=a.hi/b; dd p=two_prod(q,b); double e=((a.hi-p.hi)-p.lo+a.lo)/b; return qts(q,e); }
static const dd LN2 = {6.93147180559945286e-01, 2.31904681384629956e-17};
static const dd INV2PI = {1.59154943091895346e-01, -9.83933833759124294e-18};
static const double PI=3.14159265358979323846, TWOPI=6.28318530717958647692;
static dd dd_exp(dd a){
  double k=rint(a.hi/LN2.hi); dd r=dd_add(a, dd_muld(LN2,-k)); r.hi=ldexp(r.hi,-10); r.lo=ldexp(r.lo,-10);
  dd s=r, term=r; for(int i=2;i<=10;i++){ term=dd_divd(dd_mul(term,r),(double)i); s=dd_add(s,term); }
  for(int i=0;i<10;i++) s=dd_add(dd_muld(s,2.0), dd_mul(s,s));
  s=dd_add(s,mk(1.0,0.0)); return mk(ldexp(s.hi,(int)k), ldexp(s.lo,(int)k));
}
static dd dd_log(double x){ double y=log(x); dd e=dd_exp(mk(-y,0.0)); dd r=dd_add(dd_muld(e,x),mk(-1.0,0.0)); return dd_add(mk(y,0.0),r); }
static double theta_o(double t){ return t/2*log(t/TWOPI)-t/2-PI/8+1/(48*t); }
int main(int argc, char **argv){
  double T=atof(argv[1]), H=atof(argv[2]), h=atof(argv[3]); int P=atoi(argv[4]); double SEP=atof(argv[5]); int K=atoi(argv[6]); long CNT=atol(argv[7]); const char *tag=argv[8];
  int np=0; double fr[64], w[64]; double L[64];   /* frequencies (turns per unit t) and weights */
  dd Ld[64];
  for(int p=2;p<=P;p++){ int pr=1; for(int q=2;q*q<=p;q++) if(p%q==0){pr=0;break;} if(!pr) continue;
    Ld[np]=dd_mul(dd_log((double)p),INV2PI); L[np]=Ld[np].hi; w[np]=1/sqrt((double)p); np++;
    if(p<=7){ Ld[np]=dd_muld(Ld[np-1],2.0); L[np]=Ld[np].hi; w[np]=0.5/p; np++; } }
  long nsteps=(long)(H/h); long nbins=(long)(H/SEP)+1;
  double *bv=malloc(sizeof(double)*nbins); double *bt=malloc(sizeof(double)*nbins); for(long b=0;b<nbins;b++){bv[b]=-1e9;bt[b]=0;}
  double zr[64], zi[64], wr[64], wi[64], ph0[64];
  for(int i=0;i<np;i++){ dd q=two_prod(T,Ld[i].hi); double f=(q.hi-rint(q.hi))+(q.lo+T*Ld[i].lo); ph0[i]=f-rint(f); wr[i]=cos(TWOPI*h*L[i]); wi[i]=-sin(TWOPI*h*L[i]); }
  for(long j=0;j<nsteps;j++){
    if((j&4095)==0){ for(int i=0;i<np;i++){ double f=ph0[i]+((double)j*h)*L[i]; f-=rint(f); zr[i]=cos(TWOPI*f); zi[i]=-sin(TWOPI*f); } }
    double F=0; for(int i=0;i<np;i++){ F+=w[i]*zr[i]; double a=zr[i]*wr[i]-zi[i]*wi[i]; zi[i]=zr[i]*wi[i]+zi[i]*wr[i]; zr[i]=a; }
    long b=(long)((double)j*h/SEP); if(F>bv[b]){bv[b]=F; bt[b]=T+(double)j*h;}
  }
  /* top K bins */
  double Fmax_tot=0; for(int i=0;i<np;i++) Fmax_tot+=w[i];
  for(int k=0;k<K;k++){ long best=-1; for(long b=0;b<nbins;b++) if(best<0||bv[b]>bv[best]) best=b;
    double t=bt[best]; long long n=(long long)floor(theta_o(t)/PI);
    printf("R%s_%02d %lld %ld\n", tag, k, n-CNT/2, CNT);
    fprintf(stderr,"%s k=%d t=%.3f F=%.4f (max %.4f) n=%lld\n", tag, k, t, bv[best], Fmax_tot, n);
    bv[best]=-1e9; if(best>0) bv[best-1]=-1e9; if(best+1<nbins) bv[best+1]=-1e9; }
  return 0;
}
