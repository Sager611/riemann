// gram_win.c -- windowed scan for Jerby's Gram-point repulsion conjecture (2406.06548 conj:rep; 2310.14415 "Repulsion G-B-G"):
// g_n bad ((-1)^n Z(g_n) < 0) with g_{n-1}, g_{n+1} good  =>  mu(g_n) = |Z'(g_n)/Z(g_n)| > 4.
// Riemann-Siegel main sum + C0 remainder (error ~ C1 (t/2pi)^(-3/4) < 1e-6 for t >= 1e8).
// Precision: the phases theta(t)/2pi - t log(k)/2pi (size ~1e12 turns at t ~ 5e11) are formed relative to a block centre tc
// (a double): theta(tc)/2pi and frac(tc log(k)/2pi) are computed once per block of BS Gram points in double-double,
// the offset delta = t - tc (|delta| < 400) in double.  Phase error ~1e-13 turns (plain double at t = 1e11: ~1e-4 turns).
// usage:  gram_win WINFILE THR [NTHR [G]] WINFILE lines "id n0 count"; prints WIN/HIST/LOW summaries and CAND lines (mu < THR or
//                                          mu/theta' < NTHR); G > 0: zero count of every isolated bad block on G points per interval (ZF line)
//         gram_win -check n [n ...]       values at single Gram points (dd method and the old plain-double method)
//         gram_win -reso T H h P SEP Q M THR NTHR TAG   resonance-prefiltered scan of t in [T, T+H] (see reso_scan)
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

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
static const dd LOG2PI = {1.83787706640934556e+00, -7.75658831613448290e-17};
static const dd INV4PI = {7.95774715459476728e-02, -4.91966916879562147e-18};
static const dd INV2PI = {1.59154943091895346e-01, -9.83933833759124294e-18};
static const double PI=3.14159265358979323846, TWOPI=6.28318530717958647692;

static dd dd_exp(dd a){
  double k=rint(a.hi/LN2.hi); dd r=dd_add(a, dd_muld(LN2,-k));
  r.hi=ldexp(r.hi,-10); r.lo=ldexp(r.lo,-10);
  dd s=r, term=r;
  for(int i=2;i<=10;i++){ term=dd_divd(dd_mul(term,r),(double)i); s=dd_add(s,term); }
  for(int i=0;i<10;i++) s=dd_add(dd_muld(s,2.0), dd_mul(s,s));      /* (1+s)^2 - 1 */
  s=dd_add(s,mk(1.0,0.0));
  return mk(ldexp(s.hi,(int)k), ldexp(s.lo,(int)k));
}
static dd dd_log(double x){ double y=log(x); dd e=dd_exp(mk(-y,0.0)); dd r=dd_add(dd_muld(e,x),mk(-1.0,0.0)); return dd_add(mk(y,0.0),r); }
static dd Theta_dd(double t){   /* theta(t)/(2 pi), asymptotic series to t^-3 */
  dd L=dd_add(dd_log(t), mk(-LOG2PI.hi,-LOG2PI.lo)); L=dd_add(L, mk(-1.0,0.0));
  dd T=dd_mul(dd_muld(INV4PI, t), L);
  double c=-0.0625 + (1.0/(48*t) + 7.0/(5760*t*t*t))/TWOPI;
  return dd_add(T, mk(c,0.0));
}

static int KMAX; static double *L2h, *L2l, *rs, *ph, *lgk;
static double C[10], S[10];
#define POLC(z) (C[0]+z*(C[1]+z*(C[2]+z*(C[3]+z*(C[4]+z*(C[5]+z*(C[6]+z*(C[7]+z*C[8]))))))))
#define POLS(z) (S[0]+z*(S[1]+z*(S[2]+z*(S[3]+z*(S[4]+z*(S[5]+z*(S[6]+z*(S[7]+z*S[8]))))))))
static void init_tables(int kmax){
  KMAX=kmax; L2h=malloc(sizeof(double)*(kmax+2)); L2l=malloc(sizeof(double)*(kmax+2)); rs=malloc(sizeof(double)*(kmax+2));
  ph=malloc(sizeof(double)*(kmax+2)); lgk=malloc(sizeof(double)*(kmax+2));
  for(int k=1;k<=kmax+1;k++){ dd l=dd_log((double)k); dd q=dd_mul(l,INV2PI); L2h[k]=q.hi; L2l[k]=q.lo; rs[k]=1/sqrt((double)k); lgk[k]=log((double)k); }
  double f=1; for(int k=0;k<9;k++){ C[k]=((k&1)?-1:1)*pow(TWOPI,2*k)/f; f*=(2*k+1); S[k]=((k&1)?-1:1)*pow(TWOPI,2*k+1)/f; f*=(2*k+2); }
}
typedef struct { double tc, fc; long long Ic; double T1,T2,T3,T4; int Nmax; } Blk;
static Blk B;
static void setup_block(double tc, double tmax){
  B.tc=tc; dd Th=Theta_dd(tc); double fl=floor(Th.hi); double f=(Th.hi-fl)+Th.lo; long long I=(long long)fl;
  if(f<0){f+=1;I--;} if(f>=1){f-=1;I++;}
  B.fc=f; B.Ic=I;
  double lg=log(tc/TWOPI);
  B.T1=(0.5*lg - 1/(48*tc*tc))/TWOPI; B.T2=(0.5/tc + 1/(24*tc*tc*tc))/TWOPI; B.T3=(-0.5/(tc*tc))/TWOPI; B.T4=(1/(tc*tc*tc))/TWOPI;
  int N=(int)sqrt(tmax/TWOPI)+2; if(N>KMAX){ fprintf(stderr,"KMAX too small %d\n",N); exit(1);} B.Nmax=N;
  for(int k=1;k<=N;k++){ dd p=two_prod(tc,L2h[k]); double f2=(p.hi-rint(p.hi)) + (p.lo + tc*L2l[k]); f2-=rint(f2); ph[k]=f2; }
}
static inline double dTh(double d){ return d*(B.T1 + d*(B.T2/2 + d*(B.T3/6 + d*B.T4/24))); }
static inline double Thp(double d){ return B.T1 + d*(B.T2 + d*(B.T3/2 + d*B.T4/6)); }   /* theta'(t)/(2 pi) */
static double gram_delta(long long n){   /* delta with theta(tc + delta) = pi n */
  double Dn=0.5*(double)(n-2*B.Ic) - B.fc; double d=Dn/B.T1;
  for(int i=0;i<4;i++) d-=(dTh(d)-Dn)/Thp(d);
  return d;
}
static double C0f(double p){ return cos(TWOPI*(p*p-p-1.0/16))/cos(TWOPI*p); }
static double C0(double p){ double c=cos(TWOPI*p); if(fabs(c)<1e-3){ double h=2e-3; return 0.5*(C0f(p-h)+C0f(p+h)) - (C0f(p-2*h)+C0f(p+2*h)-C0f(p-h)-C0f(p+h))/6.0; } return C0f(p); }
static double Zval(double d){
  double t=B.tc+d, a=sqrt(t/TWOPI); int N=(int)a; double p=a-N;
  double A=B.fc+dTh(d); A-=rint(A);
  double s=0;
  for(int k=1;k<=N;k++){ double y=A-ph[k]-d*L2h[k]; y-=rint(y); double ay=fabs(y); double sg=ay>0.25?-1.0:1.0; double yr=ay>0.25?0.5-ay:ay; double z=yr*yr; s+=rs[k]*sg*POLC(z); }
  return 2*s + ((N-1)&1?-1.0:1.0)*pow(a,-0.5)*C0(p);
}
static double Zder(double d){
  double t=B.tc+d, a=sqrt(t/TWOPI); int N=(int)a; double p=a-N;
  double A=B.fc+dTh(d); A-=rint(A); double tp=Thp(d);
  double s=0;
  for(int k=1;k<=N;k++){ double y=A-ph[k]-d*L2h[k]; y-=rint(y); double ay=fabs(y); double sg=y<0?-1.0:1.0; double yr=ay>0.25?0.5-ay:ay; double z=yr*yr; s+=rs[k]*(tp-L2h[k])*sg*yr*POLS(z); }
  double sig=((N-1)&1?-1.0:1.0), h=1e-6;
  double Rp=sig*( pow(a,-0.5)*(C0(p+h)-C0(p-h))/(2*h) - 0.5*pow(a,-1.5)*C0(p) )/(4*PI*a);
  return -2*TWOPI*s + Rp;
}
/* ---- old plain-double method (ConjFals gram_fast2.c) for comparison ---- */
static double theta_o(double t){ return t/2*log(t/TWOPI)-t/2-PI/8+1/(48*t)+7/(5760*t*t*t); }
static double dtheta_o(double t){ return 0.5*log(t/TWOPI)-1/(48*t*t); }
static double gram_est(long long n){   /* plain double Gram point (accurate to ~1e-4 at 1e11) */
  double x=((double)n+0.125)/exp(1.0), w=log(x+1); for(int i=0;i<60;i++) w-=(w*exp(w)-x)/(exp(w)*(w+1));
  double t=TWOPI*exp(1.0)*exp(w); for(int i=0;i<6;i++) t-=(theta_o(t)-PI*(double)n)/dtheta_o(t);
  return t;
}
static double Zval_o(double t){
  double a=sqrt(t/TWOPI); int N=(int)a; double p=a-N; double th2=theta_o(t)/TWOPI; double s=0;
  for(int k=1;k<=N;k++){ double y=th2-t*L2h[k]; y-=rint(y); double ay=fabs(y); double sg=ay>0.25?-1.0:1.0; double yr=ay>0.25?0.5-ay:ay; double z=yr*yr; s+=rs[k]*sg*POLC(z); }
  return 2*s+((N-1)&1?-1.0:1.0)*pow(a,-0.5)*C0(p);
}
static double Zder_o(double t){
  double a=sqrt(t/TWOPI); int N=(int)a; double th=theta_o(t), thp=dtheta_o(t); double sp=0;
  for(int k=1;k<=N;k++){ double x=th-t*lgk[k]; x-=TWOPI*nearbyint(x/TWOPI); sp+=rs[k]*(thp-lgk[k])*sin(x); }
  return -2*sp;
}
#define BS 1024
static double ZZ[BS+8], DL[BS+8];
static const double AT[9]={4.0,4.25,4.5,4.75,5.0,5.5,6.0,7.0,8.0};
static const double NT[9]={0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.80};
static double NTHR=0; static int ZG=0;   /* ZG > 0: count zeros of every isolated bad block on ZG grid points per Gram interval */
#define NLOW 20
static void ins(double *a, double v){ if(v>=a[NLOW-1]) return; int j=NLOW-1; while(j>0 && a[j-1]>v){ a[j]=a[j-1]; j--; } a[j]=v; }
static int zcount(double d0, double d1, double z0, double z1, int G){ int c=0; double prev=z0; for(int q=1;q<=G;q++){ double v=(q<G)?Zval(d0+(d1-d0)*q/(double)G):z1; if(v*prev<0) c++; prev=v; } return c; }
static void window(const char *id, long long n0, long long cnt, double THR){
  clock_t c0=clock(); long long nend=n0+cnt;
  long long nbad=0, niso=0, ca[9]={0}, cn[9]={0}, nbig=0;
  double minmu=1e99, minmn=1e99, maxz=0, t_first=0, t_last=0; long long am=-1, an=-1, amz=-1;
  double am_t=0, am_Z=0, am_Zp=0, am_Zl=0, am_Zr=0, an_mu=0;
  double lowmu[NLOW], lownm[NLOW]; for(int j=0;j<NLOW;j++){ lowmu[j]=1e99; lownm[j]=1e99; }
  long long nzf=0, nros=0, noth=0, czf[9]={0}, cros[9]={0}; double minzf=1e99, minros=1e99;
  for(long long nb=n0; nb<nend; nb+=BS){
    long long ne=nb+BS<nend?nb+BS:nend, nm=(nb+ne)/2;
    double tc=gram_est(nm), sp=TWOPI/log(tc/TWOPI);
    setup_block(tc, tc+(double)(ne-nm+4)*sp+1.0);
    for(long long n=nb-2;n<ne+2;n++){ int i=(int)(n-(nb-2)); DL[i]=gram_delta(n); ZZ[i]=Zval(DL[i]); }
    if(nb==n0) t_first=B.tc+DL[2];
    t_last=B.tc+DL[ne-1-(nb-2)];
    for(long long m=nb;m<ne;m++){
      int i=(int)(m-(nb-2)); double Zm=ZZ[i], sg=(m&1)?-1.0:1.0;
      if(fabs(Zm)>maxz){maxz=fabs(Zm); amz=m;}
      if(sg*Zm<0){ nbad++;
        if((-sg)*ZZ[i-1]>0 && (-sg)*ZZ[i+1]>0){
          niso++; double zp=Zder(DL[i]); double mu=fabs(zp/Zm); double thp=TWOPI*Thp(DL[i]); double mn=mu/thp;
          if(fabs(Zm)>10) nbig++;
          for(int j=0;j<9;j++){ if(mu<AT[j]) ca[j]++; if(mn<NT[j]) cn[j]++; }
          if(mu<minmu){ minmu=mu; am=m; am_t=B.tc+DL[i]; am_Z=Zm; am_Zp=zp; am_Zl=ZZ[i-1]; am_Zr=ZZ[i+1]; }
          if(mn<minmn){ minmn=mn; an=m; an_mu=mu; }
          ins(lowmu,mu); ins(lownm,mn);
          if(ZG>0){ int zl=zcount(DL[i-1],DL[i],ZZ[i-1],Zm,ZG), zr=zcount(DL[i],DL[i+1],Zm,ZZ[i+1],ZG);
            if(zl+zr==0){ nzf++; if(mn<minzf) minzf=mn; for(int j=0;j<9;j++) if(mn<NT[j]) czf[j]++; }
            else if(zl+zr==2 && (zl==0||zr==0)){ nros++; if(mn<minros) minros=mn; for(int j=0;j<9;j++) if(mn<NT[j]) cros[j]++; }
            else noth++; }
          if(mu<THR || mn<NTHR){
            /* zero count in (g_{m-1}, g_m) and (g_m, g_{m+1}) by sign changes on a grid of 24 points per interval */
            int zl=0, zr=0; double prev=ZZ[i-1];
            for(int q=1;q<=24;q++){ double v=(q<24)?Zval(DL[i-1]+(DL[i]-DL[i-1])*q/24.0):Zm; if(v*prev<0) zl++; prev=v; }
            prev=Zm; for(int q=1;q<=24;q++){ double v=(q<24)?Zval(DL[i]+(DL[i+1]-DL[i])*q/24.0):ZZ[i+1]; if(v*prev<0) zr++; prev=v; }
            printf("CAND %s n=%lld t=%.6f mu=%.6f norm=%.5f Z=%.6f Zp=%.6f Z[-2..2]=%.4f,%.4f,%.4f,%.4f,%.4f zeros_left=%d zeros_right=%d\n",
              id,m,B.tc+DL[i],mu,mn,Zm,zp,ZZ[i-2],ZZ[i-1],Zm,ZZ[i+1],ZZ[i+2],zl,zr); fflush(stdout);
          }
        }
      }
    }
  }
  double secs=(double)(clock()-c0)/CLOCKS_PER_SEC;
  printf("WIN %s n0=%lld cnt=%lld t=[%.3f,%.3f] thetap=%.5f bad=%lld iso=%lld minmu=%.6f n=%lld t=%.6f Z=%.6f Zp=%.6f Zl=%.6f Zr=%.6f minnorm=%.5f n=%lld (mu=%.5f) maxabsZ=%.4f n=%lld iso_bigZ=%lld secs=%.1f\n",
    id,n0,cnt,t_first,t_last,0.5*log(t_first/TWOPI),nbad,niso,minmu,am,am_t,am_Z,am_Zp,am_Zl,am_Zr,minmn,an,an_mu,maxz,amz,nbig,secs);
  printf("HIST %s mu<[4,4.25,4.5,4.75,5,5.5,6,7,8]:",id); for(int j=0;j<9;j++) printf(" %lld",ca[j]);
  printf(" norm<[.35,.40,.45,.50,.55,.60,.65,.70,.80]:"); for(int j=0;j<9;j++) printf(" %lld",cn[j]); printf("\n");
  printf("LOW %s mu:",id); for(int j=0;j<NLOW;j++) printf(" %.5f",lowmu[j]); printf(" | norm:"); for(int j=0;j<NLOW;j++) printf(" %.5f",lownm[j]); printf("\n");
  if(ZG>0){ printf("ZF %s G=%d zero_free=%lld rosser2=%lld other=%lld min_norm_zf=%.5f min_norm_ros=%.5f | zf norm<[.35..0.80]:",id,ZG,nzf,nros,noth,minzf,minros);
    for(int j=0;j<9;j++) printf(" %lld",czf[j]); printf(" | ros norm<:"); for(int j=0;j<9;j++) printf(" %lld",cros[j]); printf("\n"); }
  fflush(stdout);
}

/* ---- resonance-prefiltered scan: proxy F(t) = sum_{p<=P} p^-1/2 cos(t log p) + (1/2) sum_{p<=7} p^-1 cos(2 t log p) on a grid
   t = T + j h; bins of length SEP; every bin whose max F exceeds the Q-quantile of bin maxima (pilot pass over the first
   2e5 bins) gets an exact check of the Gram points n_c-M..n_c+M around the bin's argmax (n_c = nearest Gram index). ---- */
static long long R_bins=0, R_checked=0, R_gp=0, R_bad=0, R_iso=0, R_ca[9], R_cn[9], R_last=-1;
static double R_minmu=1e99, R_minmn=1e99, R_lowmu[NLOW], R_lownm[NLOW], R_maxz=0; static long long R_am=-1, R_an=-1, R_amz=-1;
static void reso_check(double ts, int M, double THR, const char *tag){
  double sp=TWOPI/log(ts/TWOPI);
  setup_block(ts, ts+(double)(M+3)*sp+1.0);
  long long nc=2*B.Ic+(long long)llround(2*B.fc), lo=nc-M, hi=nc+M;
  if(lo<=R_last) lo=R_last+1;
  if(lo>hi) return;
  R_checked++;
  for(long long n=lo-1;n<=hi+1;n++){ int i=(int)(n-(lo-1)); DL[i]=gram_delta(n); ZZ[i]=Zval(DL[i]); }
  for(long long m=lo;m<=hi;m++){
    int i=(int)(m-(lo-1)); double Zm=ZZ[i], sg=(m&1)?-1.0:1.0; R_gp++; if(fabs(Zm)>R_maxz){R_maxz=fabs(Zm);R_amz=m;}
    if(sg*Zm<0){ R_bad++;
      if((-sg)*ZZ[i-1]>0 && (-sg)*ZZ[i+1]>0){
        R_iso++; double zp=Zder(DL[i]); double mu=fabs(zp/Zm); double thp=TWOPI*Thp(DL[i]); double mn=mu/thp;
        for(int j=0;j<9;j++){ if(mu<AT[j]) R_ca[j]++; if(mn<NT[j]) R_cn[j]++; }
        if(mu<R_minmu){R_minmu=mu;R_am=m;} if(mn<R_minmn){R_minmn=mn;R_an=m;}
        ins(R_lowmu,mu); ins(R_lownm,mn);
        if(mu<THR || mn<NTHR){
          int zl=zcount(DL[i-1],DL[i],ZZ[i-1],Zm,24), zr=zcount(DL[i],DL[i+1],Zm,ZZ[i+1],24);
          double zm2=(i>=2)?ZZ[i-2]:Zval(gram_delta(m-2)), zp2=(i+2<=(int)(hi+1-(lo-1)))?ZZ[i+2]:Zval(gram_delta(m+2));
          printf("CAND %s n=%lld t=%.6f mu=%.6f norm=%.5f Z=%.6f Zp=%.6f Z[-2..2]=%.4f,%.4f,%.4f,%.4f,%.4f zeros_left=%d zeros_right=%d\n",
            tag,m,B.tc+DL[i],mu,mn,Zm,zp,zm2,ZZ[i-1],Zm,ZZ[i+1],zp2,zl,zr); fflush(stdout);
        }
      }
    }
  }
  R_last=hi;
}
static int cmpd(const void *a, const void *b){ double x=*(const double*)a, y=*(const double*)b; return x<y?-1:x>y; }
static void reso_scan(double T, double H, double h, int P, double SEP, double Q, int M, double THR, const char *tag){
  clock_t c0=clock();
  int np=0; double L[64], w[64]; dd Ld[64];
  for(int p=2;p<=P;p++){ int pr=1; for(int q=2;q*q<=p;q++) if(p%q==0){pr=0;break;} if(!pr) continue;
    Ld[np]=dd_mul(dd_log((double)p),INV2PI); L[np]=Ld[np].hi; w[np]=1/sqrt((double)p); np++;
    if(p<=7){ Ld[np]=dd_muld(Ld[np-1],2.0); L[np]=Ld[np].hi; w[np]=0.5/p; np++; } }
  double ph0[64], zr[64], zi[64], wr[64], wi[64];
  for(int i=0;i<np;i++){ dd q=two_prod(T,Ld[i].hi); double f=(q.hi-rint(q.hi))+(q.lo+T*Ld[i].lo); ph0[i]=f-rint(f); wr[i]=cos(TWOPI*h*L[i]); wi[i]=-sin(TWOPI*h*L[i]); }
  long long nsteps=(long long)(H/h), spb=(long long)llround(SEP/h);
  for(int j=0;j<NLOW;j++){ R_lowmu[j]=1e99; R_lownm[j]=1e99; }
  for(int pass=0;pass<2;pass++){
    long long nst = pass==0 ? (nsteps<200000LL*spb?nsteps:200000LL*spb) : nsteps;
    double *pil = pass==0 ? malloc(sizeof(double)*(nst/spb+2)) : NULL; long long npil=0;
    static double FTHR=0;
    double bmax=-1e9, bt=0;
    for(long long j=0;j<nst;j++){
      if((j&4095)==0){ for(int i=0;i<np;i++){ double f=ph0[i]+((double)j*h)*L[i]; f-=rint(f); zr[i]=cos(TWOPI*f); zi[i]=-sin(TWOPI*f); } }
      double F=0; for(int i=0;i<np;i++){ F+=w[i]*zr[i]; double a=zr[i]*wr[i]-zi[i]*wi[i]; zi[i]=zr[i]*wi[i]+zi[i]*wr[i]; zr[i]=a; }
      if(F>bmax){ bmax=F; bt=T+(double)j*h; }
      if((j+1)%spb==0 || j==nst-1){
        if(pass==0) pil[npil++]=bmax; else { R_bins++; if(bmax>FTHR) reso_check(bt,M,THR,tag); }
        bmax=-1e9;
        if(pass==1 && (R_bins%5000000)==0){ printf("PROG %s t=%.1f bins=%lld checked=%lld gp=%lld iso=%lld minmu=%.5f@%lld minnorm=%.5f@%lld secs=%.0f\n",tag,bt,R_bins,R_checked,R_gp,R_iso,R_minmu,R_am,R_minmn,R_an,(double)(clock()-c0)/CLOCKS_PER_SEC); fflush(stdout); }
      }
    }
    if(pass==0){ qsort(pil,npil,sizeof(double),cmpd); FTHR=pil[(long long)(Q*(npil-1))]; free(pil);
      printf("PILOT %s bins=%lld quantile %.4f -> FTHR=%.4f (Fmax possible %.4f)\n",tag,npil,Q,FTHR,0.0); fflush(stdout); }
  }
  printf("RESO %s T=%.3f H=%.1f h=%.3f P=%d SEP=%.2f Q=%.4f M=%d bins=%lld checked=%lld gram_points=%lld bad=%lld iso=%lld minmu=%.6f n=%lld minnorm=%.5f n=%lld maxabsZ=%.4f n=%lld secs=%.1f\n",
    tag,T,H,h,P,SEP,Q,M,R_bins,R_checked,R_gp,R_bad,R_iso,R_minmu,R_am,R_minmn,R_an,R_maxz,R_amz,(double)(clock()-c0)/CLOCKS_PER_SEC);
  printf("HIST %s mu<[4,4.25,4.5,4.75,5,5.5,6,7,8]:",tag); for(int j=0;j<9;j++) printf(" %lld",R_ca[j]);
  printf(" norm<[.35,.40,.45,.50,.55,.60,.65,.70,.80]:"); for(int j=0;j<9;j++) printf(" %lld",R_cn[j]); printf("\n");
  printf("LOW %s mu:",tag); for(int j=0;j<NLOW;j++) printf(" %.5f",R_lowmu[j]); printf(" | norm:"); for(int j=0;j<NLOW;j++) printf(" %.5f",R_lownm[j]); printf("\n");
  fflush(stdout);
}
int main(int argc,char**argv){
  init_tables(300000);
  if(argc>2 && !strcmp(argv[1],"-check")){
    for(int a=2;a<argc;a++){ long long n=atoll(argv[a]); double tc=gram_est(n); setup_block(tc, tc+10);
      double d=gram_delta(n); double Z=Zval(d), Zp=Zder(d); double to=gram_est(n); double Zo=Zval_o(to), Zpo=Zder_o(to);
      printf("CHECK n=%lld tc=%.17g d=%.17g Z=%.15e Zp=%.15e old_t=%.17g old_Z=%.15e old_Zp=%.15e\n",n,tc,d,Z,Zp,to,Zo,Zpo); fflush(stdout); }
    return 0;
  }
  if(argc>11 && !strcmp(argv[1],"-reso")){   /* -reso T H h P SEP Q M THR NTHR TAG */
    NTHR=atof(argv[10]); reso_scan(atof(argv[2]),atof(argv[3]),atof(argv[4]),atoi(argv[5]),atof(argv[6]),atof(argv[7]),atoi(argv[8]),atof(argv[9]),argv[11]);
    return 0;
  }
  FILE *f=fopen(argv[1],"r"); double THR=atof(argv[2]); if(argc>3) NTHR=atof(argv[3]); if(argc>4) ZG=atoi(argv[4]); char id[64]; long long n0,cnt;
  while(fscanf(f,"%63s %lld %lld",id,&n0,&cnt)==3) window(id,n0,cnt,THR);
  return 0;
}
