/* floor.c — float64 exploration of the Polymath/Gomila Euler-mollified
   triangle floor at one Riemann–Siegel window N (first window is binding).
   floor = (1 - g - S)/modmoll - corr,  S = sum_{n>=2} W(b_n,a_n) n^{-sigma}
   b_n = sum_{d|n, d|D, n/d<=N} lam_d bt(n/d),  a_n = same with (n/d)^y
   W = |b|+g|a| (T)  or  max(|b-g a|, c|b+g a|), c=(1-g)/(1+g)  (S)
   usage: floor N t y2 T|S p1 p2 ... [@lamscale]   (lam_p = s*bt(p), s=1 default)
   Exploration only (not a certificate). */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef unsigned long long u64;
static int K; static u64 P[16]; static double lamp[16];
static double *lg, *bt, *ky; static u64 NN; static double sig, g, c; static char wt;
static double S; static u64 cnt;
static double lamd[1<<16]; static u64 dval[1<<16];
static void visit(u64 r, double lr, u64 s, double ls, int radmask){
  u64 n = s*r; if(n==1) return;
  double b=0,a=0;
  /* iterate submasks of radmask */
  int m = radmask;
  for(;;){
    u64 d = dval[m];
    if(n % d == 0){ u64 k = n/d; if(k<=NN){ double w = lamd[m]*bt[k]; b += w; a += w*ky[k]; } }
    if(m==0) break; m = (m-1)&radmask;
  }
  double W;
  if(wt=='T') W = fabs(b)+g*fabs(a);
  else { double u=fabs(b-g*a), v=c*fabs(b+g*a); W = u>v?u:v; }
  S += W*exp(-sig*(lr+ls)); cnt++;
}
/* DFS over E-smooth s with prod p^{max(e-1,0)} <= lim */
static void dfs(int i, u64 r, double lr, u64 s, double ls, int radmask, u64 sp, u64 lim){
  if(i==K){ visit(r,lr,s,ls,radmask); return; }
  /* exponent 0 */
  dfs(i+1,r,lr,s,ls,radmask,sp,lim);
  /* exponent >=1 */
  u64 s2 = s*P[i]; double ls2 = ls+log((double)P[i]); u64 sp2 = sp; 
  for(;;){
    dfs(i+1,r,lr,s2,ls2,radmask|(1<<i),sp2,lim);
    if(sp2 > lim / P[i]) break;
    sp2 *= P[i]; s2 *= P[i]; ls2 += log((double)P[i]);
    if(sp2 > lim) break;
  }
}
int main(int argc,char**argv){
  if(argc<6){fprintf(stderr,"usage\n");return 2;}
  NN = strtoull(argv[1],0,10); double t=atof(argv[2]), y2=atof(argv[3]); wt=argv[4][0];
  double y=sqrt(y2); double lamscale=1.0; K=0;
  for(int i=5;i<argc;i++){ if(argv[i][0]=='@'){lamscale=atof(argv[i]+1);continue;} P[K++]=strtoull(argv[i],0,10);} 
  lg=malloc((NN+1)*sizeof(double)); bt=malloc((NN+1)*sizeof(double)); ky=malloc((NN+1)*sizeof(double));
  for(u64 k=1;k<=NN;k++){ lg[k]=log((double)k); bt[k]=exp(0.25*t*lg[k]*lg[k]); ky[k]=exp(y*lg[k]); }
  double pi=M_PI; double xN=4*pi*(double)NN*(double)NN - pi*t/4; double L=log(xN/(4*pi));
  double inner = 1-3*y+4*y*(1+y)/(xN*xN); if(inner<0) inner=0;
  sig = (1+y)/2 + t/4*L - t/(2*xN*xN)*inner;
  g = exp(0.02*y - y/2*L); c=(1-g)/(1+g);
  double modmoll=1;
  for(int i=0;i<K;i++){ double l=log((double)P[i]); lamp[i]=lamscale*exp(0.25*t*l*l); modmoll*=1+lamp[i]*exp(-sig*l); }
  for(int m=0;m<(1<<K);m++){ u64 d=1; double lam=1; for(int i=0;i<K;i++) if(m&(1<<i)){ d*=P[i]; lam*=-lamp[i]; } dval[m]=d; lamd[m]=lam; }
  S=0; cnt=0;
  for(u64 r=1;r<=NN;r++){
    int ok=1; for(int i=0;i<K;i++) if(r%P[i]==0){ok=0;break;}
    if(!ok) continue;
    dfs(0,r,lg[r],1,0.0,0,1,NN/r);
  }
  /* kappa correction (negligible): corr = g*sum_{n=2}^N bt(n)(n^{K}-1) n^{-(sig-y)} */
  double modK = t*y/(2*(xN-6)); double corr=0;
  for(u64 n=2;n<=NN;n++) corr += bt[n]*expm1(modK*lg[n])*exp(-(sig-y)*lg[n]);
  corr*=g;
  double fl = (1-g-S)/modmoll - corr;
  printf("N=%llu t=%.9f y2=%.9f y=%.7f W=%c primes=%d lamscale=%.4f sigma=%.9f g=%.9f modmoll=%.9f S=%.12f corr=%.3e floor=%.12e terms=%llu Lambda=%.9f\n",
    NN,t,y2,y,wt,K,lamscale,sig,g,modmoll,S,corr,fl,cnt,t+y2/2);
  return 0;
}
