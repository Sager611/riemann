// Fast screening scan (double precision, Riemann-Siegel main sum + C0) for isolated bad Gram points with small
// viscosity mu = |Z'(g_n)/Z(g_n)| (Jerby 2406.06548 conj:rep, 2310.14415 Conj G-B-G: mu > 4). Candidates are certified in Arb.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
static const double PI=3.14159265358979323846, TWOPI=6.28318530717958647692;
static double theta(double t){ return t/2*log(t/TWOPI)-t/2-PI/8+1/(48*t)+7/(5760*t*t*t)+31/(80640*pow(t,5)); }
static double dtheta(double t){ return 0.5*log(t/TWOPI)-1/(48*t*t)-7/(1920*t*t*t*t); }
static double C0f(double p){ return cos(TWOPI*(p*p-p-1.0/16))/cos(TWOPI*p); }
static double C0(double p){ double c=cos(TWOPI*p); if(fabs(c)<1e-3){ double h=2e-3; return 0.5*(C0f(p-h)+C0f(p+h)) - (C0f(p-2*h)+C0f(p+2*h)-C0f(p-h)-C0f(p+h))/6.0; } return C0f(p); }
static double C[15];
#define COS2PI(y,z) (C[0]+z*(C[1]+z*(C[2]+z*(C[3]+z*(C[4]+z*(C[5]+z*(C[6]+z*(C[7]+z*C[8]))))))))   /* |y|<=1/4: error < 1e-10 */
static double *L2, *rs, *lg;
static double Zval(double t){
  double a=sqrt(t/TWOPI); int N=(int)a; double p=a-N; double th2=theta(t)/TWOPI;
  double s=0;
  #pragma clang loop vectorize(enable)
  for(int k=1;k<=N;k++){ double y=th2-t*L2[k]; y-=rint(y); double ay=fabs(y); double sg=ay>0.25?-1.0:1.0; double yr=ay>0.25?0.5-ay:ay; double z=yr*yr; s+=rs[k]*sg*COS2PI(yr,z); }
  double R=((N-1)&1?-1.0:1.0)*pow(a,-0.5)*C0(p);
  return 2*s+R;
}
static double Zder(double t){
  double a=sqrt(t/TWOPI); int N=(int)a; double th=theta(t), thp=dtheta(t); double sp=0;
  for(int k=1;k<=N;k++){ double x=th-t*lg[k]; x-=TWOPI*nearbyint(x/TWOPI); sp+=rs[k]*(thp-lg[k])*sin(x); }
  return -2*sp;
}
int main(int argc,char**argv){
  long n0=atol(argv[1]), n1=atol(argv[2]); double THR=atof(argv[3]);
  double f=1; for(int k=0;k<15;k++){ C[k]=((k&1)?-1:1)*pow(TWOPI,2*k)/f; f*=(2*k+1)*(2*k+2); }
  int KMAX=40000; L2=malloc(sizeof(double)*(KMAX+1)); rs=malloc(sizeof(double)*(KMAX+1)); lg=malloc(sizeof(double)*(KMAX+1));
  for(int k=1;k<=KMAX;k++){lg[k]=log((double)k); L2[k]=lg[k]/TWOPI; rs[k]=1/sqrt((double)k);}
  double x=(n0-1+0.125)/exp(1.0), w=log(x+1); for(int i=0;i<60;i++) w-= (w*exp(w)-x)/(exp(w)*(w+1));
  double t=TWOPI*exp(1.0)*exp(w);
  double tz[3]={0}, zz[3]={0}; int cnt=0;
  long nbad=0,niso=0; double minm=1e99, minmn=1e99; long argm=-1, argmn=-1; long ncand=0;
  for(long n=n0-1;n<=n1;n++){
    for(int i=0;i<3;i++) t-= (theta(t)-PI*(double)n)/dtheta(t);
    double Z=Zval(t);
    tz[0]=tz[1];zz[0]=zz[1]; tz[1]=tz[2];zz[1]=zz[2]; tz[2]=t; zz[2]=Z; cnt++;
    long m=n-1;
    if(cnt>=3 && m>=n0 && m<n1){
      double sg=(m&1)?-1.0:1.0;
      if(sg*zz[1]<0){ nbad++;
        if((-sg)*zz[0]>0 && (-sg)*zz[2]>0){ niso++; double zp=Zder(tz[1]); double mu=fabs(zp/zz[1]); double mn=mu/dtheta(tz[1]);
          if(mu<minm){minm=mu;argm=m;} if(mn<minmn){minmn=mn;argmn=m;}
          if(mu<THR){ ncand++; printf("CAND n=%ld t=%.9f Z=%.6e Zp=%.6e mu=%.6f mu_norm=%.5f\n",m,tz[1],zz[1],zp,mu,mn); fflush(stdout);} } } }
    t += PI/dtheta(t);
    if(n%1000000==0){ printf("progress n=%ld t=%.1f bad=%ld iso=%ld min_mu=%.5f@%ld min_norm=%.5f@%ld\n",n,t,nbad,niso,minm,argm,minmn,argmn); fflush(stdout);}
  }
  printf("DONE [%ld,%ld) bad=%ld iso=%ld min_mu=%.6f@%ld min_norm=%.6f@%ld cands=%ld\n",n0,n1,nbad,niso,minm,argm,minmn,argmn,ncand);
  return 0;
}
