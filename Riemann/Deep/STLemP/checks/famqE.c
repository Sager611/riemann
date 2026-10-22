// famq.c: family sums over all imaginary quadratic orders, X/4 <= |Delta| <= 2X, weight Phi(|Delta|/X):
//   Qnew = sum_Delta Phi * Q'_Delta(b)  (Lemma P: no coprimality to the conductor; = sum over e | f of blocks)
//   Qold = sum_Delta Phi * sum_C |S_C|^2 (ideals prime to f, as in RCApps3)
//   Qfund = fundamental Delta only.  b supported on odd squarefree n prime to 3.
// usage: famqE X N mode sigma  (modes 0-2 as famq.c; 3: lambda_E(n) n^-s phi2 (E = 37a, m=1, no mollifier); 4: Sym^2 lambda (m=2); 5: lambda_E with mollifier a_V (m=1 detector, V=sqrt(X)/2); 6: Sym^2 with mollifier V=X^{1/4})
// STLemP check, adapted from RCApps5/checks/famq.c
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>

static double *bb; static char *adm; static int NMAX;
static double sigma=5.0/6.0;

static long gcdl(long a,long b){ if(a<0)a=-a; if(b<0)b=-b; while(b){long t=a%b;a=b;b=t;} return a; }
static int isdisc(long D){ long r=((D%4)+4)%4; return r==0||r==1; }
static long conductor(long D){ long f=1,Dc=D; for(long p=2;p*p<=labs(Dc);p++){ while(Dc%(p*p)==0 && isdisc(Dc/(p*p))){Dc/=p*p; f*=p;} } return f; }
static double phi2(double v){ if(v<=0.25||v>=4.0) return 0.0; double t=log(v)/log(4.0); return exp(1.0-1.0/(1.0-t*t)); }
static double PhiW(double u){ if(u<=0.25||u>=2.0) return 0.0; if(u<0.5){double s=sin(M_PI/2*(u-0.25)/0.25); return s*s;} if(u<=1.0) return 1.0; double c=cos(M_PI/2*(u-1.0)); return c*c; }

int main(int argc,char**argv){
  alarm(1500);
  long X=atol(argv[1]); long N=atol(argv[2]); int mode=atoi(argv[3]); if(argc>4) sigma=atof(argv[4]);
  NMAX=4*N+1;
  bb=calloc(NMAX+1,sizeof(double)); adm=calloc(NMAX+1,1);
  // squarefree, odd, prime to 3
  char *sqf=malloc(NMAX+1); memset(sqf,1,NMAX+1);
  for(long p=2;p*p<=NMAX;p++) for(long m=p*p;m<=NMAX;m+=p*p) sqf[m]=0;
  for(long n=1;n<=NMAX;n++) adm[n]= sqf[n] && (n%2==1) && (n%3!=0) && (n%37!=0);
  // a_p for E = 37a: Y^2 = 4x^3 - 4x + 1, a_p = -sum_x legendre(4x^3-4x+1, p)
  double *lam=calloc(NMAX+1,sizeof(double)); for(long n=0;n<=NMAX;n++) lam[n]=1.0;
  if(mode>=3){
    char *isp=calloc(NMAX+1,1); for(long n=2;n<=NMAX;n++) isp[n]=1; for(long p=2;p*p<=NMAX;p++) if(isp[p]) for(long m=p*p;m<=NMAX;m+=p) isp[m]=0;
    signed char *qr=malloc(NMAX+1);
    for(long p=5;p<=NMAX;p++){ if(!isp[p]||p==37) continue; memset(qr,-1,p); qr[0]=0; for(long x=1;x<p;x++) qr[(x*x)%p]=1;
      long s=0; for(long x=0;x<p;x++){ long v=((4*((x*x)%p)%p*x)%p - (4*x)%p + 1 + 2*p)%p; s+=qr[v]; }
      double lp=-(double)s/sqrt((double)p); if(mode==4||mode==6) lp=lp*lp-1.0;
      for(long m=p;m<=NMAX;m+=p) lam[m]*=lp; }
    free(qr); free(isp);
  }
  // mobius up to V
  long V=(long)(sqrt((double)X)/2); if(mode==6) V=(long)pow((double)X,0.25); if(V<2)V=2;
  int *mu=malloc((V+1)*sizeof(int)); for(long k=0;k<=V;k++) mu[k]=1;
  { char *comp=calloc(V+1,1); for(long p=2;p<=V;p++) if(!comp[p]){ for(long m=p;m<=V;m+=p){ if(m>p) comp[m]=1; mu[m]=-mu[m]; } for(long m=p*p;m<=V;m+=p*p) mu[m]=0; } free(comp); }
  double *aV=calloc(NMAX+1,sizeof(double));
  for(long k=1;k<=V;k++) if(mu[k]) for(long m=k;m<=NMAX;m+=k) aV[m]+=mu[k];
  srand48(12345+mode);
  double l2=pow(4.0,1.0-sigma);
  double nb=0;
  for(long n=1;n<=NMAX;n++){ if(!adm[n]) continue; double v=(double)n/N; double val=0;
    if(mode==0) val=aV[n]*pow(n,-sigma)*(phi2(v)-l2*phi2(4*v));
    else if(mode==1) val=(drand48()<0.5?-1:1)*pow(n,-sigma)*phi2(v);
    else if(mode==2) val=pow(n,-sigma)*phi2(v);
    else if(mode==3||mode==4) val=lam[n]*pow(n,-sigma)*phi2(v);
    else val=lam[n]*aV[n]*pow(n,-sigma)*phi2(v);
    bb[n]=val; nb+=val*val; }
  double Qnew=0,Qold=0,Qfund=0, famfund=0; long nD=0, nfund=0;
  int maxforms=200000; long *FA=malloc(maxforms*sizeof(long)),*FB=malloc(maxforms*sizeof(long)),*FC=malloc(maxforms*sizeof(long));
  double *S=malloc(maxforms*sizeof(double));
  for(long Dabs=X/4; Dabs<=2*X; Dabs++){
    long D=-Dabs; if(!isdisc(D)) continue; double w=PhiW((double)Dabs/X); if(w==0) continue;
    long f=conductor(D);
    // odd squarefree divisors e of f, prime to 3 (b vanishes on multiples of 3)
    long fp[32]; int nfp=0; { long ff=f; for(long p=2;p<=ff;p++){ if(ff%p==0){ if(p%2==1 && p%3!=0) fp[nfp++]=p; while(ff%p==0) ff/=p; } } }
    double Qp=0, Q1=0; long hD=0;
    for(int mask=0; mask<(1<<nfp); mask++){
      long e=1; for(int i=0;i<nfp;i++) if(mask>>i&1) e*=fp[i];
      long Dp=D/(e*e); long Dpa=-Dp; long mmax=(4*N)/e; if(mmax<1) continue;
      int wp = (Dp==-3)?6:((Dp==-4)?4:2);
      // reduced primitive forms of Dp
      int nf=0;
      for(long a=1; 3*a*a<=Dpa; a++) for(long b=-a+1;b<=a;b++){ long num=b*b+Dpa; if(num%(4*a)) continue; long c=num/(4*a); if(c<a) continue; if(b<0 && a==c) continue; if(gcdl(gcdl(a,b),c)!=1) continue; FA[nf]=a;FB[nf]=b;FC[nf]=c;nf++; }
      if(e==1) hD=nf;
      double blk=0;
      for(int i=0;i<nf;i++){ long a=FA[i],b=FB[i],c=FC[i]; double s=0;
        long ymax=(long)sqrt(4.0*a*mmax/(double)Dpa)+1;
        for(long y=-ymax;y<=ymax;y++){ double disc=(double)b*b*y*y-4.0*a*((double)c*y*y-mmax); if(disc<0) continue; double sq=sqrt(disc);
          long x0=(long)floor((-b*y-sq)/(2.0*a))-1, x1=(long)ceil((-b*y+sq)/(2.0*a))+1;
          for(long x=x0;x<=x1;x++){ long m=a*x*x+b*x*y+c*y*y; if(m<1||m>mmax) continue; long n=e*m; if(n>NMAX||!adm[n]) continue; if(gcdl(m,f)!=1) continue; s+=bb[n]; } }
        s/=wp; blk+=s*s; }
      blk*=wp/2.0;
      Qp+=blk; if(e==1) Q1=blk;
    }
    Qnew+=w*Qp; Qold+=w*Q1; nD++;
    if(f==1){ Qfund+=w*Q1; nfund++; famfund+=w*hD*Q1; }
  }
  printf("sigma=%.3f X=%ld N=%ld mode=%d V=%ld nD=%ld nfund=%ld  ||b||^2=%.6e\n",sigma,X,N,mode,V,nD,nfund,nb);
  printf("Qnew/(X||b||^2)=%.5f  Qold/(X||b||^2)=%.5f  Qfund/(X||b||^2)=%.5f  extra=(Qnew-Qold)/Qold=%.4f  fam_fund/(X^1.5||b||^2)=%.5f\n",
    Qnew/(X*nb),Qold/(X*nb),Qfund/(X*nb),(Qnew-Qold)/Qold, famfund/(pow(X,1.5)*nb));
  return 0;
}

