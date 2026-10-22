// famq_all.c (AllDisc): copy of RCApps5/checks/famq.c that also records the LEFT side of Theorem 1^o over ALL orders:
//   fam_all = sum_Delta Phi(|Delta|/X) h(Delta) sum_C |S_C|^2  (S_C over ideals prime to f; every order, not only fundamental)
//   fam_nonfund = the same over non-fundamental Delta only; maxh = max h(Delta) in the range. Theorem 1^o: fam_all <= maxh * Qnew.
// Original header: family sums over all imaginary quadratic orders, X/4 <= |Delta| <= 2X, weight Phi(|Delta|/X):
//   Qnew = sum_Delta Phi * Q'_Delta(b)  (Lemma P: no coprimality to the conductor; = sum over e | f of blocks)
//   Qold = sum_Delta Phi * sum_C |S_C|^2 (ideals prime to f, as in RCApps3)
//   Qfund = fundamental Delta only.  b supported on odd squarefree n prime to 3.
// usage: famq X N mode   (mode 0: detector b with pole-killing weight; 1: random signs; 2: unmollified n^-s phi2, pole present)
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
  long X=atol(argv[1]); long N=atol(argv[2]); int mode=atoi(argv[3]);
  NMAX=4*N+1;
  bb=calloc(NMAX+1,sizeof(double)); adm=calloc(NMAX+1,1);
  // squarefree, odd, prime to 3
  char *sqf=malloc(NMAX+1); memset(sqf,1,NMAX+1);
  for(long p=2;p*p<=NMAX;p++) for(long m=p*p;m<=NMAX;m+=p*p) sqf[m]=0;
  for(long n=1;n<=NMAX;n++) adm[n]= sqf[n] && (n%2==1) && (n%3!=0);
  // mobius up to V
  long V=(long)(sqrt((double)X)/2); if(V<2)V=2;
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
    else val=pow(n,-sigma)*phi2(v);
    bb[n]=val; nb+=val*val; }
  double Qnew=0,Qold=0,Qfund=0, famfund=0, famall=0, famnon=0; long nD=0, nfund=0, maxh=0;
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
    Qnew+=w*Qp; Qold+=w*Q1; nD++; famall+=w*hD*Q1; if(f>1) famnon+=w*hD*Q1; if(hD>maxh) maxh=hD;
    if(f==1){ Qfund+=w*Q1; nfund++; famfund+=w*hD*Q1; }
  }
  printf("X=%ld N=%ld mode=%d V=%ld nD=%ld nfund=%ld  ||b||^2=%.6e\n",X,N,mode,V,nD,nfund,nb);
  printf("Qnew/(X||b||^2)=%.5f  Qold/(X||b||^2)=%.5f  Qfund/(X||b||^2)=%.5f  extra=(Qnew-Qold)/Qold=%.4f  fam_fund/(X^1.5||b||^2)=%.5f\n",
    Qnew/(X*nb),Qold/(X*nb),Qfund/(X*nb),(Qnew-Qold)/Qold, famfund/(pow(X,1.5)*nb));
  printf("ALL ORDERS: fam_all/(X^1.5||b||^2)=%.5f  fam_nonfund/(X^1.5||b||^2)=%.5f  fam_fund/(X^1.5||b||^2)=%.5f  maxh=%ld  fam_all/(maxh*Qnew)=%.4f (Thm 1^o needs <= 1)\n",
    famall/(pow(X,1.5)*nb), famnon/(pow(X,1.5)*nb), famfund/(pow(X,1.5)*nb), maxh, famall/((double)maxh*Qnew));
  return 0;
}

