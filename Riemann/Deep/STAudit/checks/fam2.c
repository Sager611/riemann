// fam2.c (staudit, independent of STLemP famqE.c / RC5Audit fam.c)
// F = sum_{D fund <0, X/2<|D|<=X} h(D) sum_{classes C} |S_C(b)|^2,  S_C(b) = sum_n b(n) r_C(n)
//   (r_C(n) = R_Q(n)/2 via lattice points of the reduced form Q of C), 
// G = sum_D h(D) sum_n |b(n)|^2 r_D(n)  (pure diagonal: same ideal), R1 = F/G.
// b(n) = coef(n) n^{-1/2} phi(n/N), n squarefree with all prime factors >= 19, phi bump on [1,2].
// modes: 0 lambda_E (11a1), 1 lambda_Sym2E, 2 coef = 1 (pole present).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define NM 3
static int *spf; static double *lam; static int maxn;
static double *bb[NM];
static long gcdl(long a, long b){ if(a<0)a=-a; if(b<0)b=-b; while(b){long t=a%b;a=b;b=t;} return a; }
static int sqfree_odd(long m){ for(long p=2;p*p<=m;p++){ if(m%(p*p)==0) return 0; } return 1; }
static int isfund(long D){ // D<0
  long m=-D; if((D%4+4)%4==1) return sqfree_odd(m);
  if(m%4==0){ long k=m/4; long kk=(-k%4+4)%4; if((kk==2||kk==3) && sqfree_odd(k)) return 1; } return 0; }
int main(int argc,char**argv){
  long X=atol(argv[1]); int nN=argc-2; double Nf[16]; for(int i=0;i<nN;i++) Nf[i]=atof(argv[2+i]);
  long Nmax=0; for(int i=0;i<nN;i++){ long N=(long)(Nf[i]*X); if(N>Nmax)Nmax=N; }
  maxn=2*Nmax+2; spf=calloc(maxn+1,sizeof(int)); lam=calloc(maxn+1,sizeof(double));
  for(long i=2;i<=maxn;i++) if(!spf[i]) for(long j=i;j<=maxn;j+=i) if(!spf[j]) spf[j]=i;
  FILE*f=fopen("ap.txt","r"); long p,a; while(fscanf(f,"%ld %ld",&p,&a)==2){ if(p<=maxn) lam[p]=a/sqrt((double)p); } fclose(f);
  for(int k=0;k<NM;k++) bb[k]=calloc(maxn+1,sizeof(double));
  printf("# X=%ld\n",X);
  for(int iN=0;iN<nN;iN++){
    long N=(long)(Nf[iN]*X);
    for(int k=0;k<NM;k++) memset(bb[k],0,(maxn+1)*sizeof(double));
    double nb2[NM]={0,0,0};
    for(long n=N+1;n<2*N;n++){
      long m=n; double c0=1,c1=1; int ok=1;
      while(m>1){ long q=spf[m]; m/=q; if(m%q==0||q<19){ok=0;break;} c0*=lam[q]; c1*=lam[q]*lam[q]-1; }
      if(!ok) continue;
      double v=(double)n/N, u=2*v-3, ph=exp(1-1/(1-u*u)), w=ph/sqrt((double)n);
      bb[0][n]=c0*w; bb[1][n]=c1*w; bb[2][n]=w;
      nb2[0]+=bb[0][n]*bb[0][n]; nb2[1]+=bb[1][n]*bb[1][n]; nb2[2]+=bb[2][n]*bb[2][n];
    }
    double F[NM]={0,0,0}, G[NM]={0,0,0}, H=0, Hb[NM]={0,0,0}, maxratio[NM]={0,0,0}; long nD=0;
    for(long m=X/2+1;m<=X;m++){
      long D=-m; if(!isfund(D)) continue;
      // reduced forms
      long am=(long)sqrt(m/3.0)+1; int h=0; double FD[NM]={0,0,0}, GD[NM]={0,0,0};
      for(long A=1;A<=am;A++){
        for(long B=-A+1;B<=A;B++){
          if(((B-D)%2+2)%2) continue;
          long num=B*B-D; if(num%(4*A)) continue; long C=num/(4*A);
          if(C<A) continue; if(C==A && B<0) continue; if(gcdl(gcdl(A,B),C)!=1) continue;
          h++;
          double S[NM]={0,0,0}, T[NM]={0,0,0};
          double ymax=sqrt(8.0*A*N/(double)m)+1;
          for(long y=-(long)ymax;y<=(long)ymax;y++){
            double Dy=(double)m*y*y/(4.0*A); if(Dy>=2.0*N) continue;
            double tlo=(N-Dy)>0? sqrt((N-Dy)/A):0, thi=sqrt((2.0*N-Dy)/A), sh=(double)B*y/(2.0*A);
            for(int sg=-1;sg<=1;sg+=2){
              double x0 = sg>0? sh*(-1)+tlo : sh*(-1)-thi, x1 = sg>0? -sh+thi : -sh-tlo;
              for(long x=(long)floor(x0)-1;x<=(long)ceil(x1)+1;x++){
                long Q=A*x*x+B*x*y+C*y*y; if(Q<=N||Q>=2*N) continue;
                // avoid double count when tlo==0 intervals overlap: handle via sign of t
                double t=x+sh; if(sg>0 && t<0) continue; if(sg<0 && t>=0) continue;
                for(int k=0;k<NM;k++){ double v=bb[k][Q]; S[k]+=v; T[k]+=v*v; }
              }
            }
          }
          for(int k=0;k<NM;k++){ FD[k]+=0.25*S[k]*S[k]; GD[k]+=0.5*T[k]; }
        }
      }
      nD++; H+=h;
      for(int k=0;k<NM;k++){ F[k]+=h*FD[k]; G[k]+=h*GD[k]; Hb[k]+=h*nb2[k]; if(GD[k]>0 && FD[k]/GD[k]>maxratio[k]) maxratio[k]=FD[k]/GD[k]; }
    }
    printf("N=%ld N/X=%.4f nD=%ld sum_h=%.0f",N,(double)N/X,nD,H);
    for(int k=0;k<NM;k++) printf(" | m%d R1=F/G=%.4f F/(sum_h*|b|^2)=%.4f maxD(F_D/G_D)=%.2f",k,F[k]/G[k],F[k]/Hb[k],maxratio[k]);
    printf("\n"); fflush(stdout);
  }
  return 0;
}
