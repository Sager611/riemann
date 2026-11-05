/* nflp numerics: L^p(0,1) norms of the Riesz-Moebius approximants for Dedekind zeta functions.
   F_N = 1 + sum_{k<=N} b_K(k)(1-k/N) lambda_k,  lambda_k(x) = kappa_K(1/(kx)) - kappa_K(1/x)/k,
   kappa_K(y) = c_K y - A_K(y),  A_K(y) = #{ideals a : N(a) <= y},  b_K(k) = sum_{N(a)=k} mu_K(a).
   On I_m = (1/(m+1), 1/m]:  lambda_k = A_K(m)/k - A_K(floor(m/k)), hence
   F_N(I_m) = v(m) = 1 + P_N(1) A_K(m) - sum_{n<=m} g_N(n),  g_N(n) = sum_{k|n, k<=N} b_K(k)(1-k/N) a_K(n/k),
   P_N(1) = sum_{k<=N} b_K(k)(1-k/N)/k.  ||F_N||_p^p = sum_m |v(m)|^p/(m(m+1)).
   Fields (type): 0 = Q (regression against NBLp), 1 = Q(i) (chi_{-4}), 2 = Q(sqrt 5) (chi_5).
   For quadratic K: a_K = 1 * chi, b_K multiplicative with b(p) = -(1+chi(p)), b(p^2) = chi(p), b(p^e) = 0 (e >= 3).
   N*g_N and N*sum g_N are exact 64-bit integers; v in long double.
   usage: nbk M type N1 [N2 ...]   (single-threaded; memory 2M bytes + 32 MB) */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define KMAX 100000
static int bK[KMAX+1];
static int chi(int type, long long n){ int r;
  if(type==0) return 1;
  if(type==1){ r=(int)(n%4); return r==1?1:(r==3?-1:0); }
  r=(int)(n%5); return (r==1||r==4)?1:((r==2||r==3)?-1:0); }
static void make_b(int type){ static int spf[KMAX+1]; for(int i=0;i<=KMAX;i++) spf[i]=0;
  for(int i=2;i<=KMAX;i++) if(!spf[i]) for(int j=i;j<=KMAX;j+=i) if(!spf[j]) spf[j]=i;
  bK[1]=1;
  for(int k=2;k<=KMAX;k++){ int n=k, val=1;
    while(n>1){ int p=spf[n], e=0; while(n%p==0){ n/=p; e++; }
      int f; if(type==0) f=(e==1)?-1:0; else { int c=chi(type,p); f=(e==1)?-(1+c):((e==2)?c:0); }
      val*=f; if(!val) break; }
    bK[k]=val; } }
int main(int argc,char**argv){
  long long M=atoll(argv[1]); int type=atoi(argv[2]); int nN=argc-3; make_b(type);
  short *a=malloc((M+1)*sizeof(short)); if(!a){fprintf(stderr,"oom\n");return 1;}
  /* a_K(j) = sum_{d|j} chi(d), segmented divisor-pair sieve */
  const long long B=1LL<<22;
  if(type==0){ for(long long j=1;j<=M;j++) a[j]=1; }
  else { a[0]=0; for(long long st=1;st<=M;st+=B){ long long en=st+B-1; if(en>M) en=M;
      for(long long j=st;j<=en;j++) a[j]=0;
      for(long long d=1; d*d<=en; d++){ int cd=chi(type,d); long long j0=((st+d-1)/d)*d; if(j0<d*d) j0=d*d;
        for(long long j=j0;j<=en;j+=d){ long long e=j/d; if(e==d) a[j]+=cd; else a[j]+=cd+chi(type,e); } } } }
  double cK = (type==0)?1.0:((type==1)?M_PI/4:2*log((1+sqrt(5.0))/2)/sqrt(5.0));
  { long long A=0; double mx=0; long long arg=0; for(long long m=1;m<=M;m++){ A+=a[m]; double D=fabs((double)A-cK*m)/sqrt((double)m); if(D>mx){mx=D;arg=m;} }
    printf("type=%d M=%lld A_K(M)=%lld c_K=%.12f max|Delta(m)|/sqrt(m)=%.6f at m=%lld\n",type,M,A,cK,mx,arg); fflush(stdout); }
  long long *d=malloc(B*sizeof(long long));
  const double ps[3]={32.0/31.0,1.1,1.5};
  for(int ia=0;ia<nN;ia++){ int N=atoi(argv[3+ia]);
    long double NP=0, sabs=0; for(int k=1;k<=N;k++){ NP+=(long double)bK[k]*(N-k)/k; sabs+=(long double)abs(bK[k])*(N-k)/N; }
    long double S[4]={0,0,0,0}; long long NG=0, A=0; long double vmax=0; long long next_ck=1000000;
    for(long long st=1;st<=M;st+=B){ long long en=st+B-1; if(en>M) en=M; long long len=en-st+1;
      memset(d,0,len*sizeof(long long));
      for(int k=1;k<N;k++){ if(!bK[k]) continue; long long w=(long long)bK[k]*(N-k);
        long long j0=(st+k-1)/k, j1=en/k; for(long long j=j0;j<=j1;j++){ long long n=k*j; d[n-st]+=w*a[j]; } }
      for(long long i=0;i<len;i++){ long long m=st+i; NG+=d[i]; A+=a[m];
        long double v=((long double)N+NP*(long double)A-(long double)NG)/N; double av=fabsl(v); if(av>vmax) vmax=av;
        long double wm=1.0L/((long double)m*(m+1));
        for(int q=0;q<3;q++) S[q]+=(long double)pow(av,ps[q])*wm; S[3]+=(long double)av*av*wm;
        if(m==next_ck || m==M){ printf("type=%d N=%d m=%lld S32/31=%.12Lf S1.1=%.12Lf S1.5=%.12Lf S2=%.12Lf | norms %.10Lf %.10Lf %.10Lf %.10Lf max|v|=%.4Lf\n",
              type,N,m,S[0],S[1],S[2],S[3],powl(S[0],31.0L/32.0L),powl(S[1],1/1.1L),powl(S[2],1/1.5L),sqrtl(S[3]),vmax); fflush(stdout);
          if(m==next_ck) next_ck*=10; }
      } }
    printf("SUMMARY type=%d N=%d NP1=%.12Lf sum|c_k|=%.6Lf\n",type,N,NP/N,sabs); fflush(stdout);
  }
  free(d); free(a); return 0; }
