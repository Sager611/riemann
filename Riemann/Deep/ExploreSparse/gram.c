// ceilaudit6: Nyman-Beurling Gram matrices in C (double). Modes:
//  per q N out S0..S_{q-1}: periodic-coefficient Dirichlet series, A(y)=S[floor(y) mod q] (DH: q=5, chi_-4: q=4),
//      same formula as explore38 s7/s9 gdh_coprime: g(h,k)=sum over merged breakpoints c*(psi(b/P)-psi(a/P))/P, P=q h k.
//  zeta N out: Vasyunin V1 Gram (explore38 s1_dN.py formula).
// Writes full N x N row-major float64 to out. G[h][k] = g(h/d,k/d)/d.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
static double psi(double x){ double r=0; while(x<10.0){ r-=1.0/x; x+=1.0; }
  double y2=1.0/(x*x);
  double s=log(x)-0.5/x - y2*(1.0/12 - y2*(1.0/120 - y2*(1.0/252 - y2*(1.0/240 - y2*(1.0/132 - y2*(691.0/32760 - y2*(1.0/12)))))));
  return r+s; }
static long gcdl(long a,long b){ while(b){ long t=a%b; a=b; b=t;} return a; }
static double *gt; static long Nn;
static inline long toff(long k){ return k*(k-1)/2; }  // gt[toff(k)+h-1], h<=k
int main(int argc,char**argv){
  clock_t t0=clock();
  if(argc<4){ fprintf(stderr,"usage\n"); return 1; }
  const char* mode=argv[1];
  if(!strcmp(mode,"per")){
    int q=atoi(argv[2]); long N=atol(argv[3]); const char* out=argv[4]; double S[16];
    for(int i=0;i<q;i++) S[i]=atof(argv[5+i]);
    Nn=N; double **T=malloc((N+1)*sizeof(double*));
    for(long h=1;h<=N;h++){ T[h]=malloc((q*h+1)*sizeof(double)); T[h][0]=0; for(long j=1;j<=q*h;j++) T[h][j]=psi((double)j/(double)(q*h)); }
    fprintf(stderr,"psi table %.1fs\n",(double)(clock()-t0)/CLOCKS_PER_SEC);
    gt=calloc(toff(N+1),sizeof(double));
    for(long k=1;k<=N;k++) for(long h=1;h<=k;h++){ if(gcdl(h,k)!=1) continue;
      long P=(long)q*h*k, a=0, ih=0, jk=0, nh=h, nk=k; int iq=0, jq=0; double psa=0, acc=0; const double *Tk=T[k], *Th=T[h];
      while(a<P){ long b = nh<nk? nh: nk; double psb = (b==nh)? Tk[ih+1] : Th[jk+1];
        double c=S[iq]*S[jq]; if(c!=0.0) acc+=c*(psb-psa);
        if(b==nh){ ih++; nh+=h; if(++iq==q) iq=0; } if(b==nk){ jk++; nk+=k; if(++jq==q) jq=0; } a=b; psa=psb; }
      gt[toff(k)+h-1]=acc/(double)P; }
    for(long h=1;h<=N;h++) free(T[h]); free(T);
  } else {
    long N=atol(argv[2]); const char* out0=argv[3]; (void)out0; Nn=N;
    double EG=0.57721566490153286061, LG=log(2*M_PI)-EG;
    double *V=calloc(toff(N+1),sizeof(double)); double *ct=malloc((N+1)*sizeof(double));
    for(long k=2;k<=N;k++){ long hk=(k-1)/2; if(hk==0) continue; for(long m=1;m<=hk;m++) ct[m]=1.0/tan(M_PI*(double)m/(double)k);
      for(long a=1;a<k;a++){ if(gcdl(a,k)!=1) continue; long r=a; double s=0; for(long m=1;m<=hk;m++){ s+=(double)(2*r-k)*ct[m]; r+=a; if(r>=k) r-=k; } V[toff(k)+a-1]=s/(double)k; } }
    fprintf(stderr,"V table %.1fs\n",(double)(clock()-t0)/CLOCKS_PER_SEC);
    gt=calloc(toff(N+1),sizeof(double));
    for(long k=1;k<=N;k++) for(long h=1;h<=k;h++){ if(gcdl(h,k)!=1) continue;
      double vhk = (k>1)? V[toff(k)+h-1] : 0.0;   // V[h,k]; h<k when k>1 (coprime)
      double vkh = (h>1)? V[toff(h)+(k%h)-1] : 0.0;
      gt[toff(k)+h-1]=(LG*(double)(h+k)+(double)(k-h)*log((double)h/(double)k)-M_PI*(vhk+vkh))/(2.0*(double)h*(double)k); }
    free(V); free(ct);
  }
  const char* out = !strcmp(mode,"per")? argv[4] : argv[3]; long N=Nn;
  FILE* f=fopen(out,"wb"); double *row=malloc(N*sizeof(double));
  for(long h=1;h<=N;h++){ for(long k=1;k<=N;k++){ long a=h<k?h:k, b=h<k?k:h; long d=gcdl(a,b); row[k-1]=gt[toff(b/d)+a/d-1]/(double)d; } fwrite(row,sizeof(double),N,f); }
  fclose(f); fprintf(stderr,"done %.1fs\n",(double)(clock()-t0)/CLOCKS_PER_SEC); return 0; }
