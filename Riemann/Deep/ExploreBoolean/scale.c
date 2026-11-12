
// scale.c — knapsack down-sets D_X = {S : prod_{i in S} q_i <= X} for generalized primes q_i, at scale.
// usage: scale <sys> <X> <seed> <p1> <p2> <mode> <K> <alpha>
//  sys 0 primes; 1 p+1; 2 smooth points (density 1/log u, u >= U0, primes below U0);
//      3 Poisson points, intensity 1/log t on [2,X]; 4 iid perturbed p*exp(U[-p1,p1]);
//      5 shift p*exp(-p1); 6 planted zero: density (1 - 2 u^{b-1} cos(g log u))/log u for u >= U0, b=p1, g=p2
//  mode 0: statistics; mode 1: K random sign vectors z (z_i = -1 w.p. alpha): F(z) = sum_{m in D} prod z
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
static double *q; static long nq; static double X;
static uint64_t rs=0x9E3779B97F4A7C15ULL;
static inline uint64_t xr(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static inline double ur(void){ return (xr()>>11)*(1.0/9007199254740992.0); }
static uint8_t *sv; // odd-only bit sieve
static long SV;
static void sieve(long lim){ SV=lim; long nb=(lim/2)/8+2; sv=calloc(nb,1); // bit i <-> 2i+1 composite
  for(long i=1; (2*i+1)*(2*i+1)<=lim; i++) if(!(sv[i>>3]>>(i&7)&1)){ long p=2*i+1; for(long j=(p*p)/2; 2*j+1<=lim; j+=p) sv[j>>3]|=1<<(j&7);} }
static int isp(long k){ if(k<2) return 0; if(k==2) return 1; if(!(k&1)) return 0; long i=k/2; return !(sv[i>>3]>>(i&7)&1); }
static int cmpd(const void*a,const void*b){ double x=*(const double*)a,y=*(const double*)b; return x<y?-1:x>y; }
static long cap; static void push(double v){ if(nq>=cap){ cap=cap?cap*2:1<<20; q=realloc(q,cap*sizeof(double)); } q[nq++]=v; }
static double U0=40.0, PB=0.8, PG=10.0; static int sysid;
static double gdens(double u){ if(sysid==6) return (1.0-2.0*pow(u,PB-1.0)*cos(PG*log(u)))/log(u); return 1.0/log(u); }
static void gen_cont(void){ // primes below U0, then points at Pi_c = i - 1/2
  long c=0; for(long k=2;k<U0;k++) if(isp(k)){ push((double)k); c++; }
  double t=U0, Pi=(double)c, target=c+0.5;
  while(1){ double D=target-Pi; double h=D/gdens(t);
    for(int it=0;it<4;it++){ double ga=(gdens(t)+4*gdens(t+h/2)+gdens(t+h))/6.0; h=D/ga; }
    t+=h; Pi=target; target+=1.0; if(t>X) break; push(t); }
}
// ---- statistics mode
static int ncp; static double cpt[20]; static long cN[20], cM[20];
static int32_t *binN,*binM; static long nbins; static double bsc;
static int32_t *bc,*df,*cc; static long A[64]; static long stk[64]; static int dep; static long nodes;
static inline long ub(double lim){ long lo=0,hi=nq; while(lo<hi){ long mid=(lo+hi)>>1; if(q[mid]>lim) hi=mid; else lo=mid+1; } return lo; }
static void visit(double P,int s){
  int c=0; while(P>cpt[c]) c++; cN[c]++; cM[c]+=s;
  long b=(P<=1.0)?0:(long)(log(P)*bsc); if(b>=nbins) b=nbins-1; binN[b]++; binM[b]+=s;
  A[dep]++; nodes++;
  long J=ub(X/P); if(J<nq){ df[J]++; for(int t=0;t<dep;t++) if(stk[t]>=J) bc[stk[t]]--; }
  for(int t=0;t<dep;t++) cc[stk[t]]++;
}
static void dfs(long st,double P,int s){ visit(P,s); for(long i=st;i<nq;i++){ double P2=P*q[i]; if(P2>X) break; stk[dep++]=i; dfs(i+1,P2,-s); dep--; } }
// ---- sampling mode
static int8_t *z; static long Fz;
static void dfz(long st,double P,int s){ Fz+=s; for(long i=st;i<nq;i++){ double P2=P*q[i]; if(P2>X) break; dfz(i+1,P2,s*z[i]); } }
static int cmpl(const void*a,const void*b){ long x=*(const long*)a,y=*(const long*)b; return x<y?-1:x>y; }
static double gp1, gp2;
static void gen(double p1,double p2){
  if(sysid==6){ PB=p1; PG=p2; U0=40.0; }
  long lim=(long)(X*exp(fabs(p1))*1.01)+100; if(sysid!=4 && sysid!=5) lim=(long)X+100; sieve(lim);
  if(sysid==0){ for(long k=2;k<=(long)X;k++) if(isp(k)) push((double)k); }
  else if(sysid==1){ for(long k=2;k+1<=(long)X;k++) if(isp(k)) push((double)(k+1)); }
  else if(sysid==2||sysid==6) gen_cont();
  else if(sysid==3){ for(long k=2;k<(long)X;k++){ double lam=1.0/log(k+0.5); double u=ur(), pr=exp(-lam), cdf=pr; int m=0; while(u>cdf){ m++; pr*=lam/m; cdf+=pr; }
      double tmp[16]; for(int j=0;j<m && j<16;j++) tmp[j]=k+ur(); qsort(tmp,m<16?m:16,sizeof(double),cmpd); for(int j=0;j<m&&j<16;j++) if(tmp[j]<=X) push(tmp[j]); } }
  else if(sysid==4){ for(long k=2;k<=lim;k++) if(isp(k)){ double v=k*exp(p1*(2*ur()-1)); if(v<=X) push(v);} qsort(q,nq,sizeof(double),cmpd); }
  else if(sysid==5){ for(long k=2;k<=lim;k++) if(isp(k)){ double v=k*exp(-p1); if(v<=X) push(v);} }
  free(sv);
}
static double *HA,*HB; static long NA,NB,only,cntv;
static void hd(long st,double PA,double PB2,int ok){ cntv++; if(!ok || PB2>X) only++;
  for(long i=st;i<NA;i++){ double P2=PA*HA[i]; if(P2>X) break; int ok2=ok && i<NB; hd(i+1,P2, ok2?PB2*HB[i]:1e300, ok2); } }
int main(int argc,char**argv){
  sysid=atoi(argv[1]); X=atof(argv[2]); rs^=(uint64_t)atol(argv[3])*0xD1B54A32D192ED03ULL; for(int i=0;i<10;i++) xr();
  double p1=atof(argv[4]), p2=atof(argv[5]); int mode=atoi(argv[6]); long K=atol(argv[7]); double alpha=atof(argv[8]);
  gen(p1,p2);
  if(mode==3){ // Hamming distance between this down-set and the prime down-set (i-th point <-> i-th prime)
    double *qs=q; long ns=nq; q=NULL; nq=0; cap=0; int keep=sysid; sysid=0; gen(0,0); sysid=keep; double *pp=q; long np=nq;
    long mm=ns<np?ns:np; double dmax=0,dmax3=0; for(long i=0;i<mm;i++){ double d=fabs(log(qs[i])-log(pp[i])); if(d>dmax) dmax=d; if(pp[i]>1000 && d>dmax3) dmax3=d; }
    HA=pp; NA=np; HB=qs; NB=ns; only=0; hd(0,1.0,1.0,1); long onlyP=only, NP=cntv; cntv=0;
    HA=qs; NA=ns; HB=pp; NB=np; only=0; hd(0,1.0,1.0,1); long onlyQ=only, NQ=cntv;
    printf("HAMMING sys=%d X=%.0f n_sys=%ld n_primes=%ld max_i|log q_i-log p_i|=%.4f (p_i>1000: %.4f) |D_p|=%ld |D_q|=%ld |D_p minus D_q|=%ld |D_q minus D_p|=%ld rel=%.3e\n",
      keep,X,ns,np,dmax,dmax3,NP,NQ,onlyP,onlyQ,(double)(onlyP+onlyQ)/NP);
    return 0; }
  if(q[0]<=1.0){ fprintf(stderr,"q0<=1\n"); return 1; }
  if(mode==1){
    z=malloc(nq); long *res=malloc(sizeof(long)*(K+1));
    for(long i=0;i<nq;i++) z[i]=-1; Fz=0; dfz(0,1.0,1); long Mtop=Fz;
    for(long k=0;k<K;k++){ for(long i=0;i<nq;i++) z[i]=(ur()<alpha)?-1:1; Fz=0; dfz(0,1.0,1); res[k]=Fz; }
    double s1=0,s2=0; long ge=0; for(long k=0;k<K;k++){ s1+=labs(res[k]); s2+=(double)res[k]*res[k]; if(labs(res[k])>=labs(Mtop)) ge++; }
    qsort(res,K,sizeof(long),cmpl);
    for(long i=0;i<nq;i++) z[i]=1; Fz=0; dfz(0,1.0,1); long Ntot=Fz;
    printf("SAMPLE sys=%d X=%.0f n=%ld N=%ld alpha=%.3f K=%ld M=%ld E|F|=%.2f rms=%.2f sqrtN=%.2f frac(|F|>=|M|)=%.4f median=%ld q05=%ld q95=%ld min=%ld max=%ld\n",
      sysid,X,nq,Ntot,alpha,K,Mtop,s1/K,sqrt(s2/K),sqrt((double)Ntot),(double)ge/K,res[K/2],res[K/20],res[K-1-K/20],res[0],res[K-1]);
    return 0;
  }
  ncp=0; for(double c=1000; c<X*0.999; c*=10) cpt[ncp++]=c; cpt[ncp++]=X; cpt[ncp]=1e300;
  nbins=4000L*(long)ceil(log10(X)); bsc=nbins/log(X)*(1-1e-15); binN=calloc(nbins,4); binM=calloc(nbins,4);
  bc=calloc(nq+1,4); df=calloc(nq+1,4); cc=calloc(nq+1,4);
  dfs(0,1.0,1);
  long run=0; for(long i=0;i<nq;i++){ run+=df[i]; bc[i]+=run; }
  long Ncum=0,Mcum=0; printf("sys=%d X=%.0f n=%ld q[0..3]=%.4f %.4f %.4f %.4f nodes=%ld\n",sysid,X,nq,q[0],q[1],q[2],q[3],nodes);
  for(int c=0;c<ncp;c++){ Ncum+=cN[c]; Mcum+=cM[c]; printf("  x=%.0f N=%ld M=%ld M/sqrt(x)=%.4f N/x=%.6f\n",cpt[c],Ncum,Mcum,Mcum/sqrt(cpt[c]),Ncum/cpt[c]); }
  long N=Ncum, M=Mcum;
  // checks and Boolean statistics at X
  long bad=0; double sb=0,sb2=0; long bmax=0, imax=0; for(long i=0;i<nq;i++){ if(bc[i]!=N-2*cc[i]) bad++; sb+=bc[i]; sb2+=(double)bc[i]*bc[i]; if(bc[i]>bmax){bmax=bc[i]; imax=i;} }
  double lna = nq*log(2.0)-log((double)N);
  double om=0,om2=0; for(int k=0;k<64;k++){ om+=k*(double)A[k]; om2+=(double)k*k*A[k]; } om/=N; om2=om2/N-om*om;
  long symtop=0; for(int k=0;k<64;k++) symtop+=(k&1?-1:1)*A[k];
  printf("  STATS N=%ld M=%ld symtop=%ld |M|/sqrtN=%.3f chowcheck_bad=%ld lev1ratio=%.4f I/mu=%.4f maxInf/mu=%.4f(q=%.2f) betti/N=%.4f omega_mean=%.4f omega_var=%.4f\n",
    N,M,symtop,labs(M)/sqrt((double)N),bad,sb2/(2.0*N*(double)N*lna),2*sb/N,2.0*bmax/N,q[imax],(double)bc[0]/N,om,om2);
  printf("  A_k:"); for(int k=0;k<64;k++) if(A[k]) printf(" %ld",A[k]); printf("\n");
  // M(x) curve from bins: per decade max |M|/sqrt(x) and max |M|
  long cumM=0; double edge; int dec=-1; double mx=0, mxr=0;
  FILE*fo=NULL; if(argc>9){ fo=fopen(argv[9],"w"); }
  for(long b=0;b<nbins;b++){ cumM+=binM[b]; edge=exp((b+1)/bsc); int d=(int)floor(log10(edge)+1e-12);
    if(d!=dec){ if(dec>=3) printf("  decade 1e%d: max|M|=%.0f max|M|/sqrt(x)=%.3f\n",dec,mx,mxr); dec=d; mx=0; mxr=0; }
    double am=fabs((double)cumM); if(am>mx) mx=am; if(am/sqrt(edge)>mxr) mxr=am/sqrt(edge);
    if(fo && edge>=1000) fprintf(fo,"%.6e %ld\n",edge,cumM); }
  if(dec>=3) printf("  decade 1e%d: max|M|=%.0f max|M|/sqrt(x)=%.3f\n",dec,mx,mxr);
  if(fo) fclose(fo);
  return 0;
}
