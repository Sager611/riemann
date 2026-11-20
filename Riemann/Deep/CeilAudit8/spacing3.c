// spacing3.c (ceilaudit8) = spacing2.c (explore49) + parity classes. Arg 10 cls: 0 odd q (= spacing2), 1 q=2 mod 4, 2 q=0 mod 4, 3 mixed (odd vs even pairs only, h-basis).
// Even q: Gauss sum G(a,k;q) vanishes on one parity of k (q=2 mod 4: k even; q=0 mod 4: k odd); in h'=h/2 basis, H'=H/2:
//  q=0 mod 4: x2=-abar/q, x1=b*x2 (b nearest EVEN integer to sqrt(KP)); conic coords (A,mod)=(a,q): B_delta(a,q)=P'-P
//  q=2 mod 4: x2=-(2a)bar/(q/2), x1=b*x2 (b nearest ODD integer); conic coords (2a,q/2): B_delta=P'-P
//  x3=2^1.5(2/3)/D, x4=sqrt2*eta/D, eta=sqrt(KP)-b in [-1,1). Mixed (cls 3): even arc -> 2 points (x1'/2 [+1/2], x2'/4, (2/3)/D, eta/D) in h-basis.
// spacing2.c - exact Bombieri-Iwaniec second-spacing counts B1 for zeta(1/2+it), Bourgain 1408.5794 Sec. 3 normalization.
// usage: spacing3 log10t nu [Qmult=1] [nshuffle=4] [seed=1] [Vroute=0] [verbose=0] [alpha] [dumpfile|-] [cls=0]
// [dumpfile]: write every 500th generic conic key (alpha beta g e n) and Q to this file
// [alpha=0.5]: M = floor(sqrt(K/2)) * t^(alpha-1/2) (block m ~ M below the top block when alpha < 1/2)
// Vroute=1: x2 window tightened by V = N/Q (Huxley's resonance-curve route, Bourgain (3.14)); verbose=1 prints 3 conic examples
// t = pi*K, K = round(10^log10t/pi); M = floor(sqrt(K/2)) (= sqrt(t/2pi)); N = round(M t^-nu); R = ceil(sqrt(2M^3/(NK)));
// Q = Qmult*R; H = N Q / R^2; intervals I_k=[kN,(k+1)N) in [M/2,M]; a/q = simplest rational in J(I) = f''(I)/2 = [K/(4(k+1)^2N^2), K/(4k^2N^2)];
// keep q odd, Q <= q < 2Q.  x1 = -b*inv(2a)/q, x2 = -inv(4a)/q (mod 1), x3 = (2/3)/D, x4 = eta/D, b+eta = sqrt(K a q) (eta in [-1/2,1/2)),
// D = q^{3/2} sqrt(K/(2 m0^3)), m0 = sqrt(Kq/(4a)).  Pair (I,I') counted when |x_j - x_j'| < 1/(2Y_j), Y = (6H, 6H^2, 6H^1.5, 6H^.5).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
typedef __int128 i128;
#define KEYMAX 1000000LL
typedef struct { int64_t a, q, A, qq, dv; int cl; double x1, x2, x3, x4, eta, D; } Rec; static int mixedonly=0;
static i128 isqrt128(i128 n){ if(n<2) return n; i128 r=(i128)sqrt((double)n); while(r*r>n) r--; while((r+1)*(r+1)<=n) r++; return r; }
static void simplest(i128 ln,i128 ld,i128 hn,i128 hd,i128*pn,i128*pd){
  i128 fl=ln/ld; if(fl*ld==ln){*pn=fl;*pd=1;return;} if((fl+1)*hd<=hn){*pn=fl+1;*pd=1;return;}
  i128 n2,d2; simplest(hd,hn-fl*hd,ld,ln-fl*ld,&n2,&d2); *pn=fl*n2+d2; *pd=n2; }
static int64_t inv(int64_t x,int64_t m){ int64_t g=m,a=x%m,u0=0,u1=1; if(a<0)a+=m; // returns x^-1 mod m
  while(a){ int64_t qq=g/a,t=g-qq*a; g=a; a=t; t=u0-qq*u1; u0=u1; u1=t; } if(g!=1) return -1; u0%=m; if(u0<0)u0+=m; return u0; }
static double fr(double x){ return x-floor(x); }
static double cd(double a,double b){ double d=fabs(a-b); d-=floor(d); return d<0.5?d:1-d; }
static uint64_t rs; static uint64_t rnd(void){ uint64_t z=(rs+=0x9E3779B97F4A7C15ULL); z=(z^(z>>30))*0xBF58476D1CE4E5B9ULL; z=(z^(z>>27))*0x94D049BB133111EBULL; return z^(z>>31); }
static int cmp3(const void*A,const void*B){ double a=((const Rec*)A)->x3,b=((const Rec*)B)->x3; return a<b?-1:a>b; }
static int cmpd(const void*A,const void*B){ const int64_t*a=A,*b=B; for(int i=0;i<5;i++){ if(a[i]<b[i]) return -1; if(a[i]>b[i]) return 1;} return 0; }
static double w1,w2,w3,w4; static int verbose=0, shown=0; static int64_t sq23=0,sq1234=0,ndef=0,ndef1234=0,ng0=0,maxabsn=0;
static int issq(i128 s){ if(s<0) return 0; i128 r=isqrt128(s); return r*r==s; }
static int cmpi(const void*A,const void*B){ int64_t a=*(const int64_t*)A,b=*(const int64_t*)B; return a<b?-1:a>b; }
// counts[mask] for pairs within the x3 window; mask bit0 = x1 ok, bit1 = x2 ok, bit2 = x4 ok
static void count(Rec*r,int64_t n,double*cnt,int conic,int64_t**keys,int64_t*nk,int64_t*bad){
  for(int m=0;m<8;m++) cnt[m]=0;
  for(int64_t i=0;i<n;i++) for(int64_t j=i+1;j<n && r[j].x3-r[i].x3<w3;j++){
    if(mixedonly && ((r[i].cl==0)==(r[j].cl==0))) continue;
    int m=(cd(r[i].x1,r[j].x1)<w1)|((cd(r[i].x2,r[j].x2)<w2)<<1)|((fabs(r[i].x4-r[j].x4)<w4)<<2); cnt[m]+=2;
    if(conic && (m&2)){ // verify A q - A'q' = -B_delta(A,q) and record conic key (delta, n)
      int64_t A=r[i].A,q=r[i].qq,A2=r[j].A,q2=r[j].qq; int64_t d=inv(A%q,q),d2=inv(A2%q2,q2);
      i128 g0=(i128)q2*d-(i128)d2*q; i128 qq=(i128)q*q2; i128 k=(g0>=0?(g0+qq/2)/qq:-((-g0+qq/2)/qq)); i128 dd2=d2+k*q2;
      i128 c=((i128)A*d-1)/q, c2=((i128)A2*dd2-1)/q2; i128 al=(i128)A2*d-c2*q, be=c2*A-(i128)A2*c, g=(i128)q2*d-dd2*q, e=dd2*A-(i128)q2*c;
      i128 lhs=(i128)A*q-(i128)A2*q2, B=al*g*A*A+2*be*g*A*q+be*e*q*q;
      if(lhs!=-B || al*e-be*g!=1) (*bad)++;
      int gen=1; { i128 s=-be*g; if(s<=0){ ndef++; gen=0; if(m==7) ndef1234++; if(g==0) ng0++; } else if(issq(s)){ sq23++; gen=0; if(m==7) sq1234++; } int64_t an=(int64_t)(lhs>=0?lhs:-lhs)/r[i].dv; if(an>maxabsn) maxabsn=an; }
      if(verbose && shown<3){ shown++; printf("EXAMPLE a/q=%lld/%lld a'/q'=%lld/%lld delta=[[%lld,%lld],[%lld,%lld]] Aq-A'q'=%lld B_delta(A,q)=%lld disc=%lld x4ok=%d x1ok=%d\n",(long long)r[i].a,(long long)q,(long long)r[j].a,(long long)q2,(long long)al,(long long)be,(long long)g,(long long)e,(long long)lhs,(long long)B,(long long)(-4*be*g),(m>>2)&1,m&1); }
      if(keys && gen && *nk<KEYMAX){ int64_t*K5=*keys+5*(*nk); K5[0]=(int64_t)al;K5[1]=(int64_t)be;K5[2]=(int64_t)g;K5[3]=(int64_t)e;K5[4]=(int64_t)(lhs/r[i].dv); (*nk)++; }
    }
  }
}
int main(int argc,char**argv){
  double L=atof(argv[1]), nu=atof(argv[2]), Qm=argc>3?atof(argv[3]):1; int ns=argc>4?atoi(argv[4]):4; rs=argc>5?atoll(argv[5]):1; int Vr=argc>6?atoi(argv[6]):0; verbose=argc>7?atoi(argv[7]):0; double alpha=argc>8?atof(argv[8]):0.5; const char*dumpf=(argc>9&&strcmp(argv[9],"-"))?argv[9]:NULL; int cls=argc>10?atoi(argv[10]):0; if(cls==3) mixedonly=1;
  double t=pow(10,L); i128 K=(i128)llround(t/M_PI/1e6)*1000000; // K ~ t/pi (rounded to 1e6 so i128 conversion is exact)
  if(t<1e13) K=(i128)llround(t/M_PI); t=M_PI*(double)K;
  int64_t M=(int64_t)floor(sqrt((double)K/2)*pow(t,alpha-0.5)); int64_t N=llround(M*pow(t,-nu)); if(N<1)N=1;
  int64_t R=(int64_t)ceil(sqrt(2.0*(double)M*M*M/((double)N*(double)K))); int64_t Q=(int64_t)llround(Qm*R);
  double H=(double)N*Q/((double)R*R); double Hc=(cls==1||cls==2)?H/2:H; w1=1/(12*Hc); w2=1/(12*Hc*Hc); double V=Vr?(double)N/(Qm*1.0*R):1; w2/=V; w3=1/(12*pow(Hc,1.5)); w4=1/(12*sqrt(Hc));
  int64_t k0=(M/2+N-1)/N, k1=M/N-1, nI=k1-k0+1; Rec*r=malloc(sizeof(Rec)*(size_t)(nI>0?2*nI:1)); int64_t ncl[3]={0,0,0}; int64_t n=0, nodd=0, nq[4]={0,0,0,0}, outside=0;
  for(int64_t k=k0;k<=k1;k++){ i128 lo_d=(i128)4*(k+1)*(k+1)*N*N, hi_d=(i128)4*k*k*N*N, pn,pd; simplest(K,lo_d,K,hi_d,&pn,&pd);
    int64_t a=(int64_t)pn,q=(int64_t)pd; if(q<Q) nq[0]++; else if(q<2*Q) nq[1]++; else nq[2]++;
    if(!(q>=Q && q<2*Q)) continue; int qc=(q&1)?0:((q%4==2)?1:2); ncl[qc]++; if(cls<3 && qc!=cls) continue; nodd++;
    double m0=sqrt((double)K*q/(4.0*a)); if(m0<k*N-1e-6||m0>(k+1)*N+1e-6) outside++;
    i128 PP=K*a*q, s=isqrt128(PP); double frac=(double)(PP-s*s)/(sqrt((double)PP)+(double)s); i128 b=s; double eta=frac;
    if(qc==0){ if(frac>=0.5){ b=s+1; eta=frac-1; } } else { int par=(qc==1); if((int)(s&1)!=par){ b=s+1; eta=frac-1; } } // nearest integer of the allowed parity
    double D=pow((double)q,1.5)*sqrt((double)K/(2*m0*m0*m0));
    if(qc==0){ int64_t bm=(int64_t)(b%q); Rec x={a,q,4*a,q,4,0, fr(-(double)((bm*inv((2*a)%q,q))%q)/q), fr(-(double)inv((4*a)%q,q)/q), (2.0/3)/D, eta/D, eta, D}; r[n++]=x; }
    else { int64_t md=(qc==1)?q/2:q, A=(qc==1)?2*a:a, ia=inv(A%md,md), bm=(int64_t)(b%md); double x2p=fr(-(double)ia/md), x1p=fr(-(double)((bm*ia)%md)/md);
      if(cls<3){ Rec x={a,q,A,md,1,qc, x1p, x2p, pow(2,1.5)*(2.0/3)/D, sqrt(2.0)*eta/D, eta, D}; r[n++]=x; }
      else { Rec x={a,q,A,md,1,qc, x1p/2, x2p/4, (2.0/3)/D, eta/D, eta, D}; r[n++]=x; Rec y=x; y.x1=fr(x1p/2+0.5); r[n++]=y; } } }
  qsort(r,n,sizeof(Rec),cmp3);
  double c[8]; int64_t nk=0,bad=0; int64_t*keys=malloc(sizeof(int64_t)*5*KEYMAX); // conic keys, bounded (KEYMAX*40 bytes = 40 MB)
  count(r,n,c,1,&keys,&nk,&bad);
  double c3=0,c23=0,c234=0,c1234=0,c34=0,c123=0; for(int m=0;m<8;m++){ c3+=c[m]; if(m&2)c23+=c[m]; if((m&6)==6)c234+=c[m]; if(m==7)c1234+=c[m]; if(m&4)c34+=c[m]; if((m&3)==3)c123+=c[m]; }
  // conic multiplicities among (x2,x3)-pairs
  if(dumpf){ FILE*fd=fopen(dumpf,"w"); fprintf(fd,"Q %lld\n",(long long)Q); for(int64_t i=0;i<nk;i+=500) fprintf(fd,"%lld %lld %lld %lld %lld\n",(long long)keys[5*i],(long long)keys[5*i+1],(long long)keys[5*i+2],(long long)keys[5*i+3],(long long)keys[5*i+4]); fclose(fd); }
  qsort(keys,nk,5*sizeof(int64_t),cmpd); int64_t ud=0,un=0,maxd=0,maxn=0,run=0,runn=0;
  for(int64_t i=0;i<nk;i++){ int newd=(i==0)||memcmp(keys+5*i,keys+5*(i-1),4*sizeof(int64_t)); int newn=newd||keys[5*i+4]!=keys[5*(i-1)+4];
    if(newd){ ud++; run=0;} run++; if(run>maxd)maxd=run; if(newn){ un++; runn=0;} runn++; if(runn>maxn)maxn=runn; }
  printf("cls=%d Hc=%.3f classes(odd,2mod4,0mod4)=%lld,%lld,%lld records=%lld ",cls,Hc,(long long)ncl[0],(long long)ncl[1],(long long)ncl[2],(long long)n); printf("alpha=%.3f ",alpha); printf("log10t=%.2f nu=%.4f K=%.6e M=%lld N=%lld R=%lld Q=%lld H=%.3f intervals=%lld q<Q:%lld Q<=q<2Q:%lld q>=2Q:%lld kept(odd)=%lld m0_outside=%lld\n",L,nu,(double)K,(long long)M,(long long)N,(long long)R,(long long)Q,H,(long long)nI,(long long)nq[0],(long long)nq[1],(long long)nq[2],(long long)n,(long long)outside);
  printf("TRUE  diag=%lld off: x3=%.0f x2x3=%.0f x2x3x4=%.0f x1x2x3=%.0f all4=%.0f x3x4=%.0f\n",(long long)n,c3,c23,c234,c123,c1234,c34);
  printf("CONIC x2x3-pairs/2=%lld identity_fail=%lld distinct_delta=%lld max_pairs_per_delta=%lld distinct_(delta,n)=%lld max_per_conic=%lld\n",(long long)nk,(long long)bad,(long long)ud,(long long)maxd,(long long)un,(long long)maxn);
  { int64_t*nn=malloc(sizeof(int64_t)*(size_t)(nk>0?nk:1)); for(int64_t i=0;i<nk;i++) nn[i]=keys[5*i+4]; qsort(nn,nk,sizeof(int64_t),cmpi); int64_t un2=0,run2=0,maxn2=0; for(int64_t i=0;i<nk;i++){ if(i==0||nn[i]!=nn[i-1]){un2++;run2=0;} run2++; if(run2>maxn2)maxn2=run2; }
    double Pmin=1e300,Pmax=0; for(int64_t i=0;i<n;i++){ double P=(double)r[i].a*r[i].q; if(P<Pmin)Pmin=P; if(P>Pmax)Pmax=P; }
    printf("DEGEN x2x3-pairs on square-disc (line-pair) conics=%lld (of %lld, all4 among them=%lld) disc<=0=%lld (of which g=0 translations a'=a+jq, q'=q: %lld; all4 among disc<=0=%lld) [keys/CONIC stats: generic conics only] | n=aq-a'q': max|n|=%lld distinct n=%lld max pairs per n=%lld, P-range=[%.3g,%.3g] => admissible-P density per n ~ %.3g\n",(long long)sq23,(long long)nk,(long long)sq1234,(long long)ndef,(long long)ng0,(long long)ndef1234,(long long)maxabsn,(long long)un2,(long long)maxn2,Pmin,Pmax,(double)nk/((double)(un2>0?un2:1)*(Pmax-Pmin))); free(nn); }
  double s23=0,s234=0,s1234=0,s123=0,s34=0;
  for(int it=0;it<ns;it++){ // random model: permute x1, x2, eta independently (x3, D stay tied to P)
    for(int64_t i=n-1;i>0;i--){ int64_t j=rnd()%(i+1); double t1=r[i].x1; r[i].x1=r[j].x1; r[j].x1=t1; }
    for(int64_t i=n-1;i>0;i--){ int64_t j=rnd()%(i+1); double t2=r[i].x2; r[i].x2=r[j].x2; r[j].x2=t2; }
    for(int64_t i=n-1;i>0;i--){ int64_t j=rnd()%(i+1); double te=r[i].eta; r[i].eta=r[j].eta; r[j].eta=te; }
    for(int64_t i=0;i<n;i++) r[i].x4=r[i].eta/r[i].D;
    count(r,n,c,0,NULL,NULL,&bad); for(int m=0;m<8;m++){ if(m&2)s23+=c[m]; if((m&6)==6)s234+=c[m]; if(m==7)s1234+=c[m]; if((m&3)==3)s123+=c[m]; if(m&4)s34+=c[m]; } }
  if(ns) printf("RAND  off(mean of %d shuffles): x2x3=%.1f x2x3x4=%.1f x1x2x3=%.1f all4=%.1f x3x4=%.1f\n",ns,s23/ns,s234/ns,s123/ns,s1234/ns,s34/ns);
  double HW=(1/(3*H*H))*((double)Q*Q/(6.0*R*R*H*H))*pow((double)M/N,2)*pow((double)Q/R,4)/V; // [H-W] shape Δ1Δ2(M/N)^2(Q/R)^4 (/V)
  printf("SHAPE V=%.3f x2x3/HWshape=%.5f all4/(HWshape*H^-2)=%.5f all4/(HWshape*Delta4^(2/3))=%.5f\n",V,c23/HW,c1234/(HW/(H*H)),c1234/(HW*pow((double)Q*Q/(6.0*R*R*H),2.0/3)));
  printf("RATIO off_all4/diag=%.4g  true/rand: x2x3=%.3f x2x3x4=%.3f all4=%.3f  saving(all4/x2x3)=%.4g  H^2=%.4g\n",c1234/n,c23/(s23/ns),c234/(s234/ns),c1234/(s1234/ns),c1234/c23,H*H);
  return 0; }
