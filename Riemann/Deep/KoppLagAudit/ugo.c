/* ugo.c — ugoaudit independent check of unit-generated-orders (KL Problem 1.7). Own code, no reuse of Riemann/Deep/KoppLag.
   mode "sieve": split-prime sieve (Cor. B) over n in [1001, 10^9+1], both N; Legendre tables by Euler's criterion; early exit.
   mode "class NMAX": for all n <= NMAX (n>=3 for N=+1, n>=1 for N=-1): reduced primitive forms, rho-cycles -> h+, h;
     2-torsion decided WITHOUT composition: genus theory (|Cl+/Cl+^2| = 2^(mu-1), principal genus = Cl+^2) gives
     narrow2 <=> h+ = 2^(mu-1);  wide2 <=> h = |Cl/Cl^2| = 2^(mu-1) if J in Cl+^2 else 2^(mu-2), J = [(-1,b0,c)],
     J in Cl+^2 <=> every assigned genus character is +1 at -1.  Also: Lemma A (Lagrange form) on the principal cycle;
     analytic class number formula h = sqrt(D) L(1,chi) prod_{p|f}(1-chi(p)/p) / (2 R) for n <= ANMAX and listed n.
   cc -O2 -o ugo ugo.c -lm */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
typedef long long ll; typedef unsigned long long ull;
static void die(const char *m, ll x){ fprintf(stderr,"ASSERT %s %lld\n",m,x); exit(3); }
static ll isq(ll x){ ll r=(ll)sqrtl((long double)x); while(r*r>x) r--; while((r+1)*(r+1)<=x) r++; return r; }
static ll gcdl(ll a, ll b){ if(a<0)a=-a; if(b<0)b=-b; while(b){ ll t=a%b; a=b; b=t; } return a; }
static ull powmod(ull b, ull e, ull m){ ull r=1%m; b%=m; while(e){ if(e&1) r=(unsigned __int128)r*b%m; b=(unsigned __int128)b*b%m; e>>=1; } return r; }
static int *spf; static ll SPFN;
static int factor(ll m, ll *ps, int *es){ int k=0; if(m<=SPFN){ while(m>1){ ll p=spf[m]; int e=0; while(m%p==0){m/=p;e++;} ps[k]=p; es[k++]=e; } return k; }
  for(ll p=2;p*p<=m;p+=(p==2?1:2)){ if(m%p==0){ int e=0; while(m%p==0){m/=p;e++;} ps[k]=p; es[k++]=e; } } if(m>1){ ps[k]=m; es[k++]=1; } return k; }
/* Jacobi (a/n), n odd > 0 */
static int jacobi(ll a, ll n){ int t=1; a%=n; if(a<0)a+=n; while(a){ while(!(a&1)){ a>>=1; ll r=n&7; if(r==3||r==5) t=-t; } ll x=a; a=n; n=x; if((a&3)==3&&(n&3)==3) t=-t; a%=n; } return n==1?t:0; }
static int kron(ll D, ll a){ /* Kronecker (D/a), a >= 1, D = 0,1 mod 4 */ int t=1; while(!(a&1)){ a>>=1; if(!(D&1)) return 0; ll r=((D%8)+8)%8; if(r==3||r==5) t=-t; } return t*jacobi(D,a); }
typedef struct { ll a,b,c; } F;
static int cmpF(const void *x, const void *y){ const F *p=x,*q=y; if(p->a!=q->a) return p->a<q->a?-1:1; if(p->b!=q->b) return p->b<q->b?-1:1; return 0; }
static F *fs; static int nf, capf;
static void addf(ll a, ll b, ll c){ if(nf==capf){ capf=capf?2*capf:1024; fs=realloc(fs,capf*sizeof(F)); } fs[nf].a=a; fs[nf].b=b; fs[nf].c=c; nf++; }
static int findf(ll a, ll b){ int lo=0,hi=nf-1; while(lo<=hi){ int m=(lo+hi)/2; if(fs[m].a==a&&fs[m].b==b) return m; if(fs[m].a<a||(fs[m].a==a&&fs[m].b<b)) lo=m+1; else hi=m-1; } return -1; }
static void rho(F f, ll D, ll s, F *g){ ll C=f.c<0?-f.c:f.c, m=2*C, r; if(C>s){ r=((-f.b)%m+m)%m; if(r>C) r-=m; } else r=s-(((s+f.b)%m+m)%m);
  ll num=r*r-D; if(num%(4*f.c)) die("rho",D); g->a=f.c; g->b=r; g->c=num/(4*f.c); }
static int isred(F f, ll s){ ll A=f.a<0?-f.a:f.a; return f.b>0 && f.b<=s && 2*A+f.b>s && 2*A-f.b<=s; }
/* divisors of M in [lo,hi] */
static ll dv[100000]; static int ndv;
static void divs(ll M){ ll ps[64]; int es[64]; int k=factor(M,ps,es); ndv=1; dv[0]=1; for(int i=0;i<k;i++){ int old=ndv; ll pp=1; for(int e=1;e<=es[i];e++){ pp*=ps[i]; for(int j=0;j<old;j++) dv[ndv++]=dv[j]*pp; } } }

static int mu_of(ll D, ll *ps, int *es, int k, int *JinSq){ /* KL l.1527-1537 mu; JinSq: all assigned characters = +1 at -1 (Cox Thm 3.15 chars) */
  int r=0, e2=0, ok=1; for(int i=0;i<k;i++){ if(ps[i]==2) e2=es[i]; else { r++; if(ps[i]%4==3) ok=0; } }
  int mu;
  if(D%4==1) mu=r;
  else { ll m=D/4; ll m4=((m%4)+4)%4, m8=((m%8)+8)%8;
    if(m4==1) mu=r;                       /* D = 4 mod 16: no 2-adic character */
    else if(m4==3){ mu=r+1; ok=0; }       /* D = 12 mod 16: delta, delta(-1) = -1 */
    else if(m8==2){ mu=r+1; }             /* D = 8 mod 32: epsilon, epsilon(-1) = +1 */
    else if(m8==6){ mu=r+1; ok=0; }       /* D = 24 mod 32: delta*epsilon, = -1 at -1 */
    else if(m8==4){ mu=r+1; ok=0; }       /* D = 16 mod 32: delta */
    else { mu=r+2; ok=0; }                /* D = 0 mod 32: delta, epsilon */
    /* cross-check with KL's D mod 16/32 table */
    ll d16=D%16, d32=D%32; int muKL = (d16==4)? r : ((d16==8||d16==12||d32==16)? r+1 : (d32==0? r+2 : -99));
    if(muKL!=mu) die("mu table",D);
  }
  (void)e2; *JinSq=ok; return mu; }

static void sieve(void){
  const ull DS=1000000000000000000ULL; int PM=23000; char *comp=calloc(PM+1,1); int pr[3000], np=0;
  for(int i=2;i<=PM;i++) if(!comp[i]){ pr[np++]=i; for(int j=2*i;j<=PM;j+=i) comp[j]=1; }
  ull p4[3000]; for(int k=0;k<np;k++){ ull p=pr[k]; p4[k]=4*p*p*p*p; }
  if(!(p4[np-1] > DS+8)) die("prime list short",np);
  for(int N=1;N>=-1;N-=2){
    char **tab=malloc(np*sizeof(char*));
    for(int k=1;k<np;k++){ ull p=pr[k]; tab[k]=malloc(p); for(ull r=0;r<p;r++){ ll v=((ll)((r*r)%p)-4*N)%(ll)p; if(v<0) v+=p; tab[k][r]= (v!=0 && powmod(v,(p-1)/2,p)==1); } }
    ull tested=0, surv=0, nlast=0;
    for(ull n=1001;;n++){
      ull D = (N==1)? n*n-4 : n*n+4;
      if(n>1000000001ULL) break;        /* n <= 10^9+1 covers every Delta_n < 10^18 + 2*10^9 */
      tested++; nlast=n;
      if(D%8==1) continue;              /* 2 splits (never in these families) */
      int ok=1; unsigned nn=(unsigned)n;
      for(int k=1;k<np;k++){ if(p4[k]>=D) break; if(tab[k][nn%(unsigned)pr[k]]){ ok=0; break; } }
      if(ok){ surv++; printf("SURV N=%+d n=%llu D=%llu\n",N,n,D); }
    }
    printf("SIEVE N=%+d tested %llu n in [1001,%llu] survivors %llu\n",N,tested,nlast,surv); fflush(stdout);
    for(int k=1;k<np;k++) free(tab[k]); free(tab);
  }
}

static double Lchi(ll q){ /* L(1,(q/.)) = -(2/sqrt q) sum_{a<q/2} chi(a) log sin(pi a/q), chi even */
  long double S=0; for(ll a=1;2*a<q;a++){ int x=kron(q,a); if(x) S+= x*logl(sinl(3.14159265358979323846264338327950288L*a/q)); }
  return (double)(-2.0L*S/sqrtl((long double)q)); }

static ll SAMPLE=0;  /* analytic check also on every n > 1000 with n % SAMPLE == 0 */
static void classes(ll NMAX, ll ANMAX){
  SPFN=(NMAX*NMAX+4)/4+2; spf=calloc(SPFN+1,sizeof(int));
  for(ll i=2;i<=SPFN;i++) if(!spf[i]){ for(ll j=i;j<=SPFN;j+=i) if(!spf[j]) spf[j]=(int)i; }
  ll norders=0, nw=0, nn2=0, lemA=0; double maxerr=0; ll nan=0;
  for(int N=1;N>=-1;N-=2) for(ll n=(N==1?3:1); n<=NMAX; n++){
    ll D=n*n-4*N, s=isq(D); if(s*s==D) die("square",D);
    nf=0;
    for(ll b=s; b>0; b--){ if(((b-D)%2+2)%2) continue; ll M=(D-b*b)/4, lo=(s-b)/2+1, hi=(s+b)/2;
      divs(M); for(int j=0;j<ndv;j++){ ll a=dv[j]; if(a<lo||a>hi) continue; ll c=M/a; if(gcdl(gcdl(a,b),c)!=1) continue; addf(a,b,-c); addf(-a,b,c); } }
    qsort(fs,nf,sizeof(F),cmpF);
    for(int i=0;i<nf;i++){ if(!isred(fs[i],s)) die("notred",D); if(fs[i].b*fs[i].b-4*fs[i].a*fs[i].c!=D) die("disc",D); if(i&&fs[i].a==fs[i-1].a&&fs[i].b==fs[i-1].b) die("dup",D); }
    int *pm=malloc(nf*sizeof(int)), *indeg=calloc(nf,sizeof(int)), *cy=malloc(nf*sizeof(int));
    for(int i=0;i<nf;i++){ F g; rho(fs[i],D,s,&g); int j=findf(g.a,g.b); if(j<0) die("rho out",D); if(fs[j].c!=g.c) die("rho c",D); pm[i]=j; indeg[j]++; cy[i]=-1; }
    for(int i=0;i<nf;i++) if(indeg[i]!=1) die("not perm",D);
    int hp=0; for(int i=0;i<nf;i++) if(cy[i]<0){ int j=i; do{ cy[j]=hp; j=pm[j]; }while(j!=i); hp++; }
    ll b0 = ((s-D)%2==0)? s : s-1; int iP=findf(1,b0), iJ=findf(-1,b0); if(iP<0||iJ<0) die("P/J",D);
    int pid=cy[iP], jid=cy[iJ]; ll h = (pid==jid)? hp : hp/2; if(pid!=jid && hp%2) die("odd h+",D);
    /* Lemma A, Lagrange form: on the principal cycle every |a| with 4a^2 < D is 1 */
    for(int i=0;i<nf;i++) if(cy[i]==pid){ ll A=fs[i].a<0?-fs[i].a:fs[i].a; if(4*A*A<D && A!=1 && D>=100){ lemA++; printf("LEMMA_A_FAIL N=%+d n=%lld D=%lld a=%lld\n",N,n,D,fs[i].a); } }
    ll ps[64]; int es[64]; int k=factor(D,ps,es); int JinSq; int mu=mu_of(D,ps,es,k,&JinSq);
    ll g=1LL<<(mu-1); if(hp%g) die("genus count",D); if(pid==jid && !JinSq) die("J=1 but not in principal genus",D);
    if(!JinSq && pid==jid) die("x",D);
    int narrow2 = (hp==g); int wide2 = JinSq? (h==g) : (2*h==g);
    /* split-prime test of Cor. B (p with 4p^4 < D) */
    int split=0; for(ll p=2; 4*p*p*p*p<D; p++){ int isp=1; for(ll d=2;d*d<=p;d++) if(p%d==0){isp=0;break;} if(!isp) continue; if(kron(D,p)==1){ split=1; break; } }
    if(wide2 && D>=100 && split) die("Cor B violated",D);
    norders++; nw+=wide2; nn2+=narrow2;
    /* conductor and fundamental discriminant */
    ll f=1, u0=1; int e2=0; for(int i=0;i<k;i++){ if(ps[i]==2){ e2=es[i]; continue; } for(int e=0;e<es[i]/2;e++) f*=ps[i]; if(es[i]%2) u0*=ps[i]; }
    ll D0; if(e2==0) D0=u0; else if(e2%2==0){ if(u0%4==1){ D0=u0; for(int e=0;e<e2/2;e++) f*=2; } else { D0=4*u0; for(int e=0;e<(e2-2)/2;e++) f*=2; } } else { D0=8*u0; for(int e=0;e<(e2-3)/2;e++) f*=2; }
    if(f*f*D0!=D || !(D0%4==1 || (D0%4==0 && ((D0/4)%4==2 || (D0/4)%4==3)))) die("D0",D);
    if(wide2 || narrow2) printf("TORS N=%+d n=%lld D=%lld f=%lld D0=%lld h=%lld h+=%d mu=%d JinSq=%d wide2=%d narrow2=%d split=%d\n",N,n,D,f,D0,h,hp,mu,JinSq,wide2,narrow2,split);
    if(!split && n>1000) printf("SURVCLASS N=%+d n=%lld D=%lld h=%lld h+=%d wide2=%d narrow2=%d\n",N,n,D,h,hp,wide2,narrow2);
    if(n<=ANMAX || (!split && n>1000) || (SAMPLE && n>1000 && n%SAMPLE==0)){
      double L=Lchi(D0), pr=1; for(int i=0;i<k;i++) if(f%ps[i]==0) pr*= 1.0 - (double)kron(D0,ps[i])/ps[i];
      double eps=(n+sqrt((double)D))/2, R = (N==1&&n==3)? log((1+sqrt(5.0))/2) : log(eps);
      double ha=sqrt((double)D)*L*pr/(2*R), err=fabs(ha-h); if(err>maxerr) maxerr=err; nan++;
      if(err>1e-3) printf("ANALYTIC_MISMATCH N=%+d n=%lld D=%lld h=%lld ha=%.6f\n",N,n,D,h,ha);
    }
    free(pm); free(indeg); free(cy);
  }
  printf("CLASS NMAX=%lld orders=%lld wide2=%lld narrow2=%lld lemmaA_fail=%lld analytic_checked=%lld max|h_an-h|=%.3g\n",NMAX,norders,nw,nn2,lemA,nan,maxerr);
}

/* pcyc: Lemma A in Lagrange form at large n: principal rho-cycle of Delta = n^2 -/+ 4 for pseudo-random n in [10^6, 10^9]
   and the boundary n = 10^9, 10^9 - 1; checks every |a| < sqrt(Delta)/2 on it is 1, prints cycle lengths seen. */
static void pcyc(int K){
  ull x=88172645463325252ULL; ll worst=0; int lens[64]={0}; ll fails=0, cnt=0;
  for(int t=0;t<K+4;t++) for(int N=1;N>=-1;N-=2){
    ll n; if(t<K){ x^=x<<13; x^=x>>7; x^=x<<17; n=1000000+(ll)(x%999000001ULL); } else n = (t==K)?1000000000LL:(t==K+1?999999999LL:(t==K+2?1001LL:1000LL));
    ll D=n*n-4*N, s=isq(D), b0=((s-D)%2==0)?s:s-1; F f={1,b0,(b0*b0-D)/4}, g=f; int L=0;
    do{ ll A=g.a<0?-g.a:g.a; if(4*A*A<D && A!=1){ fails++; printf("LEMMA_A_FAIL n=%lld N=%d a=%lld\n",n,N,g.a); } F h; rho(g,D,s,&h); g=h; L++; if(L>1000000) die("cycle too long",n); }while(!(g.a==f.a&&g.b==f.b));
    if(L<64) lens[L]++; else worst++; cnt++;
  }
  printf("PCYC orders=%lld lemmaA_fail=%lld cycle-length histogram:",cnt,fails); for(int i=0;i<64;i++) if(lens[i]) printf(" %d:%d",i,lens[i]); printf(" (>=64: %lld)\n",worst);
}

int main(int argc, char **argv){
  if(argc>1 && !strcmp(argv[1],"sieve")) sieve();
  else if(argc>3 && !strcmp(argv[1],"class")){ if(argc>4) SAMPLE=atoll(argv[4]); classes(atoll(argv[2]), atoll(argv[3])); }
  else if(argc>2 && !strcmp(argv[1],"pcyc")) pcyc(atoi(argv[2]));
  else { fprintf(stderr,"usage\n"); return 1; }
  return 0; }

