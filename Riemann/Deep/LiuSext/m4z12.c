// m4z12.c — normalised fourth moment of sextic / cubic residue-symbol sums over F = Q(zeta_12) (unit rank 1).
// elements k = a0 + a1 z + a2 z^2 + a3 z^3, z = zeta_12, z^4 = z^2 - 1.  sigma1: z -> e^{i pi/6}, sigma5: z -> e^{5 i pi/6}.
// rows: one generator per ideal prime to 6 with N k <= U: |sigma1 k|/|sigma5 k| in [eps^{-1/2}, eps^{1/2}) (eps = 2+sqrt3 = ratio of 1+z),
//       root of unity fixed by the minimal residue code mod 3 among the 12 associates (mu_12 injects into (O/3)^x).
// columns: ideals n prime to 6, weight W(N n/X), X = U^m.  (k/p)_6 = k^{(Np-1)/6} read in mu_6 via zeta_6 = z^2.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <complex.h>
typedef long long ll;
static ll pm(ll a, ll e, ll m){ ll r=1; a%=m; if(a<0)a+=m; while(e){ if(e&1) r=r*a%m; a=a*a%m; e>>=1;} return r; }
static int factor(ll m, ll *f){ int nf=0; for(ll q=2;q*q<=m;q++) if(m%q==0){ f[nf++]=q; while(m%q==0) m/=q; } if(m>1) f[nf++]=m; return nf; }
static int *pr,npr; static void sieve(int n){ char*c=calloc(n+2,1); pr=malloc(sizeof(int)*(n/2+10)); npr=0; for(int i=2;i<=n;i++) if(!c[i]){ pr[npr++]=i; for(ll j=(ll)i*i;j<=n;j+=i) c[j]=1; } free(c); }
// F_{p^2} = F_p[t]/(t^2 - c)
static ll FC; static void fmul(ll p, ll a1,ll b1,ll a2,ll b2,ll*a,ll*b){ ll x=(a1*a2 + b1*b2%p*FC)%p, y=(a1*b2+a2*b1)%p; *a=x; *b=y; }
static void fpow(ll p, ll a,ll b,ll e,ll*ra,ll*rb){ ll xa=1,xb=0; while(e){ if(e&1) fmul(p,xa,xb,a,b,&xa,&xb); fmul(p,a,b,a,b,&a,&b); e>>=1;} *ra=xa;*rb=xb; }
typedef struct { int p, norm, type, u; ll off; ll r[4]; ll ra[4], rb[4]; } PI; // type0: z->r mod p (r[k]=z^k); type1: z->(ra,rb) in F_{p^2}, powers k=0..3
static PI *pis; static int npis; static uint8_t *tab; static ll tabn;
static int cmpPI(const void*x,const void*y){ return ((const PI*)x)->norm-((const PI*)y)->norm; }
static void build_primes(ll lim){ sieve((int)lim+1); pis=malloc(sizeof(PI)*(4*npr+10)); npis=0; tabn=0;
  for(int i=0;i<npr;i++){ int p=pr[i]; if(p<5) continue; if(p%12==1) tabn+=p; else if((ll)p*p<=lim) tabn+=(ll)p*p; }
  tab=malloc(tabn+1); ll off=0;
  for(int i=0;i<npr;i++){ ll p=pr[i]; if(p<5) continue;
    if(p%12==1){ ll f[32]; int nf=factor(p-1,f); ll g; for(g=2;;g++){ int ok=1; for(int j=0;j<nf;j++) if(pm(g,(p-1)/f[j],p)==1){ok=0;break;} if(ok)break; }
      ll x=1; for(ll t=0;t<p-1;t++){ tab[off+x]=(uint8_t)(t%6); x=x*g%p; } tab[off]=255;
      ll R=pm(g,(p-1)/12,p); ll z6=pm(g,(p-1)/6,p);
      int ex[4]={1,5,7,11}; for(int s=0;s<4;s++){ ll r=pm(R,ex[s],p); ll r2=r*r%p; int u=(r2==z6)?1:((r2==pm(z6,5,p))?5:-1); if(u<0){fprintf(stderr,"u fail %lld\n",p);exit(1);}
        PI q; memset(&q,0,sizeof q); q.p=(int)p; q.norm=(int)p; q.type=0; q.u=u; q.off=off; q.r[0]=1; for(int k=1;k<4;k++) q.r[k]=q.r[k-1]*r%p; pis[npis++]=q; }
      off+=p;
    } else if(p*p<=lim){ ll Q=p*p, o=Q-1; // nonresidue c
      for(FC=2;pm(FC,(p-1)/2,p)!=p-1;FC++); ll f[32]; int nf=factor(o,f); ll ga=0,gb=0; int fnd=0;
      for(ll ca=0;ca<p&&!fnd;ca++) for(ll cb=1;cb<p&&!fnd;cb++){ int ok=1; for(int j=0;j<nf;j++){ ll xa,xb; fpow(p,ca,cb,o/f[j],&xa,&xb); if(xa==1&&xb==0){ok=0;break;} } if(ok){ga=ca;gb=cb;fnd=1;} }
      ll xa=1,xb=0; for(ll t=0;t<o;t++){ tab[off+xa+xb*p]=(uint8_t)(t%6); fmul(p,xa,xb,ga,gb,&xa,&xb); } tab[off]=255;
      ll Ra,Rb; fpow(p,ga,gb,o/12,&Ra,&Rb); ll za,zb; fpow(p,ga,gb,o/6,&za,&zb); ll z5a,z5b; fpow(p,za,zb,5,&z5a,&z5b);
      // Frobenius orbits of primitive 12th roots R^e, e in {1,5,7,11}: pick e and e*p mod 12; one ideal per orbit
      int used[12]={0}; int ex[4]={1,5,7,11};
      for(int s=0;s<4;s++){ int e=ex[s]; if(used[e]) continue; used[e]=1; used[(e*(p%12))%12]=1;
        PI q; memset(&q,0,sizeof q); q.p=(int)p; q.norm=(int)Q; q.type=1; q.off=off; ll ra,rb; fpow(p,Ra,Rb,e,&ra,&rb);
        q.ra[0]=1;q.rb[0]=0; for(int k=1;k<4;k++) fmul(p,q.ra[k-1],q.rb[k-1],ra,rb,&q.ra[k],&q.rb[k]);
        ll r2a,r2b; fmul(p,ra,rb,ra,rb,&r2a,&r2b); q.u=(r2a==za&&r2b==zb)?1:((r2a==z5a&&r2b==z5b)?5:-1); if(q.u<0){fprintf(stderr,"u fail inert %lld\n",p);exit(1);}
        pis[npis++]=q; }
      // store FC per table: encode in r[0]
      pis[npis-1].r[0]=FC; pis[npis-2].r[0]=FC;
      off+=Q; }
  }
  qsort(pis,npis,sizeof(PI),cmpPI); }
static inline int jval(int i, const ll *a){ const PI*q=&pis[i]; ll p=q->p; uint8_t t;
  if(q->type==0){ ll x=0; for(int k=0;k<4;k++){ ll c=a[k]%p; if(c<0)c+=p; x=(x+c*q->r[k])%p; } t=tab[q->off+x]; }
  else { ll xa=0,xb=0; for(int k=0;k<4;k++){ ll c=a[k]%p; if(c<0)c+=p; xa=(xa+c*q->ra[k])%p; xb=(xb+c*q->rb[k])%p; } t=tab[q->off+xa+xb*p]; }
  if(t==255) return 255; return (t*q->u)%6; }
static int nid; static int *par,*pid; static double *nrm,*wt;
static void dfs(int cur,ll cn,int mi,ll lim){ for(int i=mi;i<npis;i++){ ll nn=cn*(ll)pis[i].norm; if(nn>lim) break; int me=nid++; par[me]=cur; pid[me]=i; nrm[me]=(double)nn; dfs(me,nn,i,lim);} }
static double Wf(double x){ if(x<=0.5||x>=2.0) return 0; double t=log2(x); return exp(1.0-1.0/(1.0-t*t)); }
static void build_ideals(double X){ ll lim=(ll)floor(2*X); int cap=(int)(lim+100); par=malloc(sizeof(int)*cap); pid=malloc(sizeof(int)*cap); nrm=malloc(sizeof(double)*cap); wt=malloc(sizeof(double)*cap);
  nid=1; par[0]=-1; pid[0]=-1; nrm[0]=1; dfs(0,1,0,lim); for(int i=0;i<nid;i++) wt[i]=Wf(nrm[i]/X); }
static uint8_t add7[64];
static void rowsum(const uint8_t *jp, uint8_t *e, double *bins){ static uint8_t j7[1<<16]; static double b7[8]; for(int j=0;j<8;j++) b7[j]=0; e[0]=0;
  for(int i=0;i<npis;i++) j7[i]=(jp[i]==255)?6:jp[i];
  for(int i=1;i<nid;i++){ uint8_t v=add7[(e[par[i]]<<3)|j7[pid[i]]]; e[i]=v; b7[v]+=wt[i]; } for(int j=0;j<6;j++) bins[j]=b7[j]; }
static const double C6[6]={1,0.5,-0.5,-1,-0.5,0.5}, S6[6]={0,0.86602540378443864676,0.86602540378443864676,0,-0.86602540378443864676,-0.86602540378443864676};
static void mulz(const ll*a,const ll*b,ll*c){ ll t[7]={0}; for(int i=0;i<4;i++) for(int j=0;j<4;j++) t[i+j]+=a[i]*b[j]; // reduce z^k, k>=4: z^4=z^2-1
  for(int k=6;k>=4;k--){ t[k-2]+=t[k]; t[k-4]-=t[k]; t[k]=0; } for(int i=0;i<4;i++) c[i]=t[i]; }
static double complex emb(const ll*a,double ang){ double complex z=cexp(I*ang), s=0, w=1; for(int k=0;k<4;k++){ s+=a[k]*w; w*=z; } return s; }
static int code3(const ll*a){ int c=0; for(int k=0;k<4;k++){ int r=(int)(((a[k]%3)+3)%3); c=c*3+r; } return c; }
int main(int argc,char**argv){ for(int x=0;x<8;x++) for(int y=0;y<8;y++) add7[(x<<3)|y]=(x>=6||y>=6)?6:(uint8_t)((x+y)%6);
  clock_t t0=clock(); const double eps=2+sqrt(3.0), lo=-0.5*log(eps), hi=0.5*log(eps);
  if(argc<3) return 1; double U=atof(argv[2]); double m=argc>3?atof(argv[3]):0.5; double X=pow(U,m);
  build_primes((ll)(2*X)); build_ideals(X); uint8_t *jp=malloc(npis), *e=malloc(nid); double bins[6];
  ll zeta[4]={0,1,0,0};
  if(!strcmp(argv[1],"test")){ // multiplicativity in k and agreement of jval with direct power k^{(q-1)/6}
    srand(3); int bad=0,tot=0; for(int it=0;it<20000;it++){ ll a[4],b[4],c[4]; for(int k=0;k<4;k++){a[k]=rand()%201-100;b[k]=rand()%201-100;} mulz(a,b,c); int i=rand()%npis;
      int j1=jval(i,a),j2=jval(i,b),j=jval(i,c); tot++; if(j1==255||j2==255){ if(j!=255) bad++; } else if(j!=(j1+j2)%6) bad++; }
    // norm check: |s1|^2|s5|^2 of random elements is an integer
    double maxdev=0; for(int it=0;it<1000;it++){ ll a[4]; for(int k=0;k<4;k++) a[k]=rand()%41-20; double n=pow(cabs(emb(a,M_PI/6))*cabs(emb(a,5*M_PI/6)),2); double d=fabs(n-llround(n)); if(d>maxdev) maxdev=d; }
    // zeta_6 = z^2 must have symbol exponent 1 at every prime ideal: (z^2/p)_6 = z^{2(q-1)/6}... check j(z^2)=(q-1)/6 mod 6
    int badz=0; ll z2[4]={0,0,1,0}; for(int i=0;i<npis;i++){ int j=jval(i,z2); if(j!=(int)(((pis[i].norm-1)/6)%6)) badz++; }
    printf("test X=%g primeideals=%d ideals=%d mult_bad=%d/%d norm_maxdev=%.2e zeta6_bad=%d\n",X,npis,nid,bad,tot,maxdev,badz); return 0; }
  if(!strcmp(argv[1],"rm")){ ll draws=atoll(argv[4]); unsigned long long st=88172645463325252ULL; long double m2=0,m46=0,m43=0,m86=0;
    for(ll d=0;d<draws;d++){ for(int i=0;i<npis;i++){ st^=st<<13; st^=st>>7; st^=st<<17; double u=(double)(st>>11)/9007199254740992.0; if(u<1.0/pis[i].norm) jp[i]=255; else jp[i]=(uint8_t)((st>>3)%6);} rowsum(jp,e,bins);
      double re=0,im=0,r3=0,i3=0; for(int j=0;j<6;j++){ re+=bins[j]*C6[j]; im+=bins[j]*S6[j]; r3+=bins[j]*C6[(2*j)%6]; i3+=bins[j]*S6[(2*j)%6]; }
      double a6=(re*re+im*im)/X, a3=(r3*r3+i3*i3)/X; m2+=a6; m46+=a6*a6; m43+=a3*a3; m86+=a6*a6*a6*a6; }
    printf("z12rm U=%.0f X=%g draws=%lld E2=%.6Lf E4_6=%.6Lf E4_3=%.6Lf E8_6=%.6Lf time=%.1f\n",U,X,draws,m2/draws,m46/draws,m43/draws,m86/draws,(double)(clock()-t0)/CLOCKS_PER_SEC); return 0; }
  // rows
  double R2=sqrt(eps)*sqrt(U); // max |sigma_i|^2 in the domain: |s1|^2 = rho*|s1||s5| <= sqrt(eps)*sqrt(U)
  ll Bx=(ll)ceil(sqrt(4.0*R2/3.0))+2; // x = x1 + x2 w with |x|^2 <= R2
  long double M2_6=0,M4_6=0,M2_3=0,M4_3=0; ll rows=0,pr6=0,pr3=0,cand=0; double mx6=0;
  for(ll x1=-Bx;x1<=Bx;x1++) for(ll x2=-Bx;x2<=Bx;x2++){ double nx=(double)(x1*x1-x1*x2+x2*x2); if(nx>R2) continue;
    for(ll y1=-Bx;y1<=Bx;y1++){ for(ll y2=-Bx;y2<=Bx;y2++){ double ny=(double)(y1*y1-y1*y2+y2*y2); if(nx+ny>R2) continue;
      ll a[4]={x1-x2,-y2,x2,y1}; double complex s1=emb(a,M_PI/6), s5=emb(a,5*M_PI/6); double n1=creal(s1*conj(s1)), n5=creal(s5*conj(s5));
      if(n1>R2*1.0000001||n5>R2*1.0000001) continue; double N=n1*n5; if(N<0.5||N>U*(1+1e-12)) continue; ll Ni=llround(N); if(Ni>U) continue; if(Ni%2==0||Ni%3==0) continue;
      double lr=0.5*log(n1/n5); if(!(lr>=lo-1e-12 && lr<hi-1e-12)) continue; cand++;
      // root-of-unity normalisation: keep iff code3 minimal among 12 associates
      int c0=code3(a), ok=1; ll b[4]; memcpy(b,a,sizeof b); for(int j=1;j<12;j++){ ll t[4]; mulz(b,zeta,t); memcpy(b,t,sizeof b); if(code3(b)<c0){ok=0;break;} if(code3(b)==c0){ fprintf(stderr,"code tie\n"); exit(1);} }
      if(!ok) continue; rows++;
      int p6=1,p3=1; for(int i=0;i<npis;i++){ int j=jval(i,a); jp[i]=(uint8_t)j; if(j!=255){ if(j) p6=0; if(j%3) p3=0; } }
      if(p6&&p3){ pr6++; pr3++; continue; } rowsum(jp,e,bins);
      double re=0,im=0,r3=0,i3=0; for(int j=0;j<6;j++){ re+=bins[j]*C6[j]; im+=bins[j]*S6[j]; r3+=bins[j]*C6[(2*j)%6]; i3+=bins[j]*S6[(2*j)%6]; }
      if(!p6){ double v=(re*re+im*im)/X; M2_6+=v; M4_6+=v*v; if(v>mx6) mx6=v; } else pr6++;
      if(!p3){ double v=(r3*r3+i3*i3)/X; M2_3+=v; M4_3+=v*v; } else pr3++; } } }
  printf("z12 U=%.0f m=%g X=%.3f ideals=%d primeideals=%d cand=%lld rows=%lld principal6=%lld principal3=%lld M2_6=%.6Le M4_6=%.6Le M2_3=%.6Le M4_3=%.6Le max6=%.4f time=%.1f\n",
    U,m,X,nid,npis,cand,rows,pr6,pr3,M2_6,M4_6,M2_3,M4_3,mx6,(double)(clock()-t0)/CLOCKS_PER_SEC); return 0; }
