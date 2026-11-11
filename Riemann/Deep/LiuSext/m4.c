// m4.c — normalised fourth moment of sextic / cubic residue-symbol sums over Q(sqrt(-3)).
// rows: k = a + b*w (w = e^{2 pi i/3}), k = 1 mod 3, N k odd, N k <= U (one generator per ideal prime to 6).
// columns: ideals n prime to 6, smooth weight W(N n / X), X = U^m (default m = 1/2).
// S6_k = sum_n (k/n)_6 W(Nn/X), S3_k = sum_n (k/n)_3 W(Nn/X)  [(k/n)_3 = (k/n)_6^2].
// M4 = sum over nonprincipal rows |S|^4/X^2, M2 = sum |S|^2/X.  Modes: run U m | rm X draws seed | test X
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
typedef long long ll;
static ll pm(ll a, ll e, ll m){ ll r=1; a%=m; if(a<0)a+=m; while(e){ if(e&1) r=r*a%m; a=a*a%m; e>>=1;} return r; }
static int *pr; static int npr;
static void sieve(int n){ char *c=calloc(n+2,1); pr=malloc(sizeof(int)*(n/2+10)); npr=0;
  for(int i=2;i<=n;i++) if(!c[i]){ pr[npr++]=i; for(ll j=(ll)i*i;j<=n;j+=i) c[j]=1; } free(c); }
static int factor(ll m, ll *f){ int nf=0; for(ll q=2;q*q<=m;q++) if(m%q==0){ f[nf++]=q; while(m%q==0) m/=q; } if(m>1) f[nf++]=m; return nf; }
typedef struct { int P, norm, type, r, u; ll off; } PI; // type 0: split, k -> a+b r mod P ; type 1: inert, k -> (a mod P)+(b mod P)P
static PI *pis; static int npis; static uint8_t *tab; static ll tabn;
static void f2mul(ll P, ll a1, ll b1, ll a2, ll b2, ll *a, ll *b){ // (a1+b1 w)(a2+b2 w), w^2=-1-w
  ll x=(a1*a2 - b1*b2)%P; ll y=(a1*b2 + a2*b1 - b1*b2)%P; if(x<0)x+=P; if(y<0)y+=P; *a=x; *b=y; }
static void f2pow(ll P, ll a, ll b, ll e, ll *ra, ll *rb){ ll xa=1,xb=0; while(e){ if(e&1) f2mul(P,xa,xb,a,b,&xa,&xb); f2mul(P,a,b,a,b,&a,&b); e>>=1; } *ra=xa; *rb=xb; }
static int cmpPI(const void*x,const void*y){ const PI*a=x,*b=y; return a->norm-b->norm; }
static void build_primes(ll lim){ // prime ideals of norm <= lim, prime to 6
  sieve((int)lim+1); pis=malloc(sizeof(PI)*(2*npr+10)); npis=0; tabn=0;
  for(int i=0;i<npr;i++){ int P=pr[i]; if(P<5) continue;
    if(P%3==1) tabn+=P; else if((ll)P*P<=lim) tabn+=(ll)P*P; }
  tab=malloc(tabn+1); ll off=0;
  for(int i=0;i<npr;i++){ int P=pr[i]; if(P<5) continue;
    if(P%3==1){ ll f[32]; int nf=factor(P-1,f); int g; for(g=2;;g++){ int ok=1; for(int j=0;j<nf;j++) if(pm(g,(P-1)/f[j],P)==1){ok=0;break;} if(ok)break; }
      ll x=1; for(ll t=0;t<P-1;t++){ tab[off+x]=(uint8_t)(t%6); x=x*g%P; } tab[off]=255;
      ll z=pm(g,(P-1)/6,P); ll r=pm(g,(P-1)/3,P); // r: root of x^2+x+1
      for(int s=0;s<2;s++){ ll rr = s? (P-1-r) : r; ll one=(1+rr)%P; int u = (one==z)?1:((one==pm(z,5,P))?5:-1); if(u<0){fprintf(stderr,"u fail %d\n",P); exit(1);}
        PI q={P,P,0,(int)rr,u,off}; pis[npis++]=q; }
      off+=P;
    } else if((ll)P*P<=lim){ ll Q=(ll)P*P, o=Q-1; ll f[32]; int nf=factor(o,f); ll ga=0,gb=0; int found=0;
      for(ll ca=0;ca<P && !found;ca++) for(ll cb=1;cb<P && !found;cb++){ int ok=1; for(int j=0;j<nf;j++){ ll xa,xb; f2pow(P,ca,cb,o/f[j],&xa,&xb); if(xa==1&&xb==0){ok=0;break;} } if(ok){ga=ca;gb=cb;found=1;} }
      ll xa=1,xb=0; for(ll t=0;t<o;t++){ tab[off+xa+xb*P]=(uint8_t)(t%6); f2mul(P,xa,xb,ga,gb,&xa,&xb); } tab[off]=255;
      ll za,zb; f2pow(P,ga,gb,o/6,&za,&zb); ll z5a,z5b; f2pow(P,za,zb,5,&z5a,&z5b);
      int u = (za==1&&zb==1)?1:((z5a==1&&z5b==1)?5:-1); if(u<0){fprintf(stderr,"u fail inert %d\n",P); exit(1);}
      PI q={P,(int)Q,1,0,u,off}; pis[npis++]=q; off+=Q; }
  }
  qsort(pis,npis,sizeof(PI),cmpPI);
}
static inline int jval(int i, ll a, ll b){ const PI*q=&pis[i]; ll P=q->P; ll x;
  if(q->type==0){ x=(a%P + (b%P)*q->r)%P; if(x<0)x+=P; }
  else { ll aa=a%P, bb=b%P; if(aa<0)aa+=P; if(bb<0)bb+=P; x=aa+bb*P; }
  uint8_t t=tab[q->off+x]; if(t==255) return 255; return (t*q->u)%6; }
static int nid; static int *par,*pid; static double *nrm,*wt;
static void dfs(int cur, ll cn, int mi, ll lim){ for(int i=mi;i<npis;i++){ ll nn=cn*(ll)pis[i].norm; if(nn>lim) break; int me=nid++; par[me]=cur; pid[me]=i; nrm[me]=(double)nn; dfs(me,nn,i,lim);} }
static double Wf(double x){ if(x<=0.5||x>=2.0) return 0; double t=log2(x); return exp(1.0-1.0/(1.0-t*t)); }
static void build_ideals(double X){ ll lim=(ll)floor(2.0*X); int cap=(int)(lim+100); par=malloc(sizeof(int)*cap); pid=malloc(sizeof(int)*cap); nrm=malloc(sizeof(double)*cap); wt=malloc(sizeof(double)*cap);
  nid=1; par[0]=-1; pid[0]=-1; nrm[0]=1; dfs(0,1,0,lim);
  // reorder by norm (parents strictly earlier, far away: no store-load chains)
  int *ord=malloc(sizeof(int)*nid), *inv=malloc(sizeof(int)*nid); for(int i=0;i<nid;i++) ord[i]=i;
  for(int i=1;i<nid;i++){ int x=ord[i]; int j=i-1; while(j>=0 && nrm[ord[j]]>nrm[x]){ ord[j+1]=ord[j]; j--; } ord[j+1]=x; } // insertion sort (nid small)
  for(int i=0;i<nid;i++) inv[ord[i]]=i;
  int *p2=malloc(sizeof(int)*nid), *q2=malloc(sizeof(int)*nid); double *n2=malloc(sizeof(double)*nid);
  for(int i=0;i<nid;i++){ int o=ord[i]; p2[i]= (par[o]<0)?-1:inv[par[o]]; q2[i]=pid[o]; n2[i]=nrm[o]; }
  for(int i=0;i<nid;i++){ par[i]=p2[i]; pid[i]=q2[i]; nrm[i]=n2[i]; if(i>0 && par[i]>=i){ fprintf(stderr,"order fail\n"); exit(1);} }
  for(int i=0;i<nid;i++) wt[i]=Wf(nrm[i]/X); }
static uint8_t addm[256][256];
static const double C6[6]={1,0.5,-0.5,-1,-0.5,0.5}, S6[6]={0,0.86602540378443864676,0.86602540378443864676,0,-0.86602540378443864676,-0.86602540378443864676};
static uint8_t add7[64];
static void rowsum(const uint8_t *jp, uint8_t *e, double *bins){ // jp: 0..5 or 255 (zero); internally 6 = zero
  static uint8_t j7[1<<16]; static double b7[8]; for(int j=0;j<8;j++) b7[j]=0; e[0]=0;
  for(int i=0;i<npis;i++) j7[i]=(jp[i]==255)?6:jp[i];
  for(int i=1;i<nid;i++){ uint8_t v=add7[(e[par[i]]<<3)|j7[pid[i]]]; e[i]=v; b7[v]+=wt[i]; }
  for(int j=0;j<6;j++) bins[j]=b7[j]; }
int main(int argc,char**argv){
  for(int x=0;x<256;x++) for(int y=0;y<256;y++) addm[x][y]=(x==255||y==255)?255:(uint8_t)((x+y)%6);
  for(int x=0;x<8;x++) for(int y=0;y<8;y++) add7[(x<<3)|y]=(x>=6||y>=6)?6:(uint8_t)((x+y)%6);
  if(argc<2){ fprintf(stderr,"usage\n"); return 1; }
  clock_t t0=clock();
  if(!strcmp(argv[1],"test")){ double X=atof(argv[2]); build_primes((ll)(2*X)); srand(7); int bad=0,tot=0,badm=0,totm=0;
    // multiplicativity in k at every prime ideal
    for(int it=0;it<20000;it++){ ll a1=rand()%2001-1000,b1=rand()%2001-1000,a2=rand()%2001-1000,b2=rand()%2001-1000; ll a=a1*a2-b1*b2, b=a1*b2+a2*b1-b1*b2; int i=rand()%npis;
      int j1=jval(i,a1,b1), j2=jval(i,a2,b2), j=jval(i,a,b); if(j1==255||j2==255){ if(j!=255) badm++; } else if(j!=(j1+j2)%6) badm++; totm++; }
    // cubic reciprocity for prime elements k = 1 mod 3 of prime norm (split)
    for(int i1=0;i1<npis;i1++) for(int i2=0;i2<npis;i2++){ if(pis[i1].type||pis[i2].type) continue; if(pis[i1].P==pis[i2].P) continue; if(tot>200000) break;
      // element generating ideal (P, w - r): find a + b w with norm P, a+b r = 0 mod P, normalised to 1 mod 3
      ll g[2][2]; int ok=1; for(int s=0;s<2;s++){ const PI*q=&pis[s?i2:i1]; ll P=q->P; int fnd=0; ll B=(ll)sqrt(4.0*P/3)+2;
        for(ll b=-B;b<=B&&!fnd;b++){ ll disc=4*P-3*b*b; if(disc<0) continue; ll sq=(ll)llround(sqrt((double)disc)); for(ll dd=sq-1;dd<=sq+1&&!fnd;dd++){ if(dd<0||dd*dd!=disc) continue; for(int sg=-1;sg<=1;sg+=2){ ll num=b+sg*dd; if(num%2) continue; ll a=num/2; if(a*a-a*b+b*b!=P) continue; if(((a%P+(b%P)*q->r)%P+P)%P!=0) continue;
          // normalise: multiply by units until = 1 mod 3: units 1,w,w^2 and negatives; a+bw times w = -b + (a-b)w
          ll ca=a,cb=b; for(int k=0;k<6;k++){ if((((ca-1)%3)+3)%3==0 && ((cb%3)+3)%3==0){ g[s][0]=ca; g[s][1]=cb; fnd=1; break;} ll na=-cb, nb=ca-cb; ca=na; cb=nb; if(k==2){ca=-ca;cb=-cb;} }
          if(fnd) break; } } }
        if(!fnd) ok=0; }
      if(!ok) continue; int x=jval(i1,g[1][0],g[1][1]), y=jval(i2,g[0][0],g[0][1]); if(x==255||y==255) continue; tot++; if(x%3!=y%3) bad++; }
    printf("test X=%g primeideals=%d mult_bad=%d/%d cubic_recip_bad=%d/%d\n",X,npis,badm,totm,bad,tot); return 0; }
  if(!strcmp(argv[1],"rm")){ double X=atof(argv[2]); ll draws=atoll(argv[3]); unsigned long long st=strtoull(argv[4],0,10)*2654435761ULL+1; build_primes((ll)(2*X)); build_ideals(X);
    uint8_t *jp=malloc(npis), *e=malloc(nid); double bins[6]; long double m2=0,m4_6=0,m4_3=0,m8_6=0;
    for(ll d=0;d<draws;d++){ for(int i=0;i<npis;i++){ st^=st<<13; st^=st>>7; st^=st<<17; double u=(double)(st>>11)/9007199254740992.0; if(u<1.0/pis[i].norm) jp[i]=255; else jp[i]=(uint8_t)((st>>3)%6);} rowsum(jp,e,bins);
      double re=0,im=0,r3=0,i3=0; for(int j=0;j<6;j++){ re+=bins[j]*C6[j]; im+=bins[j]*S6[j]; r3+=bins[j]*C6[(2*j)%6]; i3+=bins[j]*S6[(2*j)%6]; }
      double a6=(re*re+im*im)/X, a3=(r3*r3+i3*i3)/X; m2+=a6; m4_6+=a6*a6; m4_3+=a3*a3; m8_6+=a6*a6*a6*a6; }
    printf("rm X=%g ideals=%d primeideals=%d draws=%lld E2=%.6Lf E4_6=%.6Lf E4_3=%.6Lf E8_6=%.6Lf time=%.1f\n",X,nid,npis,draws,m2/draws,m4_6/draws,m4_3/draws,m8_6/draws,(double)(clock()-t0)/CLOCKS_PER_SEC); return 0; }
  if(!strcmp(argv[1],"prin")){ // principal rows: k = j^6 (sextic) and k = j^3 (cubic), j = 1 mod 3, N j odd; their excluded contribution
    double U=atof(argv[2]); double m=argc>3?atof(argv[3]):0.5; double X=pow(U,m); build_primes((ll)(2*X)); build_ideals(X);
    uint8_t *jp=malloc(npis), *e=malloc(nid); double bins[6]; long double P6=0,P3=0; ll n6=0,n3=0; ll B=(ll)sqrt(4.0*cbrt(U)/3.0)+3;
    for(ll b=-B;b<=B;b++) for(ll a=-B-2;a<=B+2;a++){ if((((a-1)%3)+3)%3||((b%3)+3)%3) continue; ll N=a*a-a*b+b*b; if(N<=0||!(N&1)) continue; if((double)N*N*N>U) continue;
      ll ca=1,cb=0; for(int t=0;t<3;t++){ ll na=ca*a-cb*b, nb=ca*b+a*cb-cb*b; ca=na; cb=nb; } // j^3
      for(int pass=0;pass<2;pass++){ ll ka=ca,kb=cb; if(pass==1){ if((double)N*N*N*N*N*N>U) break; ll na=ca*ca-cb*cb, nb=2*ca*cb-cb*cb; ka=na; kb=nb; } // j^6
        for(int i=0;i<npis;i++) jp[i]=(uint8_t)jval(i,ka,kb); rowsum(jp,e,bins);
        double re=0,im=0,r3=0,i3=0; for(int j=0;j<6;j++){ re+=bins[j]*C6[j]; im+=bins[j]*S6[j]; r3+=bins[j]*C6[(2*j)%6]; i3+=bins[j]*S6[(2*j)%6]; }
        if(pass==0){ double a3=(r3*r3+i3*i3)/X; P3+=a3*a3; n3++; } else { double a6=(re*re+im*im)/X; P6+=a6*a6; n6++; } } }
    printf("prin U=%.0f m=%g X=%.3f sixth_powers=%lld P4_6=%.6Le cubes=%lld P4_3=%.6Le time=%.1f\n",U,m,X,n6,P6,n3,P3,(double)(clock()-t0)/CLOCKS_PER_SEC); return 0; }
  if(!strcmp(argv[1],"run")){ double U=atof(argv[2]); double m=argc>3?atof(argv[3]):0.5; double X=pow(U,m); build_primes((ll)(2*X)); build_ideals(X);
    uint8_t *jp=malloc(npis), *e=malloc(nid); double bins[6]; long double M2_6=0,M4_6=0,M2_3=0,M4_3=0; ll rows=0,pr6=0,pr3=0; double mx6=0,mx3=0; ll mxa6=0,mxb6=0,mxa3=0,mxb3=0;
    ll bmax=(ll)sqrt(4.0*U/3.0)+3;
    // dyadic shells by N k for per-shell M4 (U/2^j < N k <= U/2^{j-1})
    long double sh6[64]={0}, sh3[64]={0}; ll shn[64]={0};
    int *cur=malloc(sizeof(int)*npis), *bb=malloc(sizeof(int)*npis);
    // conjugation k -> conj(k) maps rows to rows (b -> -b) and S_k -> conj(S_k): compute b >= 0, weight 2 for b > 0
    for(ll b=0;b<=bmax;b++){ if(((b%3)+3)%3) continue; const double cw=(b==0)?1.0:2.0; double disc=U-0.75*(double)b*b; if(disc<0) continue; double s=sqrt(disc); ll alo=(ll)ceil(b/2.0-s), ahi=(ll)floor(b/2.0+s);
      while((((alo-1)%3)+3)%3) alo++;
      for(int i=0;i<npis;i++){ ll P=pis[i].P; if(pis[i].type==0){ ll x=(alo%P+(b%P)*pis[i].r)%P; if(x<0)x+=P; cur[i]=(int)x; } else { ll x=alo%P; if(x<0)x+=P; cur[i]=(int)x; ll y=b%P; if(y<0)y+=P; bb[i]=(int)(y*P); } }
      for(ll a=alo;a<=ahi;a+=3){
        int p6=1,p3=1; for(int i=0;i<npis;i++){ const PI*q=&pis[i]; int x=cur[i]; uint8_t t=tab[q->off + x + (q->type?bb[i]:0)]; int j = (t==255)?255:(q->u==1? t : (6-t)%6); jp[i]=(uint8_t)j; if(j!=255){ if(j) p6=0; if(j%3) p3=0; } x+=3; if(x>=q->P) x-=q->P; cur[i]=x; }
        ll N=a*a-a*b+b*b; if(N>U||N<=0) continue; if(!(N&1)) continue;
        rows+=(ll)cw; if(p6&&p3){pr6+=(ll)cw;pr3+=(ll)cw;continue;} rowsum(jp,e,bins);
        double re=0,im=0,r3=0,i3=0; for(int j=0;j<6;j++){ re+=bins[j]*C6[j]; im+=bins[j]*S6[j]; r3+=bins[j]*C6[(2*j)%6]; i3+=bins[j]*S6[(2*j)%6]; }
        int sh=(int)floor(log2(U/(double)N)); if(sh>63) sh=63;
        if(!p6){ double a6=(re*re+im*im)/X; M2_6+=cw*a6; M4_6+=cw*a6*a6; sh6[sh]+=cw*a6*a6; if(a6>mx6){mx6=a6;mxa6=a;mxb6=b;} } else pr6+=(ll)cw;
        if(!p3){ double a3=(r3*r3+i3*i3)/X; M2_3+=cw*a3; M4_3+=cw*a3*a3; sh3[sh]+=cw*a3*a3; if(a3>mx3){mx3=a3;mxa3=a;mxb3=b;} } else pr3+=(ll)cw;
        shn[sh]+=(ll)cw; } }
    printf("run U=%.0f m=%g X=%.3f ideals=%d primeideals=%d rows=%lld principal6=%lld principal3=%lld M2_6=%.6Le M4_6=%.6Le M2_3=%.6Le M4_3=%.6Le max6=%.4f@(%lld,%lld) max3=%.4f@(%lld,%lld) time=%.1f\n",
      U,m,X,nid,npis,rows,pr6,pr3,M2_6,M4_6,M2_3,M4_3,mx6,mxa6,mxb6,mx3,mxa3,mxb3,(double)(clock()-t0)/CLOCKS_PER_SEC);
    printf("shells(j: U/2^(j+1)<Nk<=U/2^j, rows, M4_6, M4_3):"); for(int j=0;j<6;j++) printf(" [%d %lld %.4Le %.4Le]",j,shn[j],sh6[j],sh3[j]); printf("\n"); return 0; }
  return 1; }
