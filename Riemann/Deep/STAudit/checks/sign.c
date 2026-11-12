// sign.c (staudit): for fundamental D<0, |D| in [Dlo,Dhi], every class C of Q(sqrt D):
// least degree-one prime p (p not | 11 D) with a 'prime above p' in C and a_E(p)>0, resp. <0, resp. > -sqrt(p)/2 (E=11a1).
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
typedef long long ll;
static int np; static int *P, *AP;
static ll pw(ll b, ll e, ll m){ ll r=1; b%=m; if(b<0)b+=m; while(e){ if(e&1) r=(__int128)r*b%m; b=(__int128)b*b%m; e>>=1;} return r; }
static ll sqrtmod(ll a, ll p){ a%=p; if(a<0)a+=p; if(a==0) return 0; if(p%4==3) return pw(a,(p+1)/4,p);
  ll q=p-1; int s=0; while(q%2==0){q/=2;s++;} ll z=2; while(pw(z,(p-1)/2,p)!=p-1) z++;
  ll M=s, c=pw(z,q,p), t=pw(a,q,p), R=pw(a,(q+1)/2,p);
  while(t!=1){ int i=0; ll tt=t; while(tt!=1){ tt=(__int128)tt*tt%p; i++; } ll b=c; for(int j=0;j<M-i-1;j++) b=(__int128)b*b%p; M=i; c=(__int128)b*b%p; t=(__int128)t*c%p; R=(__int128)R*b%p; }
  return R; }
static void reduce(ll *a, ll *b, ll *c){ for(;;){ // normalize b in (-a,a]
    if(*b> *a || *b<= -*a){ ll k = (*b>0)? (*b + *a - 1)/(2* *a) : -((-*b + *a)/(2* *a)); // b-2ka in (-a,a]
      ll nb=*b-2*k* *a; // c' = (nb^2 - D)/4a computed from old: c - k b + k^2 a
      *c = *c - k* *b + k*k* *a; *b=nb; if(*b<= -*a){ *b+=2* *a; *c += *b - *a; } }
    if(*c < *a){ ll t=*a; *a=*c; *c=t; *b=-*b; continue; }
    if(*c==*a && *b<0) *b=-*b; break; } }
static int sqfree(ll m){ for(ll p=2;p*p<=m;p++) if(m%(p*p)==0) return 0; return 1; }
static int isfund(ll m){ ll D=-m; if(((D%4)+4)%4==1) return sqfree(m); if(m%4==0){ ll k=m/4, kk=((-k)%4+4)%4; return (kk==2||kk==3)&&sqfree(k);} return 0; }
#define AMAX 400
static int idx[AMAX+1][2*AMAX+2];
int main(int argc,char**argv){
  ll lo=atoll(argv[1]), hi=atoll(argv[2]); double A0=atof(argv[3]);
  P=malloc(200000*sizeof(int)); AP=malloc(200000*sizeof(int)); FILE*f=fopen("ap.txt","r"); np=0; while(fscanf(f,"%d %d",&P[np],&AP[np])==2) np++; fclose(f);
  for(int i=0;i<=AMAX;i++) for(int j=0;j<2*AMAX+2;j++) idx[i][j]=-1;
  ll *fa=malloc(20000*sizeof(ll)), *fb=malloc(20000*sizeof(ll)); ll *mp=malloc(20000*sizeof(ll)), *mn=malloc(20000*sizeof(ll)), *mh=malloc(20000*sizeof(ll));
  double emaxp=0,emaxn=0,emaxh=0; ll Dp=0,Dn=0,Dh=0; long nD=0, fail=0, failh=0, incomplete=0; double rworst=1e9; ll Dr=0;
  for(ll m=lo;m<=hi;m++){ if(m<=4||!isfund(m)) continue; ll D=-m; int h=0; ll am=(ll)sqrt(m/3.0)+1; if(am>AMAX){fprintf(stderr,"AMAX\n");return 1;}
    for(ll a=1;a<=am;a++) for(ll b=-a+1;b<=a;b++){ if(((b-D)%2+2)%2) continue; ll num=b*b-D; if(num%(4*a)) continue; ll c=num/(4*a); if(c<a) continue; if(c==a&&b<0) continue;
        ll g=a; ll x=llabs(b); while(x){ll t=g%x;g=x;x=t;} x=c; ll g2=g; while(x){ll t=g2%x;g2=x;x=t;} if(g2!=1) continue;
        idx[a][b+AMAX]=h; fa[h]=a; fb[h]=b; h++; }
    for(int i=0;i<h;i++){ mp[i]=mn[i]=mh[i]=0; } int need=3*h, k;
    for(k=0;k<np && need>0;k++){ ll p=P[k]; int ap=AP[k]; if(p==11 || m%p==0) continue;
      ll a,b,c;
      if(p==2){ if(((D%8)+8)%8!=1) continue; a=2; b=1; c=(1-D)/8; }
      else { if(pw(((D%p)+p)%p,(p-1)/2,p)!=1) continue; ll r=sqrtmod(D,p); if(((r-D)%2+2)%2) r+=p; a=p; b=r; c=(b*b-D)/(4*p); }
      ll a2=a,b2=-b,c2=c; reduce(&a,&b,&c); reduce(&a2,&b2,&c2);
      int i1=idx[a][b+AMAX], i2=idx[a2][b2+AMAX]; if(i1<0||i2<0){ fprintf(stderr,"lookup fail D=%lld p=%lld (%lld,%lld,%lld)\n",D,p,a,b,c); return 1; }
      int ii[2]={i1,i2};
      for(int t=0;t<2;t++){ int i=ii[t];
        if(ap>0 && !mp[i]){ mp[i]=p; need--; } if(ap<0 && !mn[i]){ mn[i]=p; need--; }
        if(ap > -0.5*sqrt((double)p) && !mh[i]){ mh[i]=p; need--; } }
    }
    for(int i=0;i<h;i++) idx[fa[i]][fb[i]+AMAX]=-1;
    if(need>0){ incomplete++; continue; }
    ll Mp=0,Mn=0,Mh=0; for(int i=0;i<h;i++){ if(mp[i]>Mp)Mp=mp[i]; if(mn[i]>Mn)Mn=mn[i]; if(mh[i]>Mh)Mh=mh[i]; }
    double L=log((double)m), ep=log((double)Mp)/L, en=log((double)Mn)/L, eh=log((double)Mh)/L; nD++;
    if(ep>emaxp){emaxp=ep;Dp=D;} if(en>emaxn){emaxn=en;Dn=D;} if(eh>emaxh){emaxh=eh;Dh=D;}
    double bound=2*pow((double)m,A0); double r=bound/(double)(Mp>Mn?Mp:Mn); if(r<rworst){rworst=r;Dr=D;}
    if(Mp>bound||Mn>bound) fail++; if(Mh>2*pow((double)m,1.0)) failh++;
  }
  printf("|D| in [%lld,%lld]: nD=%ld incomplete(ap table too short)=%ld | max_D log(p+max)/log|D| = %.3f (D=%lld), log(p-max)/log|D| = %.3f (D=%lld), a_p>-sqrt(p)/2: %.3f (D=%lld) | fails of p+-max <= 2|D|^%.4f: %ld; tightest ratio 2|D|^A/max(p+,p-) = %.2f at D=%lld | fails of p_(>-sqrt p/2) max <= 2|D|: %ld\n",
    lo,hi,nD,incomplete,emaxp,Dp,emaxn,Dn,emaxh,Dh,A0,fail,rworst,Dr,failh);
  return 0; }
