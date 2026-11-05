// classify_sumset.c -- exhaustive test of the Wang-Sun bounds (L), (C) (EJC 33(2) #P2.64, Conj. 1; OEIS A357130 Conj. 2).
// Build: cc -O2 -o classify_sumset classify_sumset.c.  Usage: ./classify_sumset m1 m2 n umode [eqsize] [pcap] [show] [fixA1]
// Expected (2026-11-05), "L/C violating tuples":  Z/5 n=3: 20/0, n=4: 0/0, n=5: 470/0, n=6: 0/0;  Z/7 n=3: 203/0,
//   n=4: 0/336, n=5: 9471/0;  Z/11 n=3: 495/0;  Z/13 n=3: 858/0;  Z/11 n=4 with A_1={0,1}: 0/0;
//   Z ([0,11), n=3): only ({0,1,2},{0,2},{0,1,2}) up to affine maps;  Z ([0,7), n=4): 0/0;  Z ([0,6), n=5): 492/0;
//   Z/4, Z/2^2 (n<=6), Z/6 (n<=5), Z/8, Z/2xZ/4, Z/9, Z/3^2 (n<=4): 0/0;  equal sizes (Z/5 n<=8, Z/7 n<=7, Z/11 n=3,4): 0/0.
// Exhaustive test of Sun's bounds (L) and (C) (Wang-Sun, EJC 33(2) #P2.64, Conj. 1; OEIS A357130 Conj. 2).
// G = Z/m1 x Z/m2 (N = m1*m2 <= 64), sets A_i taken from a universe U of u <= 12 elements of G.
// L = {a_1+..+a_n : a_i in A_i, a_i != a_{i+1} (1<=i<n)},  C = L with also a_n != a_1.
// (L): |L| >= min{p(G), S-2n+1+(n mod 2)},  (C): |C| >= min{p(G), S-2n+(-1)^n(1+(n mod 2))},  S = sum |A_i|.
// usage: sumset m1 m2 n umode [eqsize] [pcap]
//   umode = 0: U = G (N <= 12; A_1 reduced mod x -> cx+v, weights = orbit sizes);  umode = M > 0: U = {0..M-1} in Z/m1
//   eqsize = 0 (all sizes) or k (all |A_i| = k);  pcap = -1: no p(G) cap (simulate Z with m1 > n(M-1))
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
typedef uint64_t u64;
static int m1,m2,N,n,u,eqk; static int ad[64][64]; static int U[16];
static int pG;
static u64 shtab[64][8][256];
static inline u64 shift(u64 x, int b){ u64 r=0; for(int j=0;j<8 && (x>>(8*j));j++) r|=shtab[b][j][(x>>(8*j))&255]; return r; }
static int nsub; static int subm[65536], subsz[65536];
static long long cntL=0,cntC=0,tot=0; static int minSL=999,minSC=999; static int nshowL=0,nshowC=0;
static int choice[8]; static long long wt[65536]; static int isrep[65536]; static long long curw=1;
static long long hist[2][64];  // violations by S
static u64 S[8][16];
static u64 SC[8][16][16];
static int boundL(int s){ int b=s-2*n+1+(n%2); return (pG>0 && pG<b)? pG : b; }
static int boundC(int s){ int b=s-2*n+((n%2)?-1:1)*(1+(n%2)); return (pG>0 && pG<b)? pG : b; }
static int SHOW=6;
static void printex(const char*tag,int s,int sz){
  printf("%s S=%d |set|=%d bound=%d :", tag, s, sz, tag[0]=='L'?boundL(s):boundC(s));
  for(int i=0;i<n;i++){ printf(" {"); int m=subm[choice[i]],first=1; for(int a=0;a<u;a++) if(m>>a&1){ printf(first?"%d":",%d",U[a]); first=0;} printf("}"); }
  printf("\n");
}
static void rec(int i, int s){
  if(i==n-1){
    int pm=subm[choice[i-1]];
    static u64 RL[16], RC[16]; static u64 Lm[65536], Cm[65536];
    for(int b=0;b<u;b++){
      u64 T=0; for(int a=0;a<u;a++) if((pm>>a&1) && a!=b) T|=S[i-1][a];
      RL[b]=shift(T,U[b]); u64 acc=0;
      for(int f=0;f<u;f++){ if(f==b) continue; u64 TC=0; for(int a=0;a<u;a++) if((pm>>a&1) && a!=b) TC|=SC[i-1][f][a]; acc|=shift(TC,U[b]); }
      RC[b]=acc;
    }
    Lm[0]=0; Cm[0]=0;
    for(int mm=1; mm<(1<<u); mm++){ int lb=__builtin_ctz(mm); Lm[mm]=Lm[mm&(mm-1)]|RL[lb]; Cm[mm]=Cm[mm&(mm-1)]|RC[lb]; }
    for(int t2=0;t2<nsub;t2++){
      if(eqk && subsz[t2]!=eqk) continue;
      int mm=subm[t2]; int s3=s+subsz[t2]; choice[i]=t2; tot+=curw;
      int l=__builtin_popcountll(Lm[mm]), c=__builtin_popcountll(Cm[mm]);
      if(l<boundL(s3)){ cntL+=curw; hist[0][s3]+=curw; if(s3<minSL){minSL=s3; nshowL=0;} if(SHOW<0 || (s3==minSL && nshowL<SHOW)){ printex("L-viol",s3,l); nshowL++; } }
      if(c<boundC(s3)){ cntC+=curw; hist[1][s3]+=curw; if(s3<minSC){minSC=s3; nshowC=0;} if(SHOW<0 || (s3==minSC && nshowC<SHOW)){ printex("C-viol",s3,c); nshowC++; } }
    }
    return;
  }
  for(int t=0;t<nsub;t++){
    if(eqk && subsz[t]!=eqk) continue;
    if(i==0 && !isrep[t]) continue;
    int m=subm[t]; choice[i]=t; int s2=s+subsz[t];
    if(i==0){
      curw=wt[t];
      for(int a=0;a<u;a++){ S[0][a]= (m>>a&1)? ((u64)1<<U[a]) : 0; for(int f=0;f<u;f++) SC[0][f][a]= (f==a && (m>>a&1))? ((u64)1<<U[a]) : 0; }
      rec(1,s2); continue;
    }
    int pm=subm[choice[i-1]];
    for(int b=0;b<u;b++){
      if(!(m>>b&1)){ S[i][b]=0; for(int f=0;f<u;f++) SC[i][f][b]=0; continue; }
      u64 T=0; for(int a=0;a<u;a++) if((pm>>a&1) && a!=b) T|=S[i-1][a];
      S[i][b]=shift(T,U[b]);
      for(int f=0;f<u;f++){ u64 TC=0; for(int a=0;a<u;a++) if((pm>>a&1) && a!=b) TC|=SC[i-1][f][a]; SC[i][f][b]=shift(TC,U[b]); }
    }
    rec(i+1,s2);
  }
}
int main(int argc,char**argv){
  m1=atoi(argv[1]); m2=atoi(argv[2]); n=atoi(argv[3]); int umode=atoi(argv[4]); eqk= argc>5? atoi(argv[5]):0; int pcap= argc>6? atoi(argv[6]):0;
  if(argc>7) SHOW=atoi(argv[7]);
  int fix01 = argc>8? atoi(argv[8]) : 0;
  N=m1*m2;
  for(int a=0;a<N;a++) for(int b=0;b<N;b++){ int a1=a%m1,a2=a/m1,b1=b%m1,b2=b/m1; ad[a][b]=((a1+b1)%m1)+m1*((a2+b2)%m2); }
  pG=1<<30; for(int a=1;a<N;a++){ int o=1,x=a; while(x){ x=ad[x][a]; o++; } if(o<pG) pG=o; }
  if(pcap<0) pG=0;
  for(int b=0;b<N;b++) for(int j=0;j<8;j++) for(int v=0;v<256;v++){ u64 r=0; for(int t=0;t<8;t++) if(v>>t&1){ int x=8*j+t; if(x<N) r|=(u64)1<<ad[x][b]; } shtab[b][j][v]=r; }
  if(umode==0){ u=N; for(int a=0;a<N;a++) U[a]=a; } else { u=umode; for(int a=0;a<u;a++) U[a]=a; }
  nsub=0; for(int m=0;m<(1<<u);m++){ int c=__builtin_popcount(m); if(c>=2){ subm[nsub]=m; subsz[nsub]=c; nsub++; } }
  for(int t=0;t<nsub;t++){ isrep[t]=1; wt[t]=1; }
  if(umode==0 && n>1){
    static int idx[65536]; for(int m=0;m<(1<<u);m++) idx[m]=-1; for(int t=0;t<nsub;t++) idx[subm[t]]=t;
    static int seen[65536]; memset(seen,0,sizeof seen);
    int nunits=0, units[64];
    if(m2==1){ for(int c=1;c<m1;c++){ int g=c,h=m1; while(h){int tt=g%h; g=h; h=tt;} if(g==1) units[nunits++]=c; } } else { units[0]=1; nunits=1; }
    for(int t=0;t<nsub;t++){ if(seen[t]) continue; long long osz=0;
      for(int ci=0;ci<nunits;ci++) for(int v=0;v<N;v++){ int mm=0; for(int a=0;a<u;a++) if(subm[t]>>a&1){ int x=U[a]; int y=(m2==1)? (units[ci]*x)%m1 : x; y=ad[y][v]; mm|=1<<y; }
        int t2=idx[mm]; if(!seen[t2]){ seen[t2]=1; osz++; if(t2!=t) isrep[t2]=0; } }
      wt[t]=osz; }
  }
  if(fix01){ for(int t=0;t<nsub;t++){ isrep[t]=(subm[t]==3); wt[t]=1; } }
  printf("G=Z/%d x Z/%d p(G)=%d n=%d universe=%d subsets=%d eqsize=%d fixA1={0,1}:%d\n", m1,m2,pG,n,u,nsub,eqk,fix01);
  rec(0,0);
  printf("tuples=%lld  L-violations=%lld (min S %d)  C-violations=%lld (min S %d)\n", tot,cntL,cntL?minSL:-1,cntC,cntC?minSC:-1);
  for(int w=0;w<2;w++){ int any=0; for(int s=0;s<64;s++) if(hist[w][s]){ if(!any) printf("%s violations by S:", w?"C":"L"); any=1; printf(" %d:%lld",s,hist[w][s]); } if(any) printf("\n"); }
  return 0;
}

