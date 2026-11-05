// ls2.c - WalkSAT/tabu local search: variables col[v] in Z3; hyperedge items (v,s) show colour col[v]+s;
// each hyperedge must show all 3 colours.  input: "m NE" then NE lines of items.
// usage: ls2 <file> <seconds> <seed> <shiftmode 0|1> [nvars]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
static double wall_since(struct timespec*t0){ struct timespec t; clock_gettime(CLOCK_MONOTONIC,&t); return (t.tv_sec-t0->tv_sec)+1e-9*(t.tv_nsec-t0->tv_nsec); }
static int m, NE, nv; static int *eoff, *ev, *es; static int *voff, *vedg, *vsh;
static int *col; static unsigned char (*cnt)[3];
static int *viol, *vpos, nviol;
static uint64_t rs; static inline uint64_t rnd(void){ rs^=rs<<13; rs^=rs>>7; rs^=rs<<17; return rs; }
static void addv(int code){ vpos[code]=nviol; viol[nviol++]=code; }
static void remv(int code){ int p=vpos[code]; int last=viol[--nviol]; viol[p]=last; vpos[last]=p; vpos[code]=-1; }
static int zeros(const int*t){ return (t[0]==0)+(t[1]==0)+(t[2]==0); }
static int delta(int v, int b){
  int a=col[v]; if(a==b) return 0; int d=0;
  for(int j=voff[v]; j<voff[v+1]; ){
    int e=vedg[j]; int t[3]={cnt[e][0],cnt[e][1],cnt[e][2]}; int z0=zeros(t);
    while(j<voff[v+1] && vedg[j]==e){ int s=vsh[j]; t[(a+s)%3]--; t[(b+s)%3]++; j++; }
    d+=zeros(t)-z0;
  }
  return d;
}
static void apply(int v, int b){
  int a=col[v]; if(a==b) return;
  for(int j=voff[v]; j<voff[v+1]; j++){ int e=vedg[j], s=vsh[j]; int ca=(a+s)%3, cb=(b+s)%3;
    if(--cnt[e][ca]==0) addv(e*3+ca);
    if(cnt[e][cb]++==0) remv(e*3+cb); }
  col[v]=b;
}
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); double secs=atof(argv[2]); rs=0x9E3779B97F4A7C15ULL ^ ((uint64_t)atoll(argv[3])+1)*0xD1B54A32D192ED03ULL; int shiftmode=atoi(argv[4]);
  if(fscanf(f,"%d %d",&m,&NE)!=2) return 1;
  nv = shiftmode? (m+2)/3 : m; if(argc>5) nv=atoi(argv[5]);
  int cap=NE*16+16; ev=malloc(cap*sizeof(int)); es=malloc(cap*sizeof(int)); eoff=malloc((NE+1)*sizeof(int)); int pos=0;
  char line[8192]; if(!fgets(line,sizeof line,f)) return 1;
  for(int e=0;e<NE;e++){ if(!fgets(line,sizeof line,f)) return 1; eoff[e]=pos; char*p=line,*q;
    while(1){ long x=strtol(p,&q,10); if(p==q) break; if(pos>=cap){cap*=2; ev=realloc(ev,cap*sizeof(int)); es=realloc(es,cap*sizeof(int));}
      if(shiftmode){ ev[pos]=(int)(x/3); es[pos]=(int)(x%3);} else { ev[pos]=(int)x; es[pos]=0; } pos++; p=q; } }
  eoff[NE]=pos; fclose(f);
  int *deg=calloc(nv,sizeof(int)); for(int i=0;i<pos;i++) deg[ev[i]]++;
  voff=malloc((nv+1)*sizeof(int)); voff[0]=0; for(int v=0;v<nv;v++) voff[v+1]=voff[v]+deg[v];
  vedg=malloc(pos*sizeof(int)); vsh=malloc(pos*sizeof(int)); int*fl=calloc(nv,sizeof(int));
  for(int e=0;e<NE;e++) for(int i=eoff[e];i<eoff[e+1];i++){ int v=ev[i]; vedg[voff[v]+fl[v]]=e; vsh[voff[v]+fl[v]]=es[i]; fl[v]++; }
  col=malloc(nv*sizeof(int)); cnt=calloc(NE,3); viol=malloc(NE*3*sizeof(int)); vpos=malloc(NE*3*sizeof(int));
  long long *tabu=calloc(nv*3,sizeof(long long)); int bestc=1<<30;
  struct timespec ts0; clock_gettime(CLOCK_MONOTONIC,&ts0); long long it=0; int restarts=0;
  while(1){
    for(int v=0;v<nv;v++) col[v]=rnd()%3;
    memset(cnt,0,(size_t)NE*3); nviol=0; for(int i=0;i<NE*3;i++) vpos[i]=-1;
    for(int e=0;e<NE;e++){ for(int i=eoff[e];i<eoff[e+1];i++) cnt[e][(col[ev[i]]+es[i])%3]++; for(int c=0;c<3;c++) if(!cnt[e][c]) addv(e*3+c); }
    long long steps=0, maxsteps=200000LL+2000LL*nv; int lb=nviol; long long lastimp=0;
    while(1){
      it++; steps++;
      if(nviol<bestc){ bestc=nviol; }
      if(nviol==0) goto done;
      if(nviol<lb){ lb=nviol; lastimp=steps; }
      if(steps-lastimp>maxsteps) break;
      if((it&255)==0 && wall_since(&ts0)>secs) goto out;
      int code=viol[rnd()%nviol]; int e=code/3, c=code%3;
      int bestd=1<<30, nc=0; int cv[64], cb[64];
      for(int i=eoff[e];i<eoff[e+1];i++){ int v=ev[i]; int b=((c-es[i])%3+3)%3; if(col[v]==b) continue;
        int d=delta(v,b); int tb=tabu[v*3+b]>it; if(tb && nviol+d>=bestc) continue;
        if(d<bestd){bestd=d;nc=0;} if(d==bestd && nc<64){cv[nc]=v;cb[nc]=b;nc++;} }
      int v,b;
      if(nc==0 || (rnd()%1000)<80){ int i=eoff[e]+rnd()%(eoff[e+1]-eoff[e]); v=ev[i]; b=((c-es[i])%3+3)%3; if(col[v]==b) continue; }
      else { int k=rnd()%nc; v=cv[k]; b=cb[k]; }
      int a=col[v]; apply(v,b); tabu[v*3+a]=it+2+rnd()%(4+nv/8);
    }
    restarts++;
  }
done:
  printf("SOLVED it=%lld restarts=%d t=%.1f\nSOL", it, restarts, wall_since(&ts0)); for(int v=0;v<nv;v++) printf(" %d",col[v]); printf("\n"); fflush(stdout); return 0;
out:
  printf("FAIL best=%d it=%lld restarts=%d\n", bestc, it, restarts); fflush(stdout); return 2;
}

