// tmsim.c - independent simulator for the RH 120-state table with register-trace instrumentation.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXS 512
static int W[MAXS][2], D[MAXS][2], T[MAXS][2];
static int nst = 0;
static int leafReg[MAXS], leafOp[MAXS];
static int stRet0=-1, stRet1=-1, stRoot2=-1, stRootN=-1;
static unsigned char *tape; static long TS=1L<<25, ORG=1L<<16;

static void die(const char*m){fprintf(stderr,"ERR %s\n",m);exit(2);}

static void load_table(const char*fn){
  FILE*f=fopen(fn,"r"); if(!f) die("table");
  char line[256]; for(int i=0;i<MAXS;i++){T[i][0]=T[i][1]=-2;}
  while(fgets(line,sizeof line,f)){
    if(line[0]=='#'||line[0]=='\n') continue;
    int q,s,w; char d, n[64];
    if(sscanf(line,"%d %d %d %c %63s",&q,&s,&w,&d,n)!=5) die("parse");
    W[q][s]=w; D[q][s]=(d=='R')?1:-1; T[q][s]=(n[0]=='H')?-1:atoi(n);
    if(q+1>nst) nst=q+1;
  }
  fclose(f);
  for(int q=0;q<nst;q++) for(int s=0;s<2;s++) if(T[q][s]==-2) die("missing transition");
}
static void load_roles(const char*fn){
  for(int i=0;i<MAXS;i++){leafReg[i]=-1;leafOp[i]=-1;}
  FILE*f=fopen(fn,"r"); if(!f) die("roles");
  char role[64]; int id;
  while(fscanf(f,"%63s %d",role,&id)==2){
    if(!strcmp(role,"ret0")) stRet0=id; else if(!strcmp(role,"ret1")) stRet1=id;
    else if(!strcmp(role,"root2")) stRoot2=id; else if(!strcmp(role,"rootN")) stRootN=id;
    else if(!strncmp(role,"leaf",4)){ // leaf<reg><op>  e.g. leaf20 = reg2 inc, leaf21 = reg2 dec
      leafReg[id]=role[4]-'0'; leafOp[id]=role[5]-'0'; }
  }
  fclose(f);
}
// decode 5 registers starting at coordinate 5; returns end coordinate (cell after 5th delimiter)
static long decode(long long r[5]){
  long c=5;
  for(int k=0;k<5;k++){ long long n=0; while(tape[c+ORG]==1){n++;c++; if(c+ORG>=TS-2) die("decode overflow");}
    if(n==0) return -1; r[k]=n-1; c++; }
  return c;
}
static long tail_blocks(long c,int*ok){ // count "10" blocks after c until "00"
  long nb=0; *ok=1;
  while(1){ if(tape[c+ORG]==1){ if(tape[c+ORG+1]!=0){*ok=0;return nb;} nb++; c+=2;} else { if(tape[c+ORG+1]!=0) *ok=0; return nb; } }
}
static void rewrite(long long nr[5]){
  long long r[5]; long end=decode(r); if(end<0) die("rewrite decode");
  // find end of used tape (after tail's 00)
  long e=end; while(!(tape[e+ORG]==0&&tape[e+ORG+1]==0&&tape[e+ORG+2]==0&&tape[e+ORG+3]==0)) e++;
  long tl=e-end; unsigned char*tmp=malloc(tl+1); memcpy(tmp,tape+end+ORG,tl);
  long c=5; long long need=0; for(int k=0;k<5;k++) need+=nr[k]+2;
  if(5+need+tl+16+ORG>=TS) die("rewrite too big");
  for(long z=5; z<e+16; z++) tape[z+ORG]=0;
  for(int k=0;k<5;k++){ for(long long i=0;i<nr[k]+1;i++) tape[ORG+c++]=1; tape[ORG+c++]=0; }
  memcpy(tape+ORG+c,tmp,tl); free(tmp);
}
static const char*RN="mjxfv";
int main(int argc,char**argv){
  if(argc<3) die("usage: tmsim table roles [opts]");
  load_table(argv[1]); load_roles(argv[2]);
  long long maxsteps=4000000000000LL; int maxstage=4; long long saveK=-1; const char*savefn=0; const char*loadfn=0;
  long long nr[5]={0}; int checkevery=0; long long untilops=-1;
  for(int i=3;i<argc;i++){
    if(!strcmp(argv[i],"--maxsteps")) maxsteps=atoll(argv[++i]);
    else if(!strcmp(argv[i],"--maxstage")) maxstage=atoi(argv[++i]);
    else if(!strcmp(argv[i],"--save")){ saveK=atoll(argv[++i]); savefn=argv[++i]; }
    else if(!strcmp(argv[i],"--load")){ loadfn=argv[++i]; for(int k=0;k<5;k++) nr[k]=atoll(argv[++i]); }
    else if(!strcmp(argv[i],"--check")) checkevery=atoi(argv[++i]);
    else die("opt");
  }
  tape=calloc(TS,1); if(!tape) die("mem");
  int q=0; long h=0; long long steps=0; int started=0; long long gops=0, sops=0; uint64_t hs=0;
  int pending=-1; long long tr[5]={0}; int stage=2; long minh=0,maxh=0;
  if(loadfn){
    FILE*f=fopen(loadfn,"rb"); if(!f) die("loadfile"); long lo,hi; 
    if(fread(&q,sizeof q,1,f)!=1||fread(&h,sizeof h,1,f)!=1||fread(&steps,sizeof steps,1,f)!=1||fread(&lo,sizeof lo,1,f)!=1||fread(&hi,sizeof hi,1,f)!=1) die("rd");
    if(fread(tape+lo+ORG,1,hi-lo+1,f)!=(size_t)(hi-lo+1)) die("rd2"); fclose(f);
    rewrite(nr); started=1; for(int k=0;k<5;k++) tr[k]=nr[k]; stage=-1; steps=0;
    printf("LOADED q=%d h=%ld regs=%lld %lld %lld %lld %lld\n",q,h,nr[0],nr[1],nr[2],nr[3],nr[4]);
  }
  while(steps<maxsteps){
    unsigned char s=tape[h+ORG];
    if(!started && q==stRootN){
      started=1; long long r[5]; long e=decode(r); int ok; long nb=tail_blocks(e,&ok);
      long ones=0, lo=1L<<40, hi=-(1L<<40); int supp_ok=1;
      for(long c=minh-2;c<=maxh+2;c++) if(tape[c+ORG]){ones++; if(c<lo)lo=c; if(c>hi)hi=c;
         int good = (c==-3821||c==-3820) || (c>=3 && c<=3821 && ((c-3)%2==0)); if(!good) supp_ok=0; }
      long expected=2+1+1909; 
      printf("BOOT steps=%lld state=%d head=%ld ones=%ld lo=%ld hi=%ld support_exact=%d (expect ones=%ld) regs=%lld,%lld,%lld,%lld,%lld tailblocks=%ld tailok=%d\n",
        steps,q,h,ones,lo,hi,supp_ok&&ones==expected,expected,r[0],r[1],r[2],r[3],r[4],nb,ok);
      fflush(stdout);
    }
    if(started){
      if(s==1 && leafReg[q]>=0){
        if(pending>=0) die("op start while pending");
        if(saveK>=0 && gops==saveK){
          FILE*f=fopen(savefn,"wb"); long lo=minh-4, hi=maxh+4; 
          fwrite(&q,sizeof q,1,f); fwrite(&h,sizeof h,1,f); fwrite(&steps,sizeof steps,1,f); fwrite(&lo,sizeof lo,1,f); fwrite(&hi,sizeof hi,1,f);
          fwrite(tape+lo+ORG,1,hi-lo+1,f); fclose(f);
          long long r[5]; decode(r);
          printf("SAVED op=%lld state=%d head=%ld steps=%lld regs=%lld,%lld,%lld,%lld,%lld\n",gops,q,h,steps,r[0],r[1],r[2],r[3],r[4]);
          return 0;
        }
        pending=q;
        if(checkevery>0 && (gops%checkevery)==0){
          long long r[5]; if(decode(r)<0) die("decode fail");
          for(int k=0;k<5;k++) if(r[k]!=tr[k]){ printf("MISMATCH op=%lld reg %c tape=%lld tracked=%lld\n",gops,RN[k],r[k],tr[k]); return 1; }
        }
      }
      if(pending>=0 && (q!=stRet0 && q!=stRet1)){
        int t=T[q][s];
        if(t==stRet0 || t==stRet1){
          int out = (t==stRet1); int a=leafReg[pending], op=leafOp[pending];
          if(op==0 && out) die("inc returned zero-outcome");
          if(op==1){ if(out){ if(tr[a]!=0){printf("MISMATCH zero-outcome but tracked %c=%lld op=%lld\n",RN[a],tr[a],gops);return 1;} }
                     else { if(tr[a]==0){printf("MISMATCH dec success but tracked %c=0 op=%lld\n",RN[a],gops);return 1;} tr[a]--; } }
          else tr[a]++;
          hs = hs*1000003ULL + (uint64_t)(a*4+op*2+out);
          gops++; sops++; pending=-1;
        }
      }
      if(q==stRoot2 && s==1){
        long long r[5]; long e=decode(r); int ok; long nb=tail_blocks(e,&ok);
        int match=1; for(int k=0;k<5;k++) if(r[k]!=tr[k]) match=0;
        printf("RESTART stage_n=%d steps=%lld stage_ops=%lld hash=%016llx regs=%lld,%lld,%lld,%lld,%lld tracked_match=%d tailblocks=%ld tailok=%d minh=%ld maxh=%ld\n",
          stage,steps,sops,(unsigned long long)hs,r[0],r[1],r[2],r[3],r[4],match,nb,ok,minh,maxh);
        fflush(stdout);
        if(loadfn) return 0;
        sops=0; hs=0; stage++;
        if(stage>maxstage) return 0;
      }
    }
    int w=W[q][s], t=T[q][s];
    tape[h+ORG]=(unsigned char)w; h+=D[q][s]; steps++;
    if(h<minh){minh=h; if(h+ORG<4) die("left overflow");}
    if(h>maxh){maxh=h; if(h+ORG>=TS-8) die("right overflow");}
    if(t<0){
      long long r[5]; long e=decode(r);
      printf("HALT steps=%lld from_state=%d read=%d stage_n=%d stage_ops=%lld hash=%016llx pending=%d regs(after)=%lld,%lld,%lld,%lld,%lld decode_ok=%d\n",
        steps,q,s,stage,sops,(unsigned long long)hs,pending,r[0],r[1],r[2],r[3],r[4],e>0);
      return 0;
    }
    q=t;
  }
  printf("STEPLIMIT steps=%lld stage=%d stage_ops=%lld\n",steps,stage,sops);
  return 0;
}
