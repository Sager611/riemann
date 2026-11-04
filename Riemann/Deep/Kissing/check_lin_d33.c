// d33full.c : brute-force exact check of Lin's 370240-point configuration in R^33 (ERS: Cheng-Sloane [32,17,8] dense layer,
// 1852 weight-8 supports x 128 even sign patterns, 2112 roots). All points have squared norm 32; need <x,y> <= 16 for all distinct pairs.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <arm_neon.h>
#define NPT 370240
static int8_t P[NPT][48] __attribute__((aligned(16)));
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); char line[512]; uint64_t G[32]; int ng=0;
  while(fgets(line,sizeof line,f)){ if(line[0]=='#') continue; uint64_t m=0; int k=0; for(char*p=line;*p;p++) if(*p=='0'||*p=='1'){ if(*p=='1') m|=1ull<<k; k++; } if(k==32) G[ng++]=m; } fclose(f);
  if(ng!=17){ fprintf(stderr,"gen rows %d\n",ng); return 1; }
  int n=0; for(int s=0;s<(1<<17);s++){ uint64_t c=0; for(int b=0;b<17;b++) if(s>>b&1) c^=G[b]; memset(P[n],0,48); for(int i=0;i<32;i++) P[n][i]=(c>>i&1)?-1:1; n++; }
  f=fopen(argv[2],"r"); int nsup=0; while(fgets(line,sizeof line,f)){ int pos[40],np=0,k=0; for(char*p=line;*p;p++) if(*p=='0'||*p=='1'){ if(*p=='1') pos[np++]=k; k++; } if(k!=33) continue; if(np!=8){ fprintf(stderr,"bad support weight\n"); return 1; } nsup++;
    for(int sg=0;sg<256;sg++){ if(__builtin_popcount(sg)&1) continue; memset(P[n],0,48); for(int i=0;i<8;i++) P[n][pos[i]]=(sg>>i&1)?-2:2; n++; } } fclose(f);
  for(int a=0;a<33;a++) for(int b=a+1;b<33;b++) for(int s=0;s<4;s++){ memset(P[n],0,48); P[n][a]=(s&1)?-4:4; P[n][b]=(s&2)?-4:4; n++; }
  printf("points: %d (supports %d)\n",n,nsup); if(n!=NPT){ fprintf(stderr,"count mismatch\n"); return 1; }
  for(int i=0;i<n;i++){ int s=0; for(int t=0;t<48;t++) s+=P[i][t]*P[i][t]; if(s!=32){ fprintf(stderr,"norm %d at %d\n",s,i); return 1; } }
  double t0=(double)clock()/CLOCKS_PER_SEC; long viol=0; int mx=-1000; int32x4_t c16=vdupq_n_s32(16), z=vdupq_n_s32(0); int32x4_t vmax=vdupq_n_s32(-1000);
  for(int i=0;i<n;i++){ int8x16_t a0=vld1q_s8(P[i]),a1=vld1q_s8(P[i]+16),a2=vld1q_s8(P[i]+32);
    int j=i+1;
    for(;j+4<=n;j+=4){
      int32x4_t s0=vdotq_s32(vdotq_s32(vdotq_s32(z,a0,vld1q_s8(P[j])),a1,vld1q_s8(P[j]+16)),a2,vld1q_s8(P[j]+32));
      int32x4_t s1=vdotq_s32(vdotq_s32(vdotq_s32(z,a0,vld1q_s8(P[j+1])),a1,vld1q_s8(P[j+1]+16)),a2,vld1q_s8(P[j+1]+32));
      int32x4_t s2=vdotq_s32(vdotq_s32(vdotq_s32(z,a0,vld1q_s8(P[j+2])),a1,vld1q_s8(P[j+2]+16)),a2,vld1q_s8(P[j+2]+32));
      int32x4_t s3=vdotq_s32(vdotq_s32(vdotq_s32(z,a0,vld1q_s8(P[j+3])),a1,vld1q_s8(P[j+3]+16)),a2,vld1q_s8(P[j+3]+32));
      int32x4_t s=vpaddq_s32(vpaddq_s32(s0,s1),vpaddq_s32(s2,s3));
      vmax=vmaxq_s32(vmax,s);
      uint32x4_t m=vcgtq_s32(s,c16); if(vmaxvq_u32(m)){ int32_t sv[4]; vst1q_s32(sv,s); for(int t=0;t<4;t++) if(sv[t]>16) viol++; }
    }
    for(;j<n;j++){ int s=0; for(int t=0;t<48;t++) s+=P[i][t]*P[j][t]; if(s>mx) mx=s; if(s>16) viol++; }
    if((i&65535)==0){ fprintf(stderr,"  i=%d (%.0fs)\n",i,(double)clock()/CLOCKS_PER_SEC-t0); }
  }
  int vm=vmaxvq_s32(vmax); if(vm>mx) mx=vm;
  long long pairs=(long long)n*(n-1)/2;
  printf("pairs checked: %lld; max inner product %d (limit 16 = |x|^2/2; 32 would be a repeat); violations %ld\n",pairs,mx,viol);
  printf(viol? "INVALID\n" : "VALID kissing configuration of %d points in R^33\n", n); return viol?1:0; }
