#include <libproc.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
int main(int argc, char** argv){
  int pid = atoi(argv[1]);
  int ms = argc>2? atoi(argv[2]) : 50;
  int every = argc>3? atoi(argv[3]) : 10;
  unsigned long long maxfp=0, maxrss=0; long n=0;
  while (kill(pid,0)==0) {
    struct rusage_info_v2 ri;
    if (proc_pid_rusage(pid, RUSAGE_INFO_V2, (rusage_info_t*)&ri)==0){
       unsigned long long fp = ri.ri_phys_footprint, rss = ri.ri_resident_size;
       if (fp>maxfp) maxfp=fp; if (rss>maxrss) maxrss=rss;
       if (n % every == 0) { printf("t=%.2fs rss=%lluMB fp=%lluMB\n", n*ms/1000.0, rss>>20, fp>>20); fflush(stdout);}
    }
    n++;
    usleep(ms*1000);
  }
  printf("max rss=%llu MB max fp=%llu MB\n", maxrss>>20, maxfp>>20);
  return 0;
}

