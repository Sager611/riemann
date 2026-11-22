/* fpmem pid... : prints "<max physical footprint MB> <max RSS MB>" over the given pids.
   Physical footprint = what `top -l 1 -pid <pid> -stats mem` shows (man top: "Physical memory footprint"). */
#include <libproc.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
int main(int argc, char **argv) {
  unsigned long long fp = 0, rss = 0;
  for (int i = 1; i < argc; i++) {
    struct rusage_info_v2 ri;
    if (proc_pid_rusage(atoi(argv[i]), RUSAGE_INFO_V2, (rusage_info_t *)&ri) == 0) {
      if (ri.ri_phys_footprint > fp) fp = ri.ri_phys_footprint;
      if (ri.ri_resident_size > rss) rss = ri.ri_resident_size;
    }
  }
  printf("%llu %llu\n", fp >> 20, rss >> 20);
  return 0;
}
