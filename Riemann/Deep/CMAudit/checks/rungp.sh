#!/bin/bash
# usage: bash Riemann/Deep/CMAudit/checks/rungp.sh SECONDS file.gp  (from repo root; inside memcap 290 MB; perl fork-and-kill timer)
t=$1; f=$2
OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 $((t+20)) perl -e '$t=shift; $p=fork; if(!$p){exec @ARGV; exit 1} $SIG{ALRM}=sub{kill 9,$p; exit 124}; alarm $t; waitpid($p,0); exit($?>>8)' $t gp -q -s 64000000 -D parisizemax=120000000 -D nbthreads=1 $f < /dev/null
