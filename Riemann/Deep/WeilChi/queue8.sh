#!/bin/sh
# queue8.sh: verify the memory-lean code on m4 c=30 (margin must equal 2.23133452614e-36), then m3 c=47 odd (default knobs, K2=201, 320 bits), then even characters at c=40
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi
export WEILCHI_SCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/scan WEILCHI_CSCAN=/Users/adriansager/playground/riemann/.scratch/weilchi/cscan
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while pgrep -f 'cert_chi.py p5 40 even' > /dev/null; do sleep 3; done
R="uv run --no-project --with mpmath --with numpy --with python-flint python3"
perl -e 'alarm 600; exec @ARGV' $R cert_chi.py m4 30 even 0 448 201 >> logs/leancheck.log 2>&1
perl -e 'alarm 3500; exec @ARGV' /usr/bin/time -l $R cert_chi.py m3 47 odd 2^-275 320 201 >> logs/m3_47_cert.log 2>> logs/m3_47_odd2.time
perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_chi.py p8 40 even 2^-91 384 601 >> logs/p8_40_cert.log 2>> logs/p8_40.time
perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_chi.py p5 40 odd 2^-140 384 601 >> logs/p5_40_cert.log 2>> logs/p5_40.time
perl -e 'alarm 2400; exec @ARGV' /usr/bin/time -l $R cert_chi.py p8 40 odd 2^-91 384 601 >> logs/p8_40_cert.log 2>> logs/p8_40.time
echo QUEUE8 DONE >> logs/p8_40_cert.log
