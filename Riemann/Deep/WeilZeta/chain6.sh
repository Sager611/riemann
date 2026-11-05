#!/bin/sh
# chain6.sh: start queue6.sh once queue5.sh has finished (one compute process at a time)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
while pgrep -f 'queue5.sh' > /dev/null; do sleep 10; done
./queue6.sh
