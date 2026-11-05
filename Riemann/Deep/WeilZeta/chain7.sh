#!/bin/sh
# chain7.sh: start queue7.sh once queue6.sh has finished (one compute process at a time)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
while pgrep -f 'queue6.sh' > /dev/null; do sleep 10; done
./queue7.sh
