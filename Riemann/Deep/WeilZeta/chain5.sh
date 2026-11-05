#!/bin/sh
# chain5.sh: start queue5.sh (and the memory watchdog) once queue4.sh has finished (one compute process at a time)
cd /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta
while pgrep -f 'queue4.sh' > /dev/null; do sleep 10; done
./queue5.sh &
sleep 2; ./memwatch.sh &
wait
