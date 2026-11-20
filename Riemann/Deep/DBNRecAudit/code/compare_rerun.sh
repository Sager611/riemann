#!/bin/bash
# dbnrecaudit compare_rerun.sh: /tmp/dbnra rerun of DBNRec2 row 8 vs the logged files (prep.json line 1, bb.out without timings/memcap, unres.txt, cert JSON).
R=/Users/adriansager/playground/riemann/Riemann/Deep/DBNRec2/logs; V=/tmp/dbnra/logs
for j in 00 01 02 03 04 05; do t=r8_$j
  a=$(head -1 $V/$t.prep.json | cmp -s - <(head -1 $R/$t.prep.json) && echo same || echo DIFF)
  b=$(diff -q <(rg -v 'memcap|maximum resident|peak memory' $R/$t.bb.out | sed -E 's/[0-9.]+ s$//; s/tables [0-9.]+ s/tables/') <(sed -E 's/[0-9.]+ s$//; s/tables [0-9.]+ s/tables/' $V/$t.bb.out) > /dev/null && echo same || echo DIFF)
  c=$(cmp -s $V/$t.unres.txt $R/$t.unres.txt && echo same || echo DIFF)
  echo "$t prep.json:$a bb.out:$b unres.txt:$c $(rg -o 'memcap: peak [0-9]+ MB' $V/$t.bb.err) $(rg -o '[0-9]+  maximum resident' $V/$t.bb.err)"
done
python3 - <<'PY'
import json
a = json.load(open('/tmp/dbnra/logs/cert_r8.json')); b = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/DBNRec2/logs/cert_r8.json'))
print('cert identical:', a == b, '| A_new', a['A_new'] == b['A_new'], '| B_new', a['B_new'] == b['B_new'], '| A', a['A_new_float'], 'B', a['B_new_float'])
PY
