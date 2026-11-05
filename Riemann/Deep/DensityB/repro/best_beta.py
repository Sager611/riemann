
import sys, os, time, json
PY = os.environ.get("EXPDB_PY", "expdb/blueprint/src/python")  # ANTEDB python dir (github.com/teorth/expdb, Apache-2.0), not included
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(PY); sys.path.insert(0, PY)
t0 = time.time()
from literature import *
import exponent_pair as ep
import bound_beta as bbeta
hs = Hypothesis_Set()
hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Exponent pair"))
hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Exponent pair transform"))
hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Upper bound on beta"))
print("n hyps", len(hs), file=sys.stderr)
hs.add_hypotheses(ep.compute_exp_pairs(hs, search_depth=int(sys.argv[1]) if len(sys.argv)>1 else 5, prune=True))
hs.add_hypotheses(bbeta.exponent_pairs_to_beta_bounds(hs))
best = bbeta.compute_best_beta_bounds(hs)
out = []
for b in best:
    bd = b.data.bound
    out.append([str(bd.domain.x0), str(bd.domain.x1), str(bd.c), str(bd.m)])
out.sort(key=lambda r: float(eval(r[0].replace('/','*1.0/'))) if r[0] not in ('inf','-inf') else 0)
json.dump(out, open(os.path.join(HERE, "best_beta.json"),"w"), indent=0)
print("pieces", len(out), "time", round(time.time()-t0,1))
for r in out:
    print(r)
