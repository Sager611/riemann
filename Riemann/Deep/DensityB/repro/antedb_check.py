
# densityb antedb_check.py: run ANTEDB's own Corollary 11.7 pipeline (zd.lv_zlv_to_zd) with and without
# the new Montgomery range (our hypothesis), on a sigma interval. Own driver script; imports ANTEDB at runtime.
import sys, os, time
PY = os.environ.get("EXPDB_PY", "expdb/blueprint/src/python")  # ANTEDB python dir (github.com/teorth/expdb, Apache-2.0), not included
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(PY); sys.path.insert(0, PY)
t0 = time.time()
from literature import *
import exponent_pair as ep
import bound_beta as bbeta
import large_values as lv
import zero_density_estimate as zd
from reference import Reference
from fractions import Fraction as frac

B = frac(87338196661320752, 10**17)
h = frac(1,2)
# tau_B(s) = (2s - 2 + B)/(B - 1/2);  row: 2 - 2s + t - tau_B(s) = c + a s + 1 t
a = -2 - 2/(B-h)
c = 2 - (B-2)/(B-h)
new_lv = lv.literature_bound_LV_max([[2,-2,0],[c,a,1]], Reference.make("densityb half-plane Montgomery range", 2026))

def base_set(with_new):
    hs = Hypothesis_Set()
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Exponent pair"))
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Exponent pair transform"))
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Upper bound on beta"))
    hs.add_hypotheses(ep.compute_exp_pairs(hs, search_depth=5, prune=True))
    hs.add_hypotheses(bbeta.exponent_pairs_to_beta_bounds(hs))
    hs.add_hypotheses(bbeta.compute_best_beta_bounds(hs))
    hs.add_hypothesis(lv.large_value_estimate_L2)
    for k in range(2, 5):
        hs.add_hypothesis(lv.raise_to_power_hypothesis(k))
    names = ["Huxley large value estimate", "Heath-Brown large value estimate"]
    for hyp in literature.list_hypotheses(hypothesis_type="Large value estimate"):
        nm = hyp.name
        if ("Jutila" in nm and any(nm.endswith(f"k = {k}") for k in (1,2,3,4))) or "Huxley" in nm or ("Heath-Brown" in nm and "large value estimate" in nm) or "Bourgain" in nm or nm.startswith("Guth--Maynard (2024) large value estimate") :
            hs.add_hypothesis(hyp)
    hs.add_hypotheses(literature.list_hypotheses(hypothesis_type="Zeta large value estimate"))
    if with_new: hs.add_hypothesis(new_lv)
    return hs

lo, hi = frac(sys.argv[1]), frac(sys.argv[2]); tau0 = frac(sys.argv[3])
for with_new in (False, True):
    hs = base_set(with_new)
    print("with_new", with_new, "#LV", len(hs.list_hypotheses(hypothesis_type="Large value estimate")), flush=True)
    from polytope import Polytope
    from region import Region
    lvr = zd.compute_large_value_region(hs, Region.from_polytope(Polytope.rect((lo, hi), (tau0, 2*tau0), (0, Constants.LV_DEFAULT_UPPER_BOUND))), zeta=False)
    sup1 = zd.compute_sup_rho_on_tau([r.child for r in lvr.data.region.child], Interval(lo, hi))
    for item in sup1:
        f, iv = item[0], item[1]
        print("   sup LV/tau on", iv, ":", f, flush=True)
        for sg in (lo, (lo+hi)/2, hi):
            if iv.contains(sg):
                v = f.at(sg) if hasattr(f,'at') else f(sg)
                print("      sigma", float(sg), " A <=", float(v/(1-sg)), flush=True)
    print("  time", round(time.time()-t0,1), flush=True)
