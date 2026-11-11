# PYAudit: does the 54.6 certificate depend on the reading of the recursion in Lemma 1.4?  Re-evaluates SubWeyl/certify.py's
# interval bound at the binding intervals with A_k, B_k computed under (proof) the reading of Yang's proof, (alt) delta_3 of the
# third-derivative test reused in the A_3 -> A_4 step.  uv run --no-project --with python-flint --with numpy python3 reading_check.py
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SubWeyl'))
sys.dont_write_bytecode = True   # keep SubWeyl/ free of __pycache__
import certify as C
from flint import arb
orig = C.AkBk
def AkBk_alt(k, eta, h):
    A_, B_ = C.A3B3(eta, h)
    d3 = (C.fr(1, 2) + (1 + C.fr(3, 8)*C.PI.sqrt()*eta**C.fr(3, 2)).sqrt()/2).sqrt()
    for j in range(3, k):
        J = 2**(j - 1)
        d = (1 + 2/arb(2337)**(1 - C.fr(2, J))*(9*C.PI/1024*eta)**C.fr(1, J)).sqrt()
        if j == 3: d = d3
        A_ = d*(h**(-C.fr(1, J)) + arb(2)**C.fr(19, 12)*(J - 1)/arb((2*J - 1)*(4*J - 3)).sqrt()*A_.sqrt())
        B_ = d*arb(2)**C.fr(3, 2)*(J - 1)/arb((2*J - 3)*(4*J - 5)).sqrt()*B_.sqrt()
    return A_, B_
for name, fn in (('proof', orig), ('alt', AkBk_alt)):
    C.AkBk = fn
    for (L0, L1) in ((648.0, 650.0), (638.0, 640.0), (540.0, 542.0)):
        B = C.interval_bound(L0, L1)
        print('%-5s L in [%g, %g]: bound %.5f' % (name, L0, L1, C.up(B))); sys.stdout.flush()
    print('%-5s tail t >= e^2000: %.5f' % (name, C.up(C.tail_bound(*C.TAIL))))
