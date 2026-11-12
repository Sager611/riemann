import sys
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
sys.path.insert(0, HERE)
from ctest import arrays
for r in (1, 2):
    a, b, cc = arrays(3.30 * np.pi, r, 12)
    for name, h in (('a', a), ('b', b), ('cc', cc)):
        hs = h[:6, :6]
        print(r, name, 'max|h-h^H| (i,j<6) = %.2e, max|h-h^T| = %.2e, max|Im h| = %.2e, max|h| = %.2e' % (np.abs(hs - hs.conj().T).max(), np.abs(hs - hs.T).max(), np.abs(hs.imag).max(), np.abs(hs).max()))
    print(np.round(a[:4, :4], 4))
