import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreShape'))
from exact4 import check4
import numpy as np
# my J=4 LP form (a/a1, phi/a1) at kappa=3.358 (logs/pg_k3.358_J1-5.log)
a = [-0.95767, 1, 0.011945, -0.023295, 0.0025546]; p12, p32 = 0.03403, 0.02335
print('mine (5 digits):', check4(a, p12, p32, ng=400, KM=20000))
# explore53 certified J4 (kappa 3.34)
a2 = [-5899179108889603/6134092019850000, 1, 14789027537/948575124080, -4037122085/176350902122, 1982999731/771179571628]
print('explore53 cert4:', check4(a2, 3691572375143311/104390384155500000, 4037122085/176350902122, ng=400, KM=20000))
