import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreShape'))
from exact4 import check4
a = [float(x) for x in sys.argv[1:6]]; p12, p32 = float(sys.argv[6]), float(sys.argv[7])
print('check4 (exact Lemma 2, 400 g, KM=20000): pmin, (worst S, g, parity) =', check4(a, p12, p32, ng=400, KM=20000))
