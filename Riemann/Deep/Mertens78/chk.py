import json, sys, time
from cert import Cert
d = json.load(open(sys.argv[1]))
c = Cert(d); t0 = time.time()
print("Lam", c.Lam, "a", c.a)
print(c.check_C(), "time", round(time.time()-t0, 1))
