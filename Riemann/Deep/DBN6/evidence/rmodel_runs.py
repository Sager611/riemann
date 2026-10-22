# Runs of the Raymond-chain model: baseline (no derivative) and derivative term with assumed G = fac*psi2/4 + add.
import importlib.util
spec = importlib.util.spec_from_file_location('rm', 'rmodel_lib.py'); rm = importlib.util.module_from_spec(spec); spec.loader.exec_module(rm)
print('model, Raymond kernels, no derivative:', rm.run(False))
for fac, add in [(1.0, 0.0), (1.2, 0.1), (1.5, 0.2), (2.0, 0.3)]:
    print('derivative, G = %.1f*psi2/4 + %.1f:' % (fac, add), rm.run(True, fac, add))
