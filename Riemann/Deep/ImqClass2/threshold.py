# threshold.py: class-number thresholds T(N) from the certified Littlewood78 cells (cert.py, unchanged).
# Loads the definitions of Riemann/Deep/Littlewood78/cert.py (imports, constants, functions) without running
# its checks, then runs its own class_threshold(H, guess): the smallest e0 (multiple of 0.01) such that every
# cell [l, l+0.01] of l = log d from e0 to l_end certifies h(-d) > H via exp(F(l+0.01)) (arb balls), and
# claim (A) of the certificate (|L(1,chi)| >= 0.11/max(3, loglog q), all cells PASS in cert_output.txt)
# gives h(-d) >= 0.11 sqrt(d)/(pi loglog d) > H beyond l_end.
# Run: uv run --no-project --with python-flint python3 threshold.py 2000 3000 ... > thresholds.txt
import ast, math, sys, os
here = os.path.dirname(os.path.abspath(__file__))
src_path = os.path.join(here, '..', 'Littlewood78', 'cert.py')
src = open(src_path).read()
tree = ast.parse(src)
keep = []
seen_say = False
for node in tree.body:
    if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef)):
        keep.append(node)
    elif isinstance(node, ast.Assign) and not seen_say:
        keep.append(node)
    elif isinstance(node, ast.Expr) and isinstance(node.value, ast.Call) and getattr(node.value.func, 'id', '') == 'say':
        seen_say = True
mod = ast.Module(body=keep, type_ignores=[])
g = {'__name__': 'cert_defs'}
exec(compile(mod, src_path, 'exec'), g)
import hashlib
print(f"# cert.py sha256 {hashlib.sha256(src.encode()).hexdigest()}  B = {g['B_STR']}  cA = {g['CA']}")
print("# N  e0=log d threshold  d_threshold(rounded up)  l_end(claim A beyond)")
for a in sys.argv[1:]:
    H = int(a)
    guess = 2*math.log(H*math.pi/0.035)
    e0, e1 = g['class_threshold'](H, guess)
    v0 = math.exp(e0); k0 = math.floor(math.log10(v0)); up = math.ceil(v0/10**k0*1000)/1000
    print(f"{H} {e0:.2f} {up}e{k0} {e1:.2f}", flush=True)
