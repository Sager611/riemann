import glob, os, math
from decimal import Decimal, getcontext, ROUND_FLOOR
getcontext().prec = 50
base = "/Users/adriansager/playground/riemann/.scratch/bothaaudit/b/AlexGenusCore-weil-ground-state-simple-even-04f6404/data/piece4/glue"
Ls = sorted({os.path.basename(f).split('_')[1] for f in glob.glob(base + "/g_*.json")}, key=float)
lines = ["30 240 660 -"]
for L in Ls:
    c = (Decimal(L)*2).exp()
    if c <= 30: continue
    cs = str(c.quantize(Decimal('0.000001'), rounding=ROUND_FLOOR)); cf = float(cs)
    N = max(240, int(round(cf*cf/5/10))*10); bits = math.ceil(3.33*(150 + 5.46*(cf - 30))) + 160
    lines.append(f"{cs} {N} {bits} {L}")
open('logs/grid.txt', 'w').write(chr(10).join(lines) + chr(10))
print(len(lines), lines[1], lines[-1])
