
BT, BS = chr(96), chr(92)
p = "Riemann/Deep/B3Cert/LEAN.md"; s = open(p).read()
i, j = s.index("## Build"), s.index("## Review")
b = s[i:j]
b = b.replace(BS + BT, BT).replace(" " + BS + BS + "\n", " " + BS + "\n").replace("'^" + BS + BS + "s*axiom" + BS + BS + "s'", "'^" + BS + "s*axiom" + BS + "s'")
s = s[:i] + b + s[j:]
open(p, "w").write(s)
print(b.count(BS + BT), b.count(BS + BS))

