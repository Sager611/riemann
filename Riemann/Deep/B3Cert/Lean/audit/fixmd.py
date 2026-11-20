
p = "Riemann/Deep/B3Cert/LEAN.md"; s = open(p).read()
old = "| ✓; 3.35·10⁻⁵ is the brief's bound; `e1p_enclosures` adds h(371/750; κ₂) ∈ (0.000033545, 0.000033547] ∋ 3.3546·10⁻⁵ (Arb 3.354573606·10⁻⁵) |"
new = "| ✓; 3.35·10⁻⁵ is the brief's bound; `e1p_enclosures` adds h(371/750; κ₂) ∈ (0.000033545, 0.000033547] ∋ 3.354573606·10⁻⁵ (Arb, FP §8 l. 175). REPORT's '= 3.3546·10⁻⁵' and 'κ = 0.08503876039176576024' are roundings (≈) of the Arb balls h = 3.354573606042305490255·10⁻⁵ ± 1.5·10⁻³⁰, κ = 0.0850387603917657602445070238 ± 5·10⁻³² (reviewer pass 1); Lean keeps κ₂ and h symbolic and proves only inequalities |"
assert old in s
s = s.replace(old, new)
s = s.replace("| REPORT l. 39–40, FP §8 l. 177–178", "| REPORT l. 39, FP §8 l. 177–178").replace("| REPORT l. 36–38 (step 5); FP §8 l. 172", "| REPORT l. 36–37 (step 5); FP §8 l. 172")
open(p, "w").write(s)
print("ok")

