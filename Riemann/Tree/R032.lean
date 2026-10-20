/-!
# Round 32 — elementary and pretentious methods: no power saving for μ

| method | best proven | barrier |
|--------|-------------|---------|
| Selberg symmetry / Erdős–Selberg PNT (Bombieri, Diamond–Steinig) | ψ(x) = x + O(x e^{-c√log x})-type | parity, averaged correlations |
| Halász / pretentious | logarithmic savings | o(x) only |
| Matomäki–Radziwiłł, Tao log-Chowla, Tao–Teräväinen | almost-all short intervals | no uniform rate |
| Wirsing / Hall–Tenenbaum | nonnegative functions | blind to μ's signs |
| Vaughan / Heath-Brown identities | Korobov–Vinogradov saving | power saving = new strip |
| Bombieri–Vinogradov, Zhang, Maynard | level of distribution, bounded gaps | averages over moduli |
| Chowla / Sarnak | log-averaged, special systems | no power saving |
| Linnik dispersion / BDH | mean-square in progressions | exceptional moduli |

All are Beurling-sensitive (use integer structure) yet stop at o(x) or exp-type savings;
M(x) = O(x^{1-δ}) would give the quasi-RH strip (machine-checked: Fayed.R19.quasiRH_of_mertensBound).
-/
