import json, glob, math
rows = []
for f in sorted(glob.glob('logs/bt_*.json'), key=lambda f: json.load(open(f))['t']):
    r = json.load(open(f)); DL = r['DL']
    jE = min(range(len(DL)), key=lambda j: abs(DL[j] - r['Dtri']))
    j14 = DL.index(14)
    print(f"{r['name']}[{r['lo']}] t={r['t']:.3e} L={r['L']:.2f} Δ_tri={r['Dtri']:.2f} B_tri={r['Btri']:.2f} CCM¼={r['ccm']:.2f} "
          f"maxS={r['maxS']:.2f} minS={r['minS']:.2f} | ΔE={DL[jE]}: τ={r['tau'][jE]:.2f} TB={r['TB'][jE]:.2f} supOζ={r['supO_zeta'][jE]:.2f} "
          f"rand={r['supO_rand'][jE]:.2f} | Δ=14: TB={r['TB'][j14]:.1f} supOζ={r['supO_zeta'][j14]:.2f} rand={r['supO_rand'][j14]:.2f} "
          f"maxB+={r['maxBup'][j14]:.2f} | min marg={min(min(r['marg_up']), min(r['marg_lo'])):.4f} n={r['nsel']} {r['secs']}s")
