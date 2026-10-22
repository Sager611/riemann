# T2: rigorous bin masses vs direct float sums (cell 500, h = 21/25)
import numpy as np, math
from fractions import Fraction as Fr
import lib5 as L
cell = 500; h = 21/25
Psi = L.PsiTable(cell)
for (s1, s2) in ((Fr(0), Fr(1, 250)), (Fr(3, 50), Fr(1, 5))):
    R = L.Rect(Fr(21, 25), cell, s1, s2); m = R.masses(Psi)
    for t in (cell/4000, (cell + 1)/4000):
        L1 = 2*(math.log(L.N0) + float(s1)); sigU = (1 + h)/2 + t*L1/4; sigV = (1 - h)/2 + t*L1/4
        n = np.arange(1, L.N0 + 1, dtype=np.float64); ln = np.log(n)
        md = {d: math.prod(math.exp(t*math.log(p)**2/4) for p in L.PF[d]) for d in L.DIVS}
        eF = np.exp(t*ln**2/4 - sigU*ln)
        c = np.zeros(L.N0)
        for d in L.DIVS:
            k = np.arange(1, L.N0//d + 1); c[k*d - 1] += L.MU[d]*md[d]*d**(-sigU)*eF[k - 1]
        assert c[1:].min() > -1e-15
        MF = np.zeros(L.KF); MR = np.zeros(L.KRMAX - L.KRMIN + 1)
        np.add.at(MF, np.floor(ln[1:]*20).astype(int), c[1:])
        for d in L.DIVS[1:]:
            k0 = L.N0//d; nn = n[k0:]; np.add.at(MF, np.floor((ln[k0:] + math.log(d))*20).astype(int), md[d]*d**(-sigU)*eF[k0:])
        eV = math.exp(-h*L1/2)*np.exp(t*ln**2/4 - sigV*ln)
        Atot_refl = 0.0
        for d in L.DIVS:
            Atot_refl += np.sum(md[d]*d**(-sigU)*eV)
        Nm = R.Nmax
        if Nm > L.N0:
            na = np.arange(L.N0 + 1, Nm + 1, dtype=np.float64); la = np.log(na)
            eFa = np.exp(t*la**2/4 - sigU*la); eVa = math.exp(-h*L1/2)*np.exp(t*la**2/4 - sigV*la)
            for d in L.DIVS:
                np.add.at(MF, np.floor((la + math.log(d))*20).astype(int), md[d]*d**(-sigU)*eFa)
                Atot_refl += np.sum(md[d]*d**(-sigU)*eVa)
        ok = MF <= m['MF']*(1 + 1e-12)
        nzb = m['MF'] > 0
        r = MF[nzb]/m['MF'][nzb]
        print(f's [{float(s1)},{float(s2)}] t {t}: F-bins direct<=bound all: {bool(ok.all())}; ratio direct/bound min {r.min():.4f} mean {np.average(r, weights=m["MF"][nzb]):.5f} max {r.max():.6f}; total F direct {MF.sum():.6f} bound {m["MF"].sum():.6f}; refl direct {Atot_refl:.6e} bound {m["MR"].sum():.6e}', flush=True)
