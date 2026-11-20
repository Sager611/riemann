
p = "Riemann/Deep/B3Cert/Lean/Arch.lean"; s = open(p).read()
old2 = s[s.index("/-- (H2′)"):s.index("def H2P")]
new2 = "/-- (H2′) SL Thm 1.3 (Zenodo 23134787 v0.03, pdftotext l. 125–129) at every (θ, ν) where Thm 1.2's conclusion holds,\n  θ ∈ (0, 1): θ-free per FULLPROOF §5 l. 111–112, 116 (SL TeX l. 823–928: Lemma 4.2, Thm 1.3). -/\n"
old3 = s[s.index("/-- (H3′)"):s.index("def H3P")]
new3 = "/-- (H3′) SL Thm 5.1 (pdftotext l. 844–845; TeX l. 932–963), θ ∈ (0, 1): θ-free per FULLPROOF §5 l. 112–114, 116\n  (Wang 2609.07918 Thm local-pair, its TeX l. 231–232, 'Let 0 < λ < θ < 1 be fixed'). -/\n"
s = s.replace(old2, new2).replace(old3, new3)
open(p, "w").write(s)
p = "Riemann/Deep/B3Cert/Lean/Window.lean"; w = open(p).read()
w = w.replace("(FULLPROOF.md §5 l. 106–111; SL l. 730–758)", "(FULLPROOF.md §5 l. 106–110; SL TeX l. 730–758)")
w = w.replace("/-- (HCW) SL l. 730–758, compact windows (profile-free; FULLPROOF §5 l. 106–111)", "/-- (HCW) SL TeX l. 730–758, compact windows (profile-free; FULLPROOF §5 l. 106–110)")
w = w.replace("(FULLPROOF §§2–4: Lemma B, Lemma R, branches)", "(FULLPROOF §§2–4: Lemma B, Lemma R, branches; §5 l. 108–109)")
open(p, "w").write(w)

