
import re,os
os.chdir('/Users/adriansager/playground/riemann')
T='.scratch/postrh/tex/'
checks=[
('COROLLARIES.md',63,'using GRH only through Lemma lem:grh'),
('COROLLARIES.md',70,'using GRH only through θ(x,χ) ≪ x^{1/2}(log kx)²'),
('COROLLARIES.md',71,'proves it under GRH for L(s, χ_D), D < 0'),
('COROLLARIES.md',72,"GRH enters only through Montgomery–Vaughan's Σ_{n≤x} χ(n) ≪ √p log log p"),
('COROLLARIES.md',76,'GRH enters only through Lemma GRH (l.824–831: Σ_{t≤T}(t/m) ≤ T^{1/2}m^{o(1)}'),
('COROLLARIES.md',79,"using it only through Darbar–Maiti's lemma"),
('COROLLARIES.md',79,'as written its proof covers only x ≤ X^{δ/2}'),
('COROLLARIES.md',39,'GRH enters only through Prop. charsum (l. 148–173); Bebel\'s Dirichlet half-plane gives the character sums'),
('COROLLARIES.md',40,'its Assumption (l. 79) is exactly a zero-free half-plane Re s > 1 − c₂ for every Dirichlet L-function'),
('COROLLARIES.md',49,'assume only that no Dirichlet L-function'),
('COROLLARIES.md',50,'Goldfeld–Huang 2017 Thm 1.2 needs only no zero of Λ(s, χ) with 1/2 < σ < 1'),
('COROLLARIES.md',54,'Near-null spectrum of the windowed Weil form under RH'),
('COROLLARIES.md',81,'GRH for L(s,χ) implies lim_{λ→0} D_χ(λ)² log(1/λ)'),
('Riemann/Deep/Exp48/REPORT.md',28,'they need a small split prime (ERH via Bach–Sorenson)'),
('discoveries/nyman-beurling-continuous-rate.md',1,'under RH alone'),
('DISCOVERIES.md',10,'Under RH alone'),
('discoveries/li-remainder.md',7,'Under RH limsup |R(n)|/(n log n)^{1/4} ≥ 0.289'),
('discoveries/li-remainder.md',46,'Pointwise only (0.202 + o(1))·√n log n/√(log log n) is proved, under RH'),
('discoveries/odd-robin.md',3,'RH ⟺ σ(n) < (e^γ/2) n log log n for every odd n > N'),
('discoveries/nyman-beurling-lp.md',3,'ζ(s) ≠ 0 on the open half-plane Re s > 1/p'),
('discoveries/nyman-beurling-lp.md',21,'L(s, χ) ≠ 0 on Re s > 1/p'),
('COROLLARIES.md',46,'zero-freeness on Re s > 1/p is the only hypothesis'),
('discoveries/small-gap-two-thirds.md',7,'Under RH: λ > 0.6878 (Soundararajan), 0.6553 (Wu), 0.6039 (Bui–Goldston–Milinovich–Montgomery), 1/2 (Durkan–Pearce-Crump, arXiv 2610.28527, Thm B)'),
('discoveries/linear-combinations-line.md',7,'Under GRH alone: Bombieri–Hejhal (Duke 1995, Theorem C), 34.5%'),
('discoveries/linnik-numeric.md',9,'Polylogarithmic non-residue bounds: GRH (Bach)'),
('discoveries/littlewood-half-plane.md',5,'π²/(12e^γ)'),
('discoveries/prime-primitive-root.md',12,'Under GRH: ĝ(p) ≤ ((8/5)(2^{ω(p−1)} − 1) log p)² for p ≥ 10⁹ (McGown–Treviño–Trudgian, arXiv 1508.05182)'),
('discoveries/unit-group-cover.md',5,'under GRH a sketch covers every q with polylogarithmically many primes'),
('COROLLARIES.md',36,'under GRH 2 log² n (Bach 1990)'),
('COROLLARIES.md',35,'GRH gives 2 log log q + O(1)'),
('COROLLARIES.md',61,'L(s, f) = L(s, ρ_f) (Deligne–Serre) is, by Brauer induction,'),
('COROLLARIES.md',80,'if and only if RH holds and every nontrivial zero of ζ is si'),
('discoveries/chebotarev-numeric.md',9,'(4 log Δ_L + 2.5 n_L + 5)²'),
('discoveries/least-primes-forms.md',11,'A = 1'),
('COROLLARIES.md',73,'ERH enters only through effective Chebotarev (l.485–490)'),
('TREE.md',15,'zero-free Re s > 11/12 for finite-order Hecke L-functions of Q(√−3)'),
(T+'1801.05914.tex',129,'One has $\\Lambda \\geq 0$.'),
(T+'1302.5018.tex',73,'Assuming RH we have'),
(T+'2610.28527.tex',62,'Assume the Riemann Hypothesis for the Riemann zeta-function.'),
(T+'2610.28527.tex',76,'Assume the Riemann Hypothesis for $L(s,\\chi)$. Then $\\lambda>3.63846723$.'),
(T+'0908.2008.tex',100,'Assume RH.'),
(T+'2010.13307.tex',123,'Assume the Riemann Hypothesis.'),
(T+'1305.4618.tex',87,'Assume the Riemann Hypothesis is true, and let $k \\geq 0$ be fixed.'),
(T+'1410.7015.tex',60,'which is implied by the Riemann hypothesis'),
(T+'1708.04122.tex',368,'Assume the Riemann hypothesis. Then, for $x\\ge 4$, there is always a prime number in the interval'),
(T+'0705.0723.tex',62,'Assume RH.'),
(T+'math_0008177.tex',174,'Robin\'s criterion states that the'),
(T+'math_0202141.tex',47,'Riemann hypothesis is equivalent to the stronger statement that'),
(T+'1902.07321.tex',277,'he showed that the RH is equivalent to the hyperbolicity of all Jensen polynomials'),
(T+'math_0404394.tex',1096,'We can now state a general version of Li\'s criterion.'),
(T+'2005.05142.tex',119,'$\\Lambda_F \\geq 0$ for every $F \\in\\mathcal{S}^\\sharp$.'),
(T+'1309.3595.tex',218,'Assume GRH. Let $q>3$'),
(T+'1309.3595.tex',576,'Suppose that GRH holds for $L(s,\\chi)$.'),
(T+'1309.3595.tex',158,'Assume GRH.  If $q\\ge 5$ is prime, the least quadratic non-residue'),
(T+'1309.3595.tex',135,'Assume GRH. Let $q\\geq 3000$ be an integer, and let $H$ be a proper'),
(T+'1309.3595.tex',270,'requires the assumption of GRH for $L(s,\\chi_{-q})$'),
(T+'1508.05182.tex',45,'Assume GRH.'),
(T+'1508.05182.tex',86,'Let $\\chi(n)$ denote a non-principal Dirichlet character modulo $p$.'),
(T+'1803.02056.tex',254,'Assume the extended Riemann hypothesis. Then there exists a split'),
(T+'math_0412262.tex',555,'Hooley \\cite{Hooley1} proved in 1967'),
(T+'math_0305114.tex',64,'Assume that the $L$-functions of all the curves $E_{r,s}$ satisfy the Riemann Hypothesis.'),
(T+'1211.5191.tex',249,'Assume the Riemann hypothesis and that all the zeros of $\\zeta(s)$ are simple.'),
(T+'1810.08843.tex',278,'0.6792'),
(T+'2412.17199.tex',86,'Assuming the Generalised Riemann Hypothesis for Dirichlet $L$-functions (GRH)'),
(T+'2412.17199.tex',127,'Assume GRH. Then there is an effectively computable $N_0$'),
(T+'2304.13605.tex',423,'for the splitting field of $X^2 - a$'),
(T+'2509.24119.tex',174,'subject to the Extended Riemann Hypothesis (ERH) for quadratic Dirichlet ch'),
(T+'2504.15658.tex',43,'Assume GRH. Then, $B < 2.1594$.'),
(T+'2610.28527.tex',45,'with $Q=1$ corresponding to the Riemann zeta-function'),
(T+'2005.05142.tex',113,'there is a real number $\\Lambda_F$ such that all the zeros of $\\xi^F_t$ lie on the critical line'),
(T+'1810.08843.tex',281,'0.6845'),
(T+'2010.13307.tex',139,'0.5,6.992,0.252'),
]
norm=lambda s: re.sub(r'\s+',' ',s)
bad=0
for f,ln,frag in checks:
    L=open(f,encoding='utf8',errors='replace').read().split('\n')
    lo,hi=max(0,ln-4),min(len(L),ln+3)
    win=norm(' '.join(L[lo:hi]))
    ok = norm(frag) in win
    if not ok:
        # whole-file search to report actual line
        hits=[i+1 for i,l in enumerate(L) if norm(frag)[:40] in norm(l)]
        bad+=1; print('MISS',f.split('/')[-1],ln,'|',frag[:70],'| found at',hits[:5])
print('checked',len(checks),'missing',bad)
