# dpcblean NOTES
## 15:05 PROOF.md digest (Riemann/Deep/DPCBeyond/PROOF.md)
- Thm U': kappa>0 density slope; phi in (max{2/(2+k),1/2},1) => c(phi)>0 count >= c T log T. kappa=8/7-delta => phi>7/11.
- What changes: theta in (1/phi-1, min{k/2,1}); construction needs 0<phi*theta<1<phi(1+theta)<2, theta<1. Step 5a -> Lemma Z + M; 6e-6f -> 6e'; Step 7 Hoelder -> 7'.
- Lemma Z (iii): |e^{sum z_i u_i}| <= e^{max(Re z_i)^+ sum u_i} <= (X/a)^{max(Re z_i)^+}; (iv) X^{max c_i} <= sum_i X^{c_i}; sum over rho_1..rho_m X^{max} prod kappa <= m W On^{m-1}.
- M2: int <= meas{n>=1}^{1-1/w} ||n_m||_{pw}^p <= (|B_a|/L)^{1-1/w} (CT)^{1/w}; |B_a|<=C T L e^{-k a} => << T e^{-k a(1-1/w)}.
- M3: ||Off||_p << T^{1/p} sum_a e^{theta(a+1) - k a(1-1/w)/p}, converges if p theta < k (w large). q<k/theta.
- M4: count <= |B_a|^{1-1/w}(sum N^w)^{1/w} << (|m|+2) TL e^{-k a (1-1/w)}; exponent k(1-1/w)^2.
- 6e'(d): 1/v+1/u+1/w=1, u>1, 2u theta<k. H_{T,a} << T (a+2)^C [e^{2theta(a+1)} sigma_a^{1-1/w} + sigma_a^{1/v}]; sigma_a = min(2eps+o(1), C e^{-k a}).
- 6e'(e): w with k(1-1/w)>2theta, A=c log(1/eps), 0<c<(1-1/w)/(2theta): sum_{a<=A} << (A+2)^C e^{2theta(A+1)} eps^{1-1/w} ; sum_{a>A} << (A+2)^C e^{-(k(1-1/w)-2theta)A}; sigma^{1/v} series DCT, majorant (a+2)^C e^{-k a/v}.
- 7': q in (2,k/theta), s with 2/q+1/s<1: delta1 T <= C T^{2/q+1/s} meas(E)^{1-2/q-1/s} => meas(E) >= c2 T.
- Order: phi; theta; P,q,f,eta,A; q in (2,k/theta), s; u in (1,k/(2theta)), w>u/(u-1) with k(1-1/w)>2theta, 1/v = 1-1/u-1/w; r,r'; eps; T.
