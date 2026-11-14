/-!
# Round 479 — node N8 (critical-line proportions and simple zeros), round 29, threshold and input question on our Theorem E1: no paper proves a mollified or twisted second moment of ζ on windows T^θ with θ ≤ 1/2 and a fixed mollifier length ν > 0, so E1's binding Taylor step (T/H²)^J stays the method's limit; the classical linear-Q (Heath-Brown) simple-zero detector lies in E1's class but gives κ_lin(0.1325) = −0.115 < 0 (best linear Q −0.035; positive only for ν > 0.1925, θ > 0.590), so Wang's pair theorem and the parity step cannot be dropped; Santibáñez-Leal's detector counts distinct odd-multiplicity zeros; no pre-emption of E1 or E1Arb; the relative entropy method, the Eichler–Selberg trace formula and regular model sets transfer no tool

* Node: N8 (29th round; R475 before). Round type: threshold and input question on a fresh theorem of ours (lesson 213); texts
  pre-extracted (Santibáñez-Leal Zenodo 23134787 v0.03 by pdftotext; Bui–Conrey–Young 1002.4127 and Wang 2609.07918 TeX).
* [threshold] E1's θ > 1/2 is bound by the Taylor remainder O(T^{1+ν+ε}(T/H²)^J) (SantiAudit/E1.md l.53, l.87–88); search: SL
  (1/2 < θ < 1) and 2104.10243 (a ≥ 5/7 with ν < 2(7a − 5)/7; a > 1/2 + ν) — no paper reaches θ ≤ 1/2 with a fixed ν > 0
  (scout timed out after its check; partial): HOLDS.
* [input] Bui–Conrey–Young (1002.4127, TeX l.203–219, 283–287): Q real, Q(0) = 1, Q(x) + Q(1 − x) constant; 'by choosing Q(x) to
  be a linear polynomial, one obtains a lower bound on the percent of simple zeros' (κ* ≥ 0.4058 with Q = 1 − 1.03x, R = 1.12).
  Feng (1003.0059) states 41.28% critical, no simple-zero theorem.
* [input] Santibáñez-Leal (l.755–760, 843–869): the localized detector gives O ≥ N − 2N_Ṽ − O(L), O = distinct critical zeros
  of odd multiplicity, N_Ṽ with multiplicity including σ = 1/2; Wang's pair theorem and the parity inequality turn O/N into S/N;
  no obstruction to a direct simple count is stated.
* [rate, orchestrator] the direct route: Q = 1 − x is in E1's class (Q(y) + Q(1 − y) ≡ 1); with the optimal P = sinh(μx)/sinh μ,
  c = 1/2 + √(AB)·coth(ν√(B/A)), A = ∫₀¹e^{2Rv}Q², B = ∫₀¹e^{2Rv}(RQ + Q′)²: κ_lin(1/2) = 0.34857 (Levinson's regime), κ_lin(0.1325) =
  −0.11501, κ_lin(0.125) = −0.13176; best Q = 1 − ax at ν = 0.1325: a = 1.2, −0.0347; κ_lin > 0 only for ν > 0.192474, i.e.
  θ > 0.589966: Wang's pair theorem and the parity step are needed on (1/2, 0.59]; clean negative (mpmath, 30 digits, R scan).
* [check] pre-emption: SL concept 22984154 still at v0.03; arXiv 'simple zeros' / 'critical line' with 'short intervals': no
  new id below θ = 0.527 or above 0.0238 at θ = 0.51; Zenodo phrase query: nothing: HOLDS.
* [explore] relative entropy method (Gibbs variational principle) → N1: needs a positive measure; Weil's explicit-formula
  distribution is signed (test ν = δ₀ − δ₁): no transfer.
* [explore] Eichler–Selberg trace formula → N3: no map from Hecke traces to Li coefficients (λ₁ = 0.0230957 > 0 against
  τ(2) = −24): no transfer.
* [explore] regular model sets (cut-and-project) → N13: their counting function is linear, ζ's is (T/2π) log T: no transfer.
-/

namespace Fayed.R479
end Fayed.R479
