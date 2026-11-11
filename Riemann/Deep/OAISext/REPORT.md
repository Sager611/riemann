BLOCKED: background-QoS build runs 25.6 s/module (60 modules; unclamped 7.5 s); 2,099 modules left → ≈ 15 h > 13.1 h left; cause: background page eviction (800 MB cap) re-faults each lean's 4–5 GB of olean imports (≈ 300k faults/module); FundamentalDomain.lean hit 632 MB > 600 MB.

# oaisext: unconditional balanced sextic fourth moment over Q(√−3), build blocked (11 Nov 2026, 22:33–23:40 UTC)

Brief: .scratch/orca/oaisext-brief.md. State: .scratch/oaisext/STATE.md. Logs: logs/ (build.log, samples.log, chunk_00{1,2}.log trimmed to lake's Built lines with per-module times, diffcheck.log). Lean paths below are relative to Riemann/Deep/QRH/OAI/NumberTheory/DirichletL/ unless stated.

## Short answer
- **Refactor: written and statically checked, not elaborated.** 240 files under Riemann/Deep/OAISext/ (238 forked + BetaEff.lean + Wrapper.lean). Diff = betaEff substitution + 2 lines of PrimeGlobal (logs/diffcheck.log).
- **Build: blocked by the compute policy, not by Lean.** 134 of the 1,993 upstream modules are built (oleans kept in .lake/build/lib/lean/Riemann/Deep/QRH, 180 MB); the fork has not started.
- **Both stop rules fired:** projected time (≈ 15 h > 13.1 h left; the upstream part alone needs 13.2 h) and the memory cap (Eisenstein/FundamentalDomain.lean, 632 MB > 600 MB, 23:24:25 UTC).

## Build measurements
Single process: LEAN_NUM_THREADS=1 makes lake run one lean at a time (ps) and makes lean single-threaded (4-task test: real 3.06 s = user 2.83 s; unset: 1.14 s). nice 19 throughout; lake footprint 190 MB.

| run | window (UTC) | modules | s/module | lines/min | lean footprint | page faults |
|---|---|---|---|---|---|---|
| utility QoS clamp (bg removed; before the orchestrator's order) | 22:49:15–22:58:33 | 74 (avg 705 lines) | 7.5 wall; lake's own: mean 7.6, median 6.8 | 5,610 | 398 MB (1 sample) | — |
| inherited background QoS (taskpolicy -b) | 22:58:50–23:24:25 | 60 (avg 673 lines) | 25.6 wall; lake's own: mean 24.7, median 19.0; by 20s: 21.3, 21.4, 31.2 | 1,577 | 13 of 14 samples > 450 MB, max 632 MB | median 304,837 per process (max 339,589) ≈ 4.7–5.2 GB of 16 KB pages |

Controlled test (Detector/Compensation.lean, 47 lines, 5 Mathlib imports; lake env lean -j 1; runs in order):

| policy | real (s) | user (s) | sys (s) | page faults |
|---|---|---|---|---|
| background | 36.5 / 41.9 / 80.3 | 4.4 / 4.4 / 4.9 | 8.1 / 8.5 / 15.5 | 115,106 / 158,845 / 357,661 |
| background + taskpolicy -g default (I/O policy) | 36.8 / 58.8 | 4.7 / 4.1 | 8.6 / 11.3 | 151,095 / 263,385 |
| utility clamp, bg removed | 17.6 / 5.4 / 4.6 | 2.5 / 2.3 / 2.2 | 4.2 / 2.4 / 2.3 | 106,526 / 216 / 237 |

- **Cause.** vm.vm_page_background_mode = 1, vm_page_background_target = 50,000 pages (800 MB): pages faulted by background-QoS threads are reclaimed first once the background queue exceeds 800 MB. Each lean maps 2–5 GB of Mathlib + QRH oleans and walks them to build its environment, so its own pages are evicted while still needed, and every module re-reads them from SSD (sys time 2–3× user). The I/O policy does not matter (-g default: no gain); the eviction does. Machine at the time: 31 GB used, 0.2–0.5 GB free, swap 11.9 of 13.3 GB, load 7–8.6.
- **Projection.** Per-module cost is dominated by import loading (47 lines cost 37–80 s), so the per-line model (3.7 h) does not apply. Per module: 2,099 × 25.6 s = 14.9 h; at the last-20 rate (31.2 s) 18.2 h. The rate was still rising as modules got deeper. Left at the stop: 13.1 h (12.1 h before the 13 h write-up). The brief's 11 h checkpoint (CertifiedExistence by 09:33) needs 1,859 upstream × 25.6 s = 13.2 h. Same model unclamped: 2,099 × 7.5 s = 4.4 h.
- **Memory.** Under background QoS most leans exceed the 450 MB target (AnnularRecursion 558, TupleWeights 529, UnrestrictedCharacters 527, LogDensity 526, UniformTransfer 534, GlobalBinLoss 528 MB). Eisenstein/FundamentalDomain.lean (746 lines; in the closure) reached 632 MB → driver stop. Its unclamped footprint is unmeasured.

## What is ready (unbuilt)
- **Closure** (.scratch/oaisext/closure.py, closure.json). CertifiedExistence: 2,231 files; downstream of Moments/PrimeGlobal: 238 (32,151 lines; 183 Energy, 55 Moments); upstream: 1,993 (25 maximal modules: targets_up.txt). By construction no file outside the 238 imports one inside, so the fork never imports an original downstream file.
- **Fork** (fork.py → Riemann/Deep/OAISext/OAI/NumberTheory/DirichletL/{Energy,Moments}/, prefix Riemann.Deep.OAISext.).
  - Imports of forked modules rewritten; header line 3 added.
  - Substitutions: 102 HeckeZeroSupremum.beta, all inside (51/100:ℝ)≤· or 2*·−1≤κ; 45 bare beta in the 14 files that open HeckeZeroSupremum. Total 147 in 64 files.
  - Untouched: 6 local binders named beta (Energy/ReferenceDeletionCapacity, plain_kappa_eq/lower) and 7 accesses of the field OriginalData.beta.
- **PrimeGlobal, the only proof change.**
  - + import Riemann.Deep.OAISext.BetaEff.
  - zeroMaximum_le_beta: sentinel case unchanged (exact hbeta, now 51/100 ≤ betaEff); zero case (zero_re_le_beta …).trans beta_le_betaEff.
  - ray_prime_global_squared: beta_le_one → betaEff_le_one; the exponent becomes (betaEff − 1/2)r by the substitution.
- **BetaEff.lean.** betaEff := max beta (51/100); fiftyOne_le_betaEff (le_max_right), beta_le_betaEff (le_max_left), betaEff_le_one (max_le beta_le_one).
- **Wrapper.lean.**
  - zeroAt_unconditional: terminal_certificate at κ = 1 (hκ0 by norm_num, hbeta := fiftyOne_le_betaEff, hκ from betaEff_le_one), then certified_terminal, zero clause. Auxiliaries: slot type Empty, M = H = ⊤ (instance NeZero ⊤), W = 0 on [1, 1], lo = hi = 0.
  - zero_clause_unconditional: the same with ZeroAt unfolded.
- **Static checks.**
  - diffcheck.py: after undoing header, imports and substitution, 237 files equal their originals; PrimeGlobal differs in exactly the 2 lines above.
  - No def, structure, abbrev or instance in the 238 files mentions beta. Every definition in the final statement (NaturalState, plainEnergy, energy, Profiles, control, diagonalControl, internalQ) is upstream and unmodified; ZeroAt (Energy/Bands) is forked verbatim.
  - rg sorry|admit|native_decide|axiom|skipKernelTC|implemented_by|extern|unsafe: 0 hits in the 240 files (and in the 238 originals).
- **Not done:** elaboration of the 240 files, #print axioms, the fidelity table against Lemma 18.1 (o.l. 12552–12558) and Liu's Remark 12.18, the reviewer.

## The statement Wrapper.lean would prove (unbuilt)
```lean
theorem zero_clause_unconditional (a b radial Bmask L Mcap ε:ℝ)
    (ha:0<a)(haUpper:a≤1/4)(hb:1≤b)(hrad:0<radial)(hmask:0≤Bmask)(hMcap:0<Mcap)(hε:0<ε):
    ∃degree:ℕ,∃S:Finset (ℕ×ℕ),∀η₀:Character,∀Q:Ideal O,
      internalQ Q η₀≠0→internalQ Q η₀≠⊤→internalQ Q η₀≤Ideal.span {(72:O)}→
      ∃C:ℝ,0<C ∧ ∀ᶠZ:ℝ in atTop,1<Z ∧
        ∀s:NaturalState Z Bmask radial,s.fixedModulus=internalQ Q η₀→s.width≤Mcap→
        ∀p:Profiles a b,∀t X₁ X₂:ℝ,0<X₁→0<X₂→X₁≤Z^L→X₂≤Z^L→
          s.plainEnergy p t X₁ X₂≤C*diagonalControl s.radial.profile*(p.control S)^2*
            (1+‖t‖)^degree*Z^(s.width+ε)
```
- **In words** (definitions as read in LiuSext/REPORT.md §(i)). Fix 0 < a ≤ 1/4 ≤ 1 ≤ b, bΦ > 0, Bmask ≥ 0, M > 0, any L and ε > 0. Then there are degree and S such that, for every character η₀ and every ideal Q with Q ∩ mod η₀ nonzero, proper and divisible by 72, there is C > 0 with the following for all large Z.
  - For every state: character η of modulus norm ≤ Z^{characterWidth}; puncture of norm ≤ Z^{Bmask}; weight Φ supported in (−∞, bΦ] at scale Z^{rowWidth}; kept rows z ≠ 0 and not fixed-inducing; width = rowWidth + characterWidth ≤ M.
  - For every profile pair W₁, W₂ supported in [a, b], every t, and 0 < X₁, X₂ ≤ Z^L:
  - Σ'_{z∈O, kept} |(X₁X₂)^{−1/2} S_z(X₁) S_z(X₂)|² Re Φ(N z/Z^{rowWidth}) ≤ C·diag(Φ)·control_S(W)²·(1+|t|)^degree·Z^{width+ε}, where S_z(X) = Σ_I η(I) (m⁶z/I)_6 N(I)^{it} W(N I/X).
  - No hypothesis on β_*.
- **Balanced case.** η trivial, rowWidth = 1, X₁ = X₂ = Z^{1/2}, t = 0: Σ_{0<Nz≤Z, kept} |Z^{−1/4} S_z(Z^{1/2})|⁴ ≪ Z^{1+ε}.

## Resume
- **Remaining work:** 1,859 upstream + 240 fork modules. The fork may still need fix-ups.
- **Commands.**
  - rm .scratch/oaisext/STOP .scratch/oaisext/MEMSTOP; nohup bash .scratch/oaisext/driver.sh & (targets.txt: BetaEff, fork PrimeGlobal, Wrapper, then the 25 upstream tops; chunks of 1,800 s, samples every 30 s, stop at 600 MB).
  - Then #print axioms OAI.SevenEighths.OAISextWrapper.zero_clause_unconditional.
- **Options.** Each needs the user, because of the background-QoS rule and the 600 MB cap:

| option | time | cost |
|---|---|---|
| exception: one unclamped, single-threaded, nice 19 lean | ≈ 4.4 h measured, + fork fix-ups | uses a P-core |
| background QoS, cap waived for FundamentalDomain | ≈ 15 h, rising | ≈ 5 GB of SSD reads per module, ≈ 10 TB in total |
| restore a port olean archive | minutes | none on disk (the port's 1.4 GB of oleans were deleted earlier) |

- **Note for every Lean agent in this repo:** the same eviction hits every Lean build under taskpolicy -b.

RH relevance: N7 (input of number-field zero-density; Liu's B_2 uses the F ⊃ μ12 analogue).

Suggested DISCOVERIES line: none until a green build and a reviewer verdict. Then: "Balanced sextic fourth moment over Q(√−3), Σ_{Nz≤Z}|Z^{−1/4}S_z(Z^{1/2})|⁴ ≪ Z^{1+ε}, Lean-checked unconditionally (OpenAI Lemma 18.1 without 51/100 ≤ β_*)."

Files: Riemann/Deep/OAISext/{BetaEff.lean, Wrapper.lean, logs/ (incl. STATE.md), scripts/ (closure.py, fork.py, diffcheck.py, driver.sh, rate.py, fit.py, memsum.py, scan*.py; closure.json, fork_manifest.json, targets*.txt)}. The 238 forked files are not committed (orchestrator): regenerate them with scripts/fork.py.
