# NOVELTY — Koiran's HN ∈ AM without GRH (agent koiran, 2026-11-06, 11:26–12:30 UTC)

Claims checked: (A) HN ∈ AM (hence Π₂ᵖ) unconditionally, given Cureton–Thorner; (B) DIM ∈ AM unconditionally (K97);
(C) the descendants: HNZ/HNP ∈ AM (2408.13027), Schubert vanishing in coAM / AM ∩ coAM (2412.02064, 2504.03004), GW
non-vanishing ∈ AM (2508.15715), ideal primality in Σ₃ᵖ ∩ Π₃ᵖ (2503.20071), resultant vanishing ∈ AM (1210.1451),
Noether normalization in PH (1209.5993), isotropic-space problems in PH (1904.03950); (E) Bürgisser's characteristic-zero
transfer (VP_ℂ = VNP_ℂ ⇒ #P ⊆ P/poly) and its users (Koiran–Perifel cs/0610009 and 0706.1477, Fournier–Perifel–de Verclos
1304.5910, Bürgisser 2406.06217 Thm (2), 2309.07612) without GRH; (D) the numeric Lemma A
(degree-one prime ideals of every K, x ≥ e^{40}: ≥ (x/log x)[0.2772 − 0.4953(log D + 6.43)x^{−1/32} − 0.7233 e x^{−1/2}]).

## Verdict

| claim | verdict | source |
|---|---|---|
| (A) HN ∈ AM unconditionally | **KNOWN** as a remark | Andrews–Garg–Oliveira–Schost, arXiv 2611.05289v1 (3 Nov 2026), sections/introduction.tex l.18, footnote: Cureton–Thorner "claims to prove an effective version of the prime ideal theorem, which would yield an unconditional version of Koiran's result". No proof or parameters given. |
| (B) DIM ∈ AM | KNOWN implicitly (same remark; l.27 says K97 needs the same effective prime ideal theorem); one-line consequence of (A) via K97's randomised reduction | 2611.05289 l.27; K97 |
| (C) descendants | routine corollaries of (A); none stated anywhere we found. Schubert vanishing: superseded, SV ∈ coRP unconditionally since Sep 2025 (Pak–Robichaux 2509.16467, l.1141–1149) ⊆ AM ∩ coAM | TeX of each paper; arXiv versions not updated since 2025 |
| (E) Bürgisser's transfer and users without GRH | **not stated anywhere we found**; still called conditional by Bürgisser on 20 Oct 2026 (2606.25121 l.294–296, "assuming the Generalized Riemann Hypothesis. This relies on the same methods as employed in [Koiran 96]") and by 2611.07972 on 6 Nov 2026 (sections/introduction.tex l.735–738, "known only under the Generalized Riemann Hypothesis"); routine once Lemma B (count version) is available | TeX of 2406.06217, 2606.25121, 2611.07972, cs/0610009, 0706.1477, 1304.5910 |
| (D) numeric Lemma A | minor: the counting form of our chebexp machinery (Riemann/Deep/ChebExp) with L = K; FSW 2610.18491 give the same shape with an unspecified constant | count_cert.py |

Prior hypotheses weaker than GRH that already suffice: Rojas math/0301111 (2003, DZH or GIPIT ⇒ HN ∈ AM; Cureton–Thorner
implies GIPIT for K ≠ ℚ, DERIVATION §9); Rojas–Zhu 1803.04104 (2018, MDZH ⇒ DIM ∈ P^{NP^{NP}}). The existence version
(one small prime with a root, Kuperberg's lemma) is Baldwin–Sivek 2611.06705 (4 Nov 2026), who also note Agol knew it.

## Searches
- arXiv monthly listings 2026-10 and 2026-11 (.scratch/lst.py), titles: regex "Nullstellensatz|Arthur|AM\b|coAM|Koiran|
  polynomial hierarchy" over cs.CC, math.AG, math.NT (math.AG 2026-11 failed once, fetched in the second run), and the
  same plus "Riemann" over math.AG, cs.CC, cs.SC, math.AC, math.CO: 24 + 39 hits, none on HN in AM or on GRH removal
  (closest: 2611.07280 Nullstellensatz proof complexity; 2610.21039, 2610.14417 algebraic Nullstellensätze).
- arXiv API: one OR-query over au:Koiran, au:Kuperberg, au:Robichaux, au:Agol, au:Pak_Igor, au:Sivek, au:Garg_Abhibhav
  (40 newest): 2611.06705, 2509.16467 (SV ∈ coRP, unconditional), 2508.15715, 2504.03004, 2412.02064, nothing else
  relevant; au:Rojas_J_M / au:Zhu_Yuyu: nothing after 2018; abs: Nullstellensatz AND (GRH OR Riemann), "Hilbert
  Nullstellensatz" AND AM, "polynomial hierarchy" AND Riemann, coAM AND GRH, Gromov-Witten AND (vanishing OR complexity):
  only the known papers (2508.15715, 2412.02064, 2408.13027, 1803.04104, math/0501388, math/0005029, math/9811088).
  Version dates (API): 2503.20071v1 (Mar 2025), 2508.15715v1 (Aug 2025), 2408.13027v3 (Oct 2024), 2412.02064v3 (Apr 2025),
  2504.03004v1, 2509.16467v1, 1803.04104v1: no 2026 revision.
- GXL Evidence API, 2 claims ("HN is in AM unconditionally, without GRH"; "a zero-free half-plane for all Dedekind zeta
  functions suffices to place HN in AM"; 4 calls, the first 2 lost to a local parse error): supporting passages only from
  Rojas math/0301111 (conditional on DZH/GIPIT) and 2508.15715 (HNPE ∈ AM[4], conditional); refuting and neutral passages
  all say "assuming GRH" (2408.13027, 2504.03004, 2412.02064, 2602.17904, 2006.07613, 1210.1451, 2002.09472, 1209.5993,
  2610.42138 of 30 Oct 2026). They led to 2611.05289 (neutral list), whose footnote is the pre-emption.
- For (E): GXL Evidence API claim "VP = VNP over ℂ collapses the nonuniform Boolean classes, without GRH": 1 support
  (0907.2850, irrelevant: VNP-completeness of the permanent), 8 refutations all "under GRH" (2611.07972, 2406.06217,
  1604.00795, 1304.5910, 2309.07612, cs/0610009, 0710.0360); arXiv cs.CC/cs.SC listings Oct–Nov 2026, title regex
  Valiant|VNP|VP|Cook|transfer|arithmetic circuit|algebraic circuit|constant-free|counting hierarchy|Bürgisser: 4 hits
  (2610.00075, 2610.00644, 2610.42138, 2610.01440), none on the transfer; 2610.00644 and 2610.42138 read: GRH only quoted;
  Hexagon search Valiant, VNP, Burgisser, "arithmetic circuit": no hits; Zenodo ("VNP" OR "Valiant") AND ("Riemann" OR
  "GRH"): 5 unrelated; "Cook" AND "Valiant" AND "hypothesis": 3 unrelated. Bürgisser's book (2000) and Bürgisser 2009 are
  not on arXiv and were not read.
- Hexagon (hexagonmath.org/search): Koiran, Arthur-Merlin, Kuperberg: no hits; Nullstellensatz: 2610.00327 (hypercube
  Nullstellensatz counterexample), 2610.00713; polynomial hierarchy: 2611.00086 (class numbers in FP^{Σ₂}, cites K96 in its
  bibliography only); coNP: 2610.00368; quasi-Riemann: 2610.00178, 2610.00200, 2610.00216, 2610.00655, 2611.00082. TeX of
  2610.00178 ("Unconditional implications of the Quasi-Riemann Hypothesis", GPT survey, 9 Oct 2026), 2611.00086 and
  Cureton–Thorner 2610.00200v2: no Koiran, Nullstellensatz, Arthur–Merlin or Kuperberg.
- Zenodo API (AND queries): "Nullstellensatz" AND ("Arthur-Merlin" OR "polynomial hierarchy" OR "Arthur Merlin"): 0;
  "Nullstellensatz" AND ("Riemann" OR "GRH"): 1 unrelated; "Cureton" AND "Thorner": 23547020, 23373237 (number theory);
  "quasi-Riemann" AND ("complexity" OR "Arthur" OR "coNP"): 0; Schubert AND GRH AND (AM OR coAM): 0; "Koiran": 10 hits,
  none relevant (Finnish "koiran").
- HAL API: Nullstellensatz AND (Arthur OR "polynomial hierarchy" OR Riemann): the K96 notice (hal-04649685, metadata
  only) and 2408.13027; "quasi-Riemann": nothing relevant.
- Papers using Cureton–Thorner, read in TeX: 2610.18491 (FSW), 2610.37179 (sparse irreducibility in PP; no HN, no AM),
  2611.02280 (primality with one binomial congruence), 2611.06705 (Baldwin–Sivek, coNP only), 2611.05289 (the footnote).
- Our files: rg -i 'koiran|nullstellensatz' outside .scratch finds only Riemann/Deep/Harvest4/HARVEST.md (the source of
  this brief). No NOVELTY file covers Koiran.

## Closest prior work
| work | date | proves | leaves open |
|---|---|---|---|
| Koiran, J. Complexity 12 (DIMACS TR 96-27) | 1996 | HN ∈ AM under GRH | unconditional |
| Koiran, FOCS | 1997 | DIM ∈ AM under GRH | unconditional |
| Rojas, math/0301111 | 2003 | HN ∈ AM under DZH or GIPIT (weaker than GRH) | proving DZH or GIPIT |
| Rojas–Zhu, 1803.04104 | 2018 | MDZH ⇒ DIM ∈ P^{NP^{NP}} | AM; MDZH misprinted (false for ℚ) |
| Andrews–Garg–Schost, 2602.17904 | Feb 2026 | HN ∈ CH unconditionally | AM, PH |
| Cureton–Thorner, Hexagon 2610.00200v2 (Lean) | Oct 2026 | ζ_K ≠ 0 on Re s > 31/32, every K | applications |
| Frei–Sofos–Widmer, 2610.18491 | Oct 2026 | effective Chebotarev with x^{1−δ} error under QRH(δ) | numeric constant |
| Andrews–Garg–Oliveira–Schost, 2611.05289 | 3 Nov 2026 | footnote: CT "would yield an unconditional version of Koiran's result" | details; DIM and descendants not named |
| Baldwin–Sivek, 2611.06705 | 4 Nov 2026 | Kuperberg's coNP results modulo QRH(ε) (one small prime with a root) | the AM count |
| Pak–Robichaux, 2509.16467 | Sep 2025 | SV ∈ coRP unconditionally | the GW case |
| Bürgisser, book (2000) Ch. 4; survey 2406.06217 (2024); 2606.25121 (Oct 2026) | 2000–2026 | VP_ℂ = VNP_ℂ ⇒ #P ⊆ P/poly under GRH (constants elimination as in Koiran) | unconditional |
| 2611.07972 | 6 Nov 2026 | low-degree PIT ⇒ circuit lower bounds; recalls the transfer as "known only under GRH" | — |
