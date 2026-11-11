COUNTEREXAMPLE: n = 1,565,734,083, t = g_n = 568,009,639.698978187581685, μ = 3.968215684906503 (certified); two more: n = 2,565,007,238, t = 906,089,969.752693821200281, μ = 3.944872238693607, and n = 3,315,629,581, t = 1,155,454,823.107585085442780, μ = 3.922880048686967 (certified)

# JerbyFals — Jerby's Gram-point repulsion conjecture is false

Orca agent jerbyfals, 11 Nov 2026, 3 h. Every computation single-threaded, ≤ 62 MB (/usr/bin/time -l), one at a time.

**Statement** (TeX re-read). 2406.06548 (Nonlinearity 2024) conj:rep, l. 313 and l. 777: "Assume g_n is a bad Gram point with good consecutive neighbours g_{n±1}. Then |Z'(g_n)| > 4|Z(g_n)|". Gram's law (−1)^n Z(g_n) > 0, θ(g_n) = nπ (l. 164); isolated = middle point of a Gram block with N = 2 (l. 769). Same bound: 2310.14415 "Repulsion G-B-G", l. 266 (μ(g_n) = |Z'/Z| > 4), Conjecture G-B-G l. 753 ("corrupt" points never isolated); restated in 2610.21952 l. 323 as holding "to height of up to g_n of index n ~ 2.5·10^7".

**Result.** At n = 1,565,734,083 (odd) the neighbours are good, g_n is bad, and |Z'(g_n)| = 73.00385040 < 4|Z(g_n)| = 73.58859114, i.e. μ(g_n) = 3.9682156849 < 4; likewise at n = 2,565,007,238 (|Z'| = 55.17454072 < 4|Z| = 55.94557936, μ = 3.9448722387) and n = 3,315,629,581 (|Z'| = 50.23508562 < 4|Z| = 51.22265784, μ = 3.9228800487). Each Gram block [g_{n−1}, g_{n+1}] provably contains no zero: a Rosser-rule violation, not one of Jerby's "confined pairs". All three found by a resonance-targeted scan of t ∈ [1.02·10^8, 1.20·10^9] (n from 2.54·10^8, where the exhaustive scan ended, to 3.46·10^9); the 64 log-uniform windows (t ∈ [10^8, 5·10^11]) alone show none.

## 1. The counterexamples (certified)

Arb (python-flint, 256 bits; certify_cx.py, logs/certify_cx.log): interval-Newton enclosure of each g_k from θ(g_k) = kπ (θ via the continuous log Γ branch), then ζ(½ + i g_k) and its t-derivative as a 2-term acb_series; Z = Re(e^{iθ}ζ) (|Im| ≤ 2·10⁻⁶⁵). Rechecked: certify_gram.py at 200 bits (logs/certify_R2a.log, μ ∈ [3.96821568490650 ± 3.2·10⁻¹⁵]) and mpmath siegelz at 30 digits.

| k | g_k (radius ≤ 4·10⁻²²) | Z(g_k) | (−1)^k Z(g_k) | status |
|---|---|---|---|---|
| n−2 | 568,009,639.013032385630447117 | −0.01937614934137034 | +0.0194 | good |
| n−1 | 568,009,639.356005286611718227 | 0.2763400589567841 | +0.276 | good |
| n | 568,009,639.698978187581685033 | 18.39714778558964 | −18.40 | **bad** |
| n+1 | 568,009,640.041951088540347534 | 4.076959407739048 | +4.08 | good |
| n+2 | 568,009,640.384923989487705731 | −49.65483666225617 | +49.7 | good |

Z'(g_n) ∈ [73.00385040031974 ± 4.4·10⁻¹⁵], 4|Z(g_n)| ∈ [73.58859114235854 ± 3.0·10⁻¹⁵], μ(g_n) ∈ [3.968215684906503 ± 2.3·10⁻¹⁶]. mpmath: Z(g_n) = 18.3971477855896, Z'(g_n) = 73.0038504003197, Z(g_{n±1}) = 0.276340058957, 4.07695940774. Zero-freeness: Arb enclosures of Z on 31,630 adaptive subintervals of [g_{n−1}, g_{n+1}] all exclude 0, sign + (certify_zerofree.py, logs/certify_zerofree.log). So both Gram intervals of the block are empty; the hump of Z covers g_{n−1}, g_n, g_{n+1}, and g_n sits near its top (Z' small relative to Z).

Second counterexample (certify_cx.py, logs/certify_cx2.log): n = 2,565,007,238, g_n ∈ [906,089,969.752693821200280999797 ± 1.5·10⁻²²]; Z(g_{n−1}) = −4.269228364965860 (good), Z(g_n) = −13.98639483953291 (bad), Z(g_{n+1}) = −0.3215394489475987 (good); Z'(g_n) ∈ [55.17454072188091 ± 4.3·10⁻¹⁵], 4|Z(g_n)| ∈ [55.94557935813164 ± 4.7·10⁻¹⁵], μ ∈ [3.944872238693607 ± 4.9·10⁻¹⁶]; mpmath μ = 3.94487223869. Its block is zero-free too: Arb enclosures on 40,539 subintervals of [g_{n−1}, g_{n+1}], all negative (logs/certify_zerofree2.log).

Third (logs/certify_cx3.log): n = 3,315,629,581, g_n ∈ [1,155,454,823.10758508544277960703 ± 3.2·10⁻²¹]; Z(g_{n−1}) = 4.066892916417593 (good), Z(g_n) = 12.80566445995213 (bad), Z(g_{n+1}) = 0.2817344137619705 (good); Z'(g_n) ∈ [−50.23508562012596 ± 2.2·10⁻¹⁶], 4|Z(g_n)| ∈ [51.22265783980850 ± 3.9·10⁻¹⁵], μ ∈ [3.922880048686967 ± 1.1·10⁻¹⁶]; mpmath μ = 3.92288004869. Zero-free: 49,181 Arb subintervals, all positive (logs/certify_zerofree3.log).

## 2. Mechanism: the low tail is zero-free Gram blocks

If g_n is isolated bad, Z(g_{n−1}), Z(g_n), Z(g_{n+1}) have one sign, so the block holds an even number of zeros (0 or 2 here; 4 never met). Every isolated bad point with μ < 5 or μ/θ' < 0.62 was zero-counted (sign changes on 24 points per Gram interval):

| sample | zero-free blocks (count, min μ) | confined pair: 2 zeros in one interval (count, min μ) |
|---|---|---|
| 64 windows | 147, 4.2534 | 76, 5.3356 |
| resonance scan R1 | 895, 4.0527 | 135, 4.8099 |
| resonance scan R2a | 1,377, 3.9682 | 275, 5.1254 |
| resonance scan R2b | 1,726, 3.9449 | 347, 5.0845 |
| resonance scan R3 | 1,333, 3.9229 | 234, 5.2715 |

Over all runs (distinct n): every isolated bad point with μ < 4.8 (872 points) is a zero-free block, i.e. a Rosser-rule violation; of those with μ < 5, 1,745 are zero-free and 1 is a confined pair (μ = 4.8099). So the extremes lie outside the mechanism Jerby offers (2610.21952 Remark "Confined zeros and Rosser's law": in his 250,000-point range no Rosser violations occur, "every isolated bad Gram point ... is associated with a confined pair"). ConjFals' previous record μ(g_69,784,845) = 4.2155 is also a zero-free block (screening grid). Heuristic: under RH, Z'/Z(t) = π/4 + Σ_γ 1/(t − γ) + O(1/t) (sum over ±γ, paired); a confined pair within s = π/θ' of g_n contributes > 2θ'/π, only partly offset by the empty neighbour interval (pair-plus-lattice: μ ≳ θ'/π, ConjFals); a zero-free block has no near zero, and μ → π/4 when the displaced zeros balance on both sides. A refined conjecture "μ > 4 when the block holds two zeros" survives our data (min 4.81) but is untested beyond it.

## 3. Method and validation

**Scanner** (gram_win.c). Riemann–Siegel main sum + C₀ (truncation ≈ |C₁|(t/2π)^{−3/4} < 10⁻⁶ at t ≥ 10^8); Z' from the differentiated main sum plus R'(t). Precision: at t ~ 5·10^11 the phases θ(t)/2π − t log k/2π are ~10^12 turns, so plain double (ConjFals gram_fast2.c) loses ~10⁻⁴ turns per term. gram_win.c forms every phase relative to a block centre t_c (one block = 1024 Gram points): θ(t_c)/2π and frac(t_c log k/2π) in double-double (log k and log t_c by a double-double exp/Newton), the offset δ = t − t_c (|δ| < 400) and the Taylor terms of θ in double; Gram points from θ(t_c + δ) = nπ by Newton in δ. Phase error ~10⁻¹³ turns. cos/sin(2πy) by a degree-16/17 polynomial on |y| ≤ 1/4 (error < 10⁻¹²). Throughput 2.5–3.4·10^8 terms/s single-threaded (11 MB).

**Validation** (validate.py, logs/validate.log): 20 random Gram points per decade against mpmath siegelz (30 digits; derivative=1 for Z'):

| decade | max abs ΔZ, new | max abs ΔZ', new | max rel Δμ, new | max abs ΔZ, old double | max abs ΔZ', old | max rel Δμ, old |
|---|---|---|---|---|---|---|
| [10^8, 10^9) | 7.2·10⁻⁸ | 2.0·10⁻¹¹ | 2.8·10⁻⁷ | 2.0·10⁻⁶ | 1.4·10⁻⁵ | 1.8·10⁻⁵ |
| [10^9, 10^10) | 1.5·10⁻⁸ | 1.8·10⁻¹² | 2.0·10⁻⁷ | 2.9·10⁻⁵ | 2.9·10⁻⁴ | 6.2·10⁻⁵ |
| [10^10, 10^11) | 1.2·10⁻⁹ | 5.7·10⁻¹³ | 5.9·10⁻⁹ | 2.9·10⁻⁴ | 4.2·10⁻³ | 6.6·10⁻⁴ |
| [10^11, 5·10^11] | 2.3·10⁻¹⁰ | 7.3·10⁻¹² | 3.6·10⁻⁹ | 8.0·10⁻³ | 4.8·10⁻² | 1.7·10⁻² |

The new errors are the C₁ truncation (∝ t^{−3/4}); the old plain-double method degrades ∝ t (relative μ error 1.7% at 10^11). Reproduction of ConjFals' Arb-certified points: μ(g_69784845) = 4.215532 (Arb 4.21553180367418), μ(g_191116406) = 4.368210, μ(g_228357901) = 4.425059, same classification of neighbours; ConjFals had matched Jerby's own indexing and values (μ(g_730119), his first-250,000 minimum at n = 8,569).

## 4. Sampled windows: 64 log-uniform windows, t ∈ [10^8, 5·10^11]

Window j: centre t_j = 10^8·5000^{j/63}, size clamp(10^10/N(t), 10^5, 10^6) consecutive Gram points, N = ⌊√(t/2π)⌋ (make_windows.py; logs/run_L1.log, run_L2a.log, run_L2b.log; 2,349 s CPU). n from 2.48·10^8 to 1.918·10^12; total 28,990,518 Gram points, 3,076,101 isolated bad; **no μ < 4**; smallest μ = 4.253426 (L08, n = 781,947,121, t = 2.948·10^8, Arb 4.25342564992270).

| window | t (start) | Gram points | isolated bad | min μ | at n | min μ/θ' | 4/θ' |
|---|---|---|---|---|---|---|---|
| L00 | 9.9811e+07 | 1,000,000 | 105,669 | 5.3000 | 248068690 | 0.6392 | 0.4825 |
| L01 | 1.1429e+08 | 1,000,000 | 105,594 | 4.4012 | 286688825 | 0.5265 | 0.4786 |
| L02 | 1.3086e+08 | 1,000,000 | 105,876 | 5.2953 | 330627595 | 0.6284 | 0.4747 |
| L03 | 1.4983e+08 | 1,000,000 | 105,518 | 5.1069 | 381244233 | 0.6013 | 0.4709 |
| L04 | 1.7155e+08 | 1,000,000 | 105,632 | 4.8903 | 440620566 | 0.5712 | 0.4672 |
| L05 | 1.9641e+08 | 1,000,000 | 106,038 | 4.5564 | 508700884 | 0.5280 | 0.4636 |
| L06 | 2.2487e+08 | 1,000,000 | 105,622 | 5.2457 | 587605362 | 0.6031 | 0.4600 |
| L07 | 2.5745e+08 | 1,000,000 | 106,202 | 4.6584 | 677391393 | 0.5315 | 0.4564 |
| L08 | 2.9475e+08 | 1,000,000 | 105,921 | 4.2534 | 781947121 | 0.4816 | 0.4529 |
| L09 | 3.3744e+08 | 1,000,000 | 106,223 | 4.5660 | 902981076 | 0.5130 | 0.4495 |
| L10 | 3.8631e+08 | 1,000,000 | 105,915 | 5.6406 | 1042138924 | 0.6290 | 0.4461 |
| L11 | 4.4226e+08 | 1,000,000 | 106,237 | 4.9470 | 1202057770 | 0.5475 | 0.4427 |
| L12 | 5.0631e+08 | 1,000,000 | 106,034 | 5.1728 | 1386409650 | 0.5683 | 0.4394 |
| L13 | 5.7963e+08 | 1,000,000 | 105,850 | 4.9283 | 1599742061 | 0.5374 | 0.4362 |
| L14 | 6.6357e+08 | 972,956 | 103,247 | 4.8982 | 1845960970 | 0.5302 | 0.4330 |
| L15 | 7.5966e+08 | 909,361 | 96,417 | 5.2998 | 2129391128 | 0.5695 | 0.4299 |
| L16 | 8.6966e+08 | 849,922 | 89,826 | 4.9319 | 2456589911 | 0.5262 | 0.4268 |
| L17 | 9.9558e+08 | 794,369 | 84,642 | 5.7332 | 2833830962 | 0.6073 | 0.4237 |
| L18 | 1.1397e+09 | 742,447 | 79,097 | 5.0032 | 3268619395 | 0.5262 | 0.4207 |
| L19 | 1.3047e+09 | 693,918 | 73,826 | 5.3142 | 3769798207 | 0.5550 | 0.4177 |
| L20 | 1.4936e+09 | 648,562 | 68,906 | 4.9749 | 4347243136 | 0.5159 | 0.4148 |
| L21 | 1.7099e+09 | 606,170 | 64,465 | 5.2870 | 5013578613 | 0.5444 | 0.4119 |
| L22 | 1.9574e+09 | 566,549 | 60,095 | 4.5128 | 5781308681 | 0.4615 | 0.4091 |
| L23 | 2.2408e+09 | 529,518 | 56,590 | 6.1021 | 6666357631 | 0.6198 | 0.4063 |
| L24 | 2.5652e+09 | 494,907 | 52,640 | 4.6459 | 7686968138 | 0.4686 | 0.4035 |
| L25 | 2.9365e+09 | 462,559 | 49,284 | 5.2906 | 8862647022 | 0.5301 | 0.4007 |
| L26 | 3.3616e+09 | 432,325 | 46,097 | 5.1268 | 10217992112 | 0.5102 | 0.3981 |
| L27 | 3.8483e+09 | 404,067 | 43,167 | 5.2515 | 11779857236 | 0.5191 | 0.3954 |
| L28 | 4.4054e+09 | 377,656 | 39,875 | 4.8995 | 13580056320 | 0.4811 | 0.3928 |
| L29 | 5.0431e+09 | 352,971 | 37,530 | 4.9962 | 15654280031 | 0.4874 | 0.3902 |
| L30 | 5.7731e+09 | 329,900 | 35,057 | 6.0051 | 18044585613 | 0.5819 | 0.3876 |
| L31 | 6.6088e+09 | 308,337 | 32,838 | 4.8399 | 20798712107 | 0.4660 | 0.3851 |
| L32 | 7.5655e+09 | 288,183 | 30,730 | 5.0704 | 23972525400 | 0.4850 | 0.3826 |
| L33 | 8.6607e+09 | 269,346 | 28,645 | 5.3423 | 27628823507 | 0.5077 | 0.3802 |
| L34 | 9.9144e+09 | 251,741 | 26,831 | 5.4523 | 31841716824 | 0.5149 | 0.3777 |
| L35 | 1.1350e+10 | 235,287 | 24,984 | 5.8657 | 36695428216 | 0.5504 | 0.3753 |
| L36 | 1.2993e+10 | 219,908 | 23,273 | 4.9670 | 42286995624 | 0.4631 | 0.3730 |
| L37 | 1.4873e+10 | 205,534 | 21,769 | 4.5923 | 48728345068 | 0.4255 | 0.3706 |
| L38 | 1.7026e+10 | 192,100 | 20,377 | 5.4897 | 56148532776 | 0.5055 | 0.3683 |
| L39 | 1.9491e+10 | 179,543 | 19,069 | 6.0805 | 64695782172 | 0.5564 | 0.3660 |
| L40 | 2.2313e+10 | 167,808 | 18,051 | 5.5271 | 74541173402 | 0.5027 | 0.3638 |
| L41 | 2.5543e+10 | 156,840 | 16,705 | 5.3343 | 85881153789 | 0.4822 | 0.3616 |
| L42 | 2.9240e+10 | 146,588 | 15,590 | 5.5751 | 98942461770 | 0.5009 | 0.3594 |
| L43 | 3.3473e+10 | 137,007 | 14,505 | 6.0600 | 113985393439 | 0.5412 | 0.3572 |
| L44 | 3.8318e+10 | 128,051 | 13,770 | 6.6024 | 131310150788 | 0.5861 | 0.3551 |
| L45 | 4.3865e+10 | 119,682 | 12,772 | 6.0880 | 151262320574 | 0.5372 | 0.3529 |
| L46 | 5.0215e+10 | 111,859 | 11,784 | 5.7098 | 174239156239 | 0.5008 | 0.3509 |
| L47 | 5.7484e+10 | 104,547 | 11,160 | 5.9603 | 200698626316 | 0.5197 | 0.3488 |
| L48 | 6.5806e+10 | 100,000 | 10,648 | 7.0245 | 231167345481 | 0.6089 | 0.3467 |
| L49 | 7.5331e+10 | 100,000 | 10,682 | 5.6482 | 266251651545 | 0.4868 | 0.3447 |
| L50 | 8.6236e+10 | 100,000 | 10,529 | 6.0438 | 306649363305 | 0.5178 | 0.3427 |
| L51 | 9.8720e+10 | 100,000 | 10,654 | 5.8765 | 353163484998 | 0.5006 | 0.3407 |
| L52 | 1.1301e+11 | 100,000 | 10,640 | 5.5780 | 406718438665 | 0.4725 | 0.3388 |
| L53 | 1.2937e+11 | 100,000 | 10,611 | 5.6884 | 468377919698 | 0.4791 | 0.3369 |
| L54 | 1.4810e+11 | 100,000 | 10,601 | 4.9223 | 539366073986 | 0.4122 | 0.3350 |
| L55 | 1.6954e+11 | 100,000 | 10,507 | 6.5964 | 621091688584 | 0.5493 | 0.3331 |
| L56 | 1.9408e+11 | 100,000 | 10,682 | 5.9098 | 715175804165 | 0.4894 | 0.3312 |
| L57 | 2.2217e+11 | 100,000 | 10,749 | 5.6603 | 823483821018 | 0.4661 | 0.3294 |
| L58 | 2.5433e+11 | 100,000 | 10,612 | 6.1567 | 948162346775 | 0.5041 | 0.3275 |
| L59 | 2.9115e+11 | 100,000 | 10,713 | 5.8306 | 1091681272202 | 0.4748 | 0.3257 |
| L60 | 3.3329e+11 | 100,000 | 10,604 | 7.3097 | 1256882627309 | 0.5920 | 0.3240 |
| L61 | 3.8154e+11 | 100,000 | 10,635 | 6.4989 | 1447036370160 | 0.5235 | 0.3222 |
| L62 | 4.3677e+11 | 100,000 | 10,570 | 6.0132 | 1665904891170 | 0.4817 | 0.3205 |
| L63 | 5.0000e+11 | 100,000 | 10,719 | 7.0309 | 1917816785061 | 0.5602 | 0.3187 |

**Trend.** Per-window min μ = 1.31 + 0.42 log₁₀ t (r = 0.70) and min μ/θ' = 0.714 − 0.019 log₁₀ t (r = −0.42); the rise of min μ is a sample-size effect (windows above 6·10^10 hold 10^4 isolated bad points against 10^5 below 4·10^8). Pooled tail fractions per isolated bad point (analysis_L.txt) are flat in absolute μ and rise in normalised μ:

| height bin | isolated bad | P(μ < 5) | P(μ < 6) | P(μ < 7) | P(μ/θ' < 0.55) | P(μ/θ' < 0.65) |
|---|---|---|---|---|---|---|
| 10^8 | 846,403 | 8.3·10⁻⁶ | 1.9·10⁻⁴ | 2.4·10⁻³ | 4.7·10⁻⁶ | 4.0·10⁻⁵ |
| 10^8.5 | 904,391 | 7.7·10⁻⁶ | 8.7·10⁻⁵ | 1.4·10⁻³ | 7.7·10⁻⁶ | 6.4·10⁻⁵ |
| 10^9 | 504,903 | 7.9·10⁻⁶ | 9.1·10⁻⁵ | 7.3·10⁻⁴ | 3.6·10⁻⁵ | 1.6·10⁻⁴ |
| 10^9.5 | 320,770 | 9.4·10⁻⁶ | 7.5·10⁻⁵ | 3.7·10⁻⁴ | 5.3·10⁻⁵ | 2.0·10⁻⁴ |
| 10^10 | 159,818 | 1.3·10⁻⁵ | 1.3·10⁻⁴ | 2.6·10⁻⁴ | 1.3·10⁻⁴ | 3.1·10⁻⁴ |
| 10^10.5 | 106,504 | < 10⁻⁵ | 4.7·10⁻⁵ | 1.8·10⁻⁴ | 1.0·10⁻⁴ | 2.8·10⁻⁴ |
| 10^11 | 85,115 | 1.2·10⁻⁵ | 1.3·10⁻⁴ | 2.3·10⁻⁴ | 2.0·10⁻⁴ | 4.3·10⁻⁴ |

So the absolute lower edge of μ is roughly height-independent over 10^8–10^11 (P(μ < 5) ≈ 10⁻⁵), while the normalised tail grows 40-fold; the brief's normalised-margin argument does not by itself locate counterexamples. The first counterexample came after ≈ 10^8 effective isolated bad points (resonance scan, §5: 1.8·10^8 in range × ≈ 70% capture); the exhaustive 2.7·10^7 (ConjFals) and the 3·10^6 sampled here were too few.

## 5. Resonance-targeted scan (largest-|ζ| regions)

**Method** (gram_win.c -reso). Proxy F(t) = Σ_{p ≤ 29} p^{−1/2} cos(t log p) + ½ Σ_{p ≤ 7} p^{−1} cos(2t log p) (≈ log|ζ|) on the grid t = T + 0.1j (phases by rotation, re-anchored in double-double every 4096 steps); bins of length 2; every bin whose maximum exceeds the 99% quantile of bin maxima (pilot over the first 2·10^5 bins, F_thr ≈ 2.67) gets an exact check of the 7 Gram points nearest its argmax (plus neighbours). Proxy validation (reso_test.py, capture.py): of the exhaustively known window candidates inside the R1 range it re-found 2/4 with μ < 4.6, 7/9 with μ < 5, 13/19 with μ < 6, at 1% of the bins: ≈ 70% capture at ≈ 2% of the cost. Checked regions are 56% bad Gram points (18% in the windows); largest |Z(g_n)| met: 114 (R1), 124 (R2a), 129 (R2b).

| scan | t range | bins checked | Gram points | isolated bad | μ < 4 | min μ (n) | CPU |
|---|---|---|---|---|---|---|---|
| R1 | [1.0227·10^8, 4.0227·10^8] | 1,401,506 | 9,344,642 | 880,512 | 0 | 4.052662 (452,435,915) | 381 s |
| R2a | [4.0227·10^8, 7.0227·10^8] | 1,387,503 | 9,231,772 | 819,845 | **1** | **3.968216 (1,565,734,083)** | 446 s |
| R2b | [7.0227·10^8, 1.0023·10^9] | 1,379,466 | 9,165,366 | 791,711 | **1** | **3.944872 (2,565,007,238)** | 647 s |
| R3 | [1.0023·10^9, 1.2023·10^9] | 935,756 | 6,211,944 | 530,514 | **1** | **3.922880 (3,315,629,581)** | 485 s |

R1 starts where ConjFals' exhaustive scan ends (n = 2.54·10^8). Tail counts in the checked regions: R1 μ < 4.25: 11, < 4.5: 59, < 5: 421; R2a μ < 4: 1, < 4.25: 7, < 4.5: 55, < 5: 465; R2b μ < 4: 1, < 4.25: 13, < 4.5: 50, < 5: 493; R3 μ < 4: 1, < 4.25: 9, < 4.5: 49, < 5: 354. Three counterexamples in ≈ 3.4·10^8 isolated bad points in range (≈ 2.4·10^8 effective): roughly one per 10^8 at t ~ 5·10^8–1.2·10^9; largest |Z(g_n)| in R3: 130.

## 6. Certified minima

Arb 200 bits (certify_gram.py) with mpmath recheck; all listed points isolated bad (certified), radii ≤ 5·10⁻¹⁵, mpmath agrees to all printed digits (logs/certify_R1.log, certify_R2a.log, certify_R2b.log, certify_R3.log, certify_cx*.log):

| n | t | μ (Arb) | source |
|---|---|---|---|
| 3,315,629,581 | 1.1555·10^9 | **3.92288004868697** | R3 |
| 2,565,007,238 | 9.0609·10^8 | **3.94487223869361** | R2b |
| 1,565,734,083 | 5.6801·10^8 | **3.96821568490650** | R2a |
| 1,253,300,061 | 4.6026·10^8 | 4.04501154627085 | R2a |
| 452,435,915 | 1.7604·10^8 | 4.05266198533572 | R1 |
| 2,991,193,487 | 1.0669·10^9 | 4.06471489945353 | R3 |
| 1,376,902,138 | 5.0303·10^8 | 4.06565880463454 | R2a |
| 2,313,783,317 | 8.2185·10^8 | 4.08846491622984 | R2b |
| 3,080,776,521 | 1.0976·10^9 | 4.09065314376269 | R3 |
| 349,212,915 | 1.3796·10^8 | 4.12524290051701 | R1 |
| 1,714,309,720 | 6.1885·10^8 | 4.12552448596725 | R2a |
| 2,879,069,975 | 1.0302·10^9 | 4.13427590259362 | R3 |
| 1,374,792,884 | 5.0231·10^8 | 4.13886052244661 | R2a |
| 2,484,227,201 | 8.7905·10^8 | 4.15600195495161 | R2b |
| 2,139,464,629 | 7.6313·10^8 | 4.15962454966579 | R2b |
| 921,299,273 | 3.4418·10^8 | 4.16844820792381 | R1 |
| 1,944,554,913 | 6.9719·10^8 | 4.16870084932230 | R2a |
| 499,340,611 | 1.9318·10^8 | 4.20406150476804 | R1 |
| 781,947,121 | 2.9483·10^8 | 4.25342564992270 | R1 and window L08 |

All 64 window minima certified in Arb as well (SKIP_MPMATH=1, logs/certify_win.log): all isolated bad, none below 4; screened 6-digit values deviate ≤ 1.1·10⁻⁷ (relative). The other 32 certified points with μ < 4.3 (4.16–4.30) are in the logs.

## 7. Novelty

- arXiv monthly listings Oct–Nov 2026, 7 math categories (.scratch/lst.py '\bGram\b|Jerby|viscosity|Rosser|repulsion|Hardy function'): Gram-titled papers 2610.21952 (Jerby; restates the conjecture as open, data n ≤ 250,000), 2610.22107, 2610.39856 (moments at Gram points); none tests or refutes the bound. arXiv API (web_search approach arxiv, after two HTTP 429): top hit 2610.21952.
- Zenodo "Gram points", "Jerby OR bad Gram OR Gram point repulsion": only Jerby's dataset 23390338 (15 Oct 2026, 250,000 Gram points).
- GXL evidence API, one claim ("at every isolated bad Gram point μ > 4"): 5 supporting and 6 "refuting" passages, all from Jerby's own 2310.14415, 2406.06548, 2610.21952 (statements of the bound; small μ only for non-isolated points or the Davenport–Heilbronn function).
- Our files: RECORDS.md and ConjFals/REPORT.md (bound holds for n < 2.54·10^8); no NOVELTY file on Gram points.
- Verdict: no prior counterexample; the bound had been tested only to n ~ 2.5·10^7 (author) and n < 2.54·10^8 (ConjFals). The zero-free-block origin of the low tail is not stated in Jerby's three papers (rg for Rosser/zero-free/no zeros: only the 2610.21952 remark quoted in §2).

RH relevance: N8 (Gram/Rosser counting), N13 (Jerby's proposed mechanism). N8: the extremes of μ are Rosser-rule violations (zero-free Gram blocks of length 2). N13: the universal repulsion bound at isolated bad Gram points is false; all three counterexamples are zero-free blocks, so the confined-pair version stays open. RH itself is untouched.

## Files

gram_win.c (scanner: windows, -check, -reso), validate.py, make_windows.py, reso.c and reso_test.py (proxy calibration), capture.py, analyze.py, cert_list.py, certify_gram.py (copied from ConjFals), certify_cx.py, certify_zerofree.py; logs/: validate.log, windows_*.txt, run_L1.log, run_L2a.log, run_L2b.log, analysis_L.txt, reso_R1.log, reso_R2a.log, reso_R2b.log, reso_R3.log, reso_test_L24.log, certify_R1.log, certify_R2a.log, certify_R2b.log, certify_R3.log, certify_win.log, certify_cx.log, certify_cx2.log, certify_cx3.log, certify_zerofree.log, certify_zerofree2.log, certify_zerofree3.log (+ *.time: /usr/bin/time -l, max RSS ≤ 62 MB).
