| module | wall | footprint (time -l) | memcap peak (top MEM) |
|---|---|---|---|
| E1Lean Defs | 15 s | 243 MiB | 170 MB |
| E1Lean Detector | 15 s | 252 MiB | 204 MB |
| E1Lean AlgebraH | 19 s | 252 MiB | 233 MB |
| E1Lean AlgebraC | 41 s | 278 MiB | 264 MB |
| E1Lean Numeric | 16 s | 265 MiB | 248 MB |
| E1Lean Kappa | 20 s | 248 MiB | 188 MB |
| E1Lean Arch | 23 s | 249 MiB | 214 MB |
| E1Lean Parity | 49 s | 286 MiB | 274 MB |
| E1Lean Bounds2 | 38 s | 273 MiB | 251 MB |
| Trig | 17 s | 258 MiB | 189 MB |
| Exp | 74 s | 254 MiB | 240 MB |
| Num | 126 s | 273 MiB | 248 MB |
| Arch | 45 s | 255 MiB | 233 MB |
| Window | 27 s | 247 MiB | 207 MB |
| Parity | 9 s | 239 MiB | 0 MB (exited between samples) |
| DetTrig | 9 s | 255 MiB | 235 MB |
| Det | 30 s | 267 MiB | 249 MB |
| Check | 18 s | 238 MiB | 144 MB |

All compiles of the session (38 runs, peaks.log + verify.log): footprint ≤ 286 MiB, top MEM ≤ 274 MB (both E1Lean's Parity); E1PLean files alone: footprint ≤ 273 MiB, top MEM ≤ 262 MB; one failed compile (Num, a reduced-literal mismatch with Exp's c_hi, 125 s, 262 MB; fixed by writing c_hi unreduced). Scripts (Python, exact Fractions): ≤ 1 s each, below memcap's first sample. Every run under memcap 290 MB; all below the 300 MB cap.
