# Balazard–de Roton 0812.1689 intro scout

- Source: `main.tex` from arXiv e-print; extracted document lines below are source-relative.
- **Continuous D(lambda):** the introduction does not define or mention a continuous distance `D(lambda)`, dilations `[lambda,1]`, or any corresponding RH bound. Therefore no log-log exponent for `D` is stated in this paper's introduction.
- **Discrete setup:** `d_N = dist_H(chi, Vect(e_1,...,e_N))`, with `H=L^2(0,infty;t^-2 dt)`; RH iff `d_N -> 0` (lines 76–87).
- Under RH, Báez-Duarte (BD2003) is quoted: `d_N^2 << (log log N)^(-2/3)` (lines 111–115). This is the only prior log-log exponent stated for discrete `d_N`.
- **Balazard–de Roton theorem:** for every `delta>0`, RH implies `d_N^2 <<_delta (log log N)^(5/2+delta) (log N)^(-1/2)` for `N>=3` (lines 116–122; abstract gives `5/2+o(1)`).
- Short quote: “Nous améliorons ce résultat dans le présent travail.” (line 116).
- Burnol’s known lower bound: `d_N^2 >= (C+o(1))/log N`, `C=sum_rho m(rho)^2/|rho|^2` (lines 89–97).
- The paper says BDBLS2000 conjecture equality in the displayed lower-bound consequence `d_N^2 >= (2+gamma-log(4pi)+o(1))/log N`; it says this would imply RH and simple zeta zeros (lines 99–109).
- **BBLS continuous asymptotic:** no formula `D(lambda)^2 log(1/lambda) -> sum 1/|rho|^2` appears, and it is not called open here. “BBLS” itself is absent; only citation BDBLS2000 occurs (bibliography line 1162).
- Caution: do not attribute the continuous `D(lambda)` claims or their log-log exponents to this introduction; they likely belong to a related paper/source, not 0812.1689.
