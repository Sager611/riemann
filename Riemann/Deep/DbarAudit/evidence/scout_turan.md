**Source:** arXiv 1107.0039 e-print; gzip payload extracted to `/tmp/dbaraudit_turan/tex` (source line numbers below).

**Theorem 1.1 (rendered; TeX label `turan`):**
- `/tmp/dbaraudit_turan/tex:72-78`: `\begin{theorem}[\cite{N}]\label{turan}`; `Let $p(t) = \sum_{k=0}^m c_k e^{\lambda_k t}$ be an exponential polynomial, where $c_k ,\lambda_k \in \C$. Let $B \subset \R$ be an interval, and let $\Omega \subset B$ be a measurable set.`
- `/tmp/dbaraudit_turan/tex:75-78`: `$$` / `\sup_{B} |p| \le e^{\mu_1(B) \cdot \max|\Re \lambda_k|} \cdot \left( \frac{c \mu_1(B)}{\mu_1(\Omega) }\right)^{m} \cdot \sup_{\Omega} |p|` / `$$` / `where $\mu_1$ is the Lebesgue measure on $\R$ and $c>0$ is an absolute constant.`

**Interpretation:** Sum has `m+1` terms (`k=0,...,m`), not n/n+1 notation; exponent is `m` = number of terms minus 1. The real-part factor is exactly `e^{|B| max|Re lambda_k|}` (`mu_1(B)=|B|`). `c` is an absolute positive constant. E/Ω is explicitly measurable subset of interval B. Both displayed suprema are ordinary `sup`, not essential sup; no `ess sup` occurs in the TeX.

**Requested rg hits:**
- `/tmp/dbaraudit_turan/tex:175-179`: `There is a version of Tur\'an-Nazarov inequality for quasipolynomials in one or several variables due to A. Brudnyi ... this result gives an important` / `information for a wider class of quasipolynomials.`
- `/tmp/dbaraudit_turan/tex:401-405`: `A quasipolynomial is a finite sum` / `p(z) = \sum_{j=1}^k p_j(z) e^{f_j\cdot z}` / `where $p_j \in \C[z_1,\ldots,z_n]$ are polynomials in $z$ of degrees $d_j$. The degree of $p$ is $m = \deg p = \sum_{j=1}^k (d_j + 1)$.`
- `/tmp/dbaraudit_turan/tex:393-394`: `In this section we consider the version of Tur\'an-Nazarov inequality for quasipolynomials in one or several variables due to A. Brudnyi ...`

**Polynomial-coefficient extension:** No sentence in this TeX says the 1D inequality extends to `\sum p_k(t)e^{\lambda_k t}` with polynomial coefficients. The closest explicit definition is the multidimensional quasipolynomial above; its degree counts `\sum_{j=1}^k(d_j+1)` (line 405), not a 1D theorem extension. `multiplicit` has no matches. `Nazarov` hits include theorem citation at lines 70-78 and bibliography/lemma, but no polynomial-coefficient extension sentence.
