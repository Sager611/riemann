Nazarov/Turán statement (quoted verbatim in Friedland–Yomdin, arXiv:1107.0039, Thm. 1.1):
“Let p(t)=sum_{k=0}^{m} c_k e^{lambda_k t}, c_k,lambda_k in C. Let B subset R be an interval, and let Omega subset B be a measurable set. Then
 sup_B |p| <= e^{mu_1(B)·max|Re lambda_k|}·(c mu_1(B)/mu_1(Omega))^m·sup_Omega |p|
where mu_1 is Lebesgue measure on R and c>0 is an absolute constant.”
Thus user’s n is m = number of indexed terms minus 1 (m+1 terms), not number of terms. The exponential factor is exactly |B| max_k|Re lambda_k|. The cited source does not state a numerical value for c/A: only “absolute constant”. It states no L2 version of this lemma.

Independent source: O. Friedland, Y. Yomdin, “An observation on the Turán-Nazarov inequality”, arXiv:1107.0039, Thm. 1.1 (also Thm. 1.2 replaces |Omega| by metric span omega_D(Omega)).

Erdelyi: T. Erdélyi, “Inequalities for exponential sums”, arXiv:1602.02315. It defines T_n=sum_{j=1}^n a_j exp(i lambda_j t), real ordered frequencies, and quotes/proves |f(0)| <= (pi/2)n ||f||_{L2[0,1]} for f in E_n^+ (Eq. 1.5; the T_n version is Eq. 1.4 with asymptotic constant 1.8401...). This gives endpoint sup-square <= (pi^2/4)n^2 times the normalized L2 average after translation/scaling; general sup on one fixed interval needs the standard translated-interval argument, so report C n^2, not a claimed universal numerical C unless hypotheses are specified.

Chen–Kane–Price–Song, “Fourier-sparse interpolation without a frequency gap”, arXiv:1609.01361, gives the exact quoted interpolation theorem instead: for k-Fourier-sparse x*=sum_{j=1}^k v_j e^{2 pi i f_j t} on [0,T], ||x~−x*||_T <= C||g||_T + delta||x*||_T (probability >=1−2^{-Omega(k)}), with poly(k,log(1/delta)) log(FT) samples. The paper’s polynomial special case is ||P−Q||_T^2 <= C||g||_T^2, O(d) samples (Thm. 1.2); it is not itself a Nazarov lemma citation.
