# Burnol math/0103058 audit

Source: `/tmp/dbaraudit_burnol/paper.tex`, made from `https://arxiv.org/e-print/math/0103058` (gzip single TeX; read-only fetch).

## (1) D(lambda), space, dilations, constraints
- Exact setup: `paper.tex:169-175`: “Let $K = L^2(]0,\infty[,dt)$ (over the complex numbers), let $\chi$ be the indicator function of the interval $]0,1]$... Let $0<\lambda<1$ and let $\cB_\lambda$ be the sub-vector space of $K$ consisting of the finite linear combinations of the functions $t\mapsto \rho({\theta\over:t})$, for $\lambda\leq\theta\leq1$.”
- Exact distance: `paper.tex:189-193`: “Let us write $D(\la)$ for the Hilbert-space distance $\inf_{f\in\cB_\lambda} \| \chi -f \|$. We have $$\liminf_{\la\to0} D(\la)\sqrt{\log({1\over\la})} \geq \sqrt{\sum_\rho {1\over |\rho|^2}}$$”
- Thus this D uses $L^2(0,\infty)$, dilations/parameters $\theta\in[\lambda,1]$, and finite linear combinations. Burnol imposes no displayed constraint such as $\sum c_k\theta_k=0$ in this definition.
- Later equivalent reformulation uses a different-looking model: `paper.tex:380-386`: “Let $\cC_\la$ ($0<\la<1$) be the sub-vector space of $\HH^2$ of linear combinations of the contractions $D_\theta(A)$ for $\la\leq\theta\leq1$... The quantity $D(\lambda)$ ... is thus the Hilbert space distance between $\chi_1(t)$ and $\cC_\la$.”

## (2) Burnol main theorem and logical status
- Main theorem exact: `paper.tex:203-206`: “We have: $$\liminf_{\la\to0} D(\la)\sqrt{\log({1\over\la})} \geq \sqrt{\sum_\rho {m_\rho^2\over |\rho|^2}}$$”
- Immediately before it Burnol writes `paper.tex:197-201`: “If the Riemann Hypothesis fails this result is true but trivial as the left-hand side then takes the value $+\infty$. So we will assume that the Riemann Hypothesis holds. The sum on the right-hand side is over all non-trivial zeros $\rho$ of the zeta function, counted \emph{only once} independently of their multiplicities $m_\rho$.”
- Therefore the displayed improved theorem is proved in the RH-assumed regime; Burnol does not state the sum as “only zeros with $\Re\rho=1/2$”. Under the explicit RH assumption, all nontrivial zeros are on that line. Multiplicity is squared: $m_\rho^2$, not $m_\rho$.
- Proof confirms the square: `paper.tex:663-670`: “the top-left element of the inverse matrix is $m_\rho^2$... [the squared projection norm] is asymptotically equivalent ... to $\sum_{\rho\in R} {m_\rho^2\over |\rho|^2}\over\log({1\over\la})$.”

## (3) BBLS result as Burnol states it
- `paper.tex:189-194` attributes to Báez-Duarte, Balazard, Landreau and Saias: $\liminf_{\lambda\to0}D(\lambda)\sqrt{\log(1/\lambda)}\geq\sqrt{\sum_\rho1/|\rho|^2}$.
- `paper.tex:197-201` says the statement is true but trivial if RH fails (left side $+\infty$), and defines its sum over all nontrivial zeros, each counted once independently of multiplicity.
- Abstract: `paper.tex:83-85`: “We slightly improve the lower bound of B\'aez-Duarte, Balazard, Landreau and Saias in the Nyman-Beurling formulation of the Riemann Hypothesis as an approximation problem.”
