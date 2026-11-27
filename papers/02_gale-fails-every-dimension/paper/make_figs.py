# /// script
# dependencies = ["numpy", "matplotlib"]
# ///
"""Figures of gale.tex from the repository data (points7.txt, points8.txt, points7_sym69.txt, family_params.txt).

Run from anywhere:  uv run --no-project --with numpy --with matplotlib python make_figs.py
Writes figs/disk.pdf, figs/sets.pdf, figs/hgraph.pdf next to this file.
"""
from fractions import Fraction
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, Rectangle, Wedge

HERE = Path(__file__).resolve().parent
DATA = HERE.parents[2] / "Riemann" / "Deep" / "Gale7"
OUT = HERE / "figs"

plt.rcParams.update({
    "font.family": "serif", "font.serif": ["cmr10"], "mathtext.fontset": "cm", "axes.unicode_minus": False,
    "axes.formatter.use_mathtext": True, "font.size": 8, "pdf.fonttype": 42, "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2,
})
# orbit -> (colour, marker); P/A share a hue (row blocks), Q/B share a hue (column blocks); Okabe-Ito palette
STYLE = {"P": ("#0072B2", "o"), "A": ("#56B4E9", "^"), "Q": ("#D55E00", "s"), "B": ("#E69F00", "v"),
         "C": ("#009E73", "D")}
SIZE = {"P": 13, "Q": 11, "A": 11, "B": 11, "C": 9}


def read_points(name):
    rows = []
    with (DATA / name).open() as fh:
        for raw in fh:
            line = raw.split("#")[0].strip()
            if line:
                rows.append([float(Fraction(s)) for s in line.replace(",", " ").split()])
    return np.array(rows)


def params(n):
    with (DATA / "family_params.txt").open() as fh:
        for line in fh:
            w = line.split("#")[0].split()
            if w and int(w[0]) == n:
                return int(w[1]), int(w[2]), [float(x) for x in w[3:11]]
    raise KeyError(n)


def simplices(X, m1, rA, rB):
    """the unit simplex vectors u_i, v_j of a family member: its points with x'' = 0 (resp. x' = 0) are r_A u_i (r_B v_j)"""
    n1, n2 = np.linalg.norm(X[:, :m1], axis=1), np.linalg.norm(X[:, m1:], axis=1)
    U = X[(n2 < 1e-6) & (np.abs(n1 - rA) < 1e-6), :m1] / rA
    W = X[(n1 < 1e-6) & (np.abs(n2 - rB) < 1e-6), m1:] / rB
    return U, W


def family(n, U, W):
    """V_n: points and labels (orbit, i, j)"""
    m1, m2, (rA, rB, a, b, c, d, e, f) = params(n)
    pts, lab = [], []
    for i in range(m1 + 1):
        pts.append(np.r_[rA * U[i], np.zeros(m2)]); lab.append(("P", i, None))
    for j in range(m2 + 1):
        pts.append(np.r_[np.zeros(m1), rB * W[j]]); lab.append(("Q", None, j))
    for orb, (x, y) in (("A", (a, -b)), ("B", (-c, d)), ("C", (e, f))):
        for i in range(m1 + 1):
            for j in range(m2 + 1):
                pts.append(np.r_[x * U[i], y * W[j]]); lab.append((orb, i, j))
    return np.array(pts), lab


def label_by_nearest(X, Y, lab):
    d = np.linalg.norm(X[:, None, :] - Y[None, :, :], axis=2)
    return [lab[k] for k in d.argmin(axis=1)], d.min(axis=1)


def polygon_map(U, phase):
    """2 x m matrix with orthonormal rows sending the unit simplex vectors U (p x m) to a regular p-gon"""
    p, m = U.shape
    ang = phase + 2 * np.pi * np.arange(p) / p
    Z = np.sqrt(2 / m) * np.c_[np.cos(ang), np.sin(ang)]
    Pi = (m / p) * Z.T @ U
    assert np.allclose(Pi @ Pi.T, np.eye(2), atol=1e-6)
    return Pi


def load():
    sets = {}
    X8 = read_points("points8.txt")
    m1, _, prm = params(8)
    U, W = simplices(X8, m1, prm[0], prm[1])
    Y8, lab8 = family(8, U, W)
    l8, err8 = label_by_nearest(X8, Y8, lab8)
    assert err8.max() < 1e-7 and len(X8) == 85
    sets[8] = {"X": X8, "lab": l8, "U": U, "W": W, "m1": m1}
    S7 = read_points("points7_sym69.txt")
    m1, _, prm = params(7)
    U, W = simplices(S7, m1, prm[0], prm[1])
    Y7, lab7 = family(7, U, W)
    X7 = read_points("points7.txt")
    l7, _ = label_by_nearest(X7, Y7, lab7)
    assert len(X7) == 84
    sets[7] = {"X": X7, "lab": l7, "U": U, "W": W, "m1": m1}
    return sets


def scatter(ax, P2, lab, emph=None, big=1.0):
    for orb in "CBAQP":
        idx = [k for k, L in enumerate(lab) if L[0] == orb]
        col, mk = STYLE[orb]
        s = SIZE[orb] * (big if orb in "PQ" else 1.0)
        ax.scatter(P2[idx, 0], P2[idx, 1], s=s, c=col, marker=mk, linewidths=0.3, edgecolors="white", zorder=3)
    if emph:
        for k, L in enumerate(lab):
            if L in emph:
                col, mk = STYLE[L[0]]
                ax.scatter(P2[k, 0], P2[k, 1], s=3.2 * SIZE[L[0]], c=col, marker=mk, linewidths=0.6,
                           edgecolors="black", zorder=5)


def clean(ax, lim):
    ax.set_aspect("equal")
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def axes_style(ax, lo, hi, xlabel, ylabel):
    ax.set_aspect("equal"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xticks([0, 0.3, 0.6]); ax.set_yticks([0, 0.3, 0.6])
    ax.tick_params(labelsize=6.5, pad=1.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=7.5, labelpad=1)
    ax.set_ylabel(ylabel, fontsize=7.5, labelpad=0)


def fig_sets(sets):
    fig, axes = plt.subplots(2, 3, figsize=(6.3, 4.15), gridspec_kw={"wspace": 0.22, "hspace": 0.14})
    cfg = {7: {"ph1": np.pi / 2, "ph2": np.pi / 2 + np.pi / 10, "alpha": np.deg2rad(50)},
           8: {"ph1": np.pi / 2, "ph2": np.pi / 2 + np.pi / 5, "alpha": np.pi / 4}}
    for row, n in enumerate((7, 8)):
        S = sets[n]; X, lab, U, W, m1 = S["X"], S["lab"], S["U"], S["W"], S["m1"]
        x1, x2 = X[:, :m1], X[:, m1:]
        _, _, (rA, rB, a, b, c_, d, e, f) = params(n)
        # (a) polygon view: u_i -> regular p-gon, v_j -> regular q-gon; dotted: ideal pole polygons of V_n
        c = cfg[n]
        M1, M2 = polygon_map(U, c["ph1"]), polygon_map(W, c["ph2"])
        Pa = np.cos(c["alpha"]) * x1 @ M1.T + np.sin(c["alpha"]) * x2 @ M2.T
        ax = axes[row, 0]
        for V, r, wt, orb in ((U, rA, np.cos(c["alpha"]), "P"), (W, rB, np.sin(c["alpha"]), "Q")):
            M = M1 if orb == "P" else M2
            poly = wt * r * V @ M.T
            poly = np.r_[poly, poly[:1]]
            ax.plot(poly[:, 0], poly[:, 1], ls=(0, (1, 1.5)), lw=0.6, color=STYLE[orb][0], alpha=0.7, zorder=1)
        scatter(ax, Pa, lab, big=1.9)
        clean(ax, 0.62)
        # (b) cell view: the plane spanned by (u_0, 0) and (0, v_0)
        Pb = np.c_[x1 @ U[0], x2 @ W[0]]
        ax = axes[row, 1]
        ax.axhline(0, color="0.85", lw=0.4, zorder=0); ax.axvline(0, color="0.85", lw=0.4, zorder=0)
        th = np.linspace(-0.25, np.pi / 2 + 0.25, 100)
        ax.plot(0.62 * np.cos(th), 0.62 * np.sin(th), ls=(0, (2, 2)), lw=0.4, color="0.6", zorder=0)
        cell = [("P", 0, None), ("Q", None, 0), ("A", 0, 0), ("B", 0, 0), ("C", 0, 0)]
        pos = {}
        for k, L in enumerate(lab):
            if L in cell and (L not in pos or np.linalg.norm(Pb[k]) > np.linalg.norm(pos[L])):
                pos[L] = Pb[k]
        if n == 8:
            for s, t in ((("P", 0, None), ("B", 0, 0)), (("Q", None, 0), ("A", 0, 0))):
                ax.plot(*zip(pos[s], pos[t]), color="0.2", lw=0.7, zorder=2)
        scatter(ax, Pb, lab, emph=cell)
        names = {("P", 0, None): ("$P_0$", (5, 2)), ("Q", None, 0): ("$Q_0$", (5, 1)),
                 ("A", 0, 0): ("$A_{00}$", (5, -6)), ("B", 0, 0): ("$B_{00}$", (-16, -9)),
                 ("C", 0, 0): ("$C_{00}$", (5, 2))}
        for L, (txt, off) in names.items():
            ax.annotate(txt, pos[L], xytext=off, textcoords="offset points", fontsize=7, zorder=6)
        axes_style(ax, -0.3, 0.78, r"$\langle x',u_0\rangle$" if row == 1 else None, r"$\langle x'',v_0\rangle$")
        # (c) block norms
        ax = axes[row, 2]
        Pc = np.c_[np.linalg.norm(x1, axis=1), np.linalg.norm(x2, axis=1)]
        th = np.linspace(0, np.pi / 2, 100)
        ax.plot(0.62 * np.cos(th), 0.62 * np.sin(th), ls=(0, (2, 2)), lw=0.4, color="0.6", zorder=0)
        scatter(ax, Pc, lab)
        axes_style(ax, -0.03, 0.72, r"$\Vert x'\Vert$" if row == 1 else None, r"$\Vert x''\Vert$")
        if n == 8:
            for (x, y), txt, off in (((rA, 0), r"$(r_A,0)$", (-30, 3)), ((a, b), r"$(a,b)$", (-26, 1)),
                                     ((e, f), r"$(e,f)$", (-27, -3)), ((c_, d), r"$(c,d)$", (2, -10)),
                                     ((0, rB), r"$(0,r_B)$", (4, 2))):
                ax.annotate(txt, (x, y), xytext=off, textcoords="offset points", fontsize=6.5, color="0.25")
        axes[row, 0].text(-0.62, 0.62, f"$X_{n}$", fontsize=10, va="top")
    for j, t in enumerate(("(a) polygon view", "(b) cell view", "(c) block norms")):
        axes[0, j].set_title(t, fontsize=8, pad=3)
    handles = [plt.Line2D([], [], ls="", marker=STYLE[o][1], color=STYLE[o][0], markeredgecolor="white",
                          markeredgewidth=0.3, markersize=5, label=f"${o}$") for o in "PAQBC"]
    fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False, fontsize=8, handletextpad=0.2,
               columnspacing=1.2, bbox_to_anchor=(0.5, -0.02))
    fig.savefig(OUT / "sets.pdf", bbox_inches="tight", pad_inches=0.02, metadata={"CreationDate": None})
    plt.close(fig)


def fig_disk():
    """the disk of diameter D cut by the rays through the vertices of an inscribed equilateral triangle"""
    fig, ax = plt.subplots(figsize=(1.75, 1.75))
    cols = ["#56B4E9", "#E69F00", "#009E73"]
    for k in range(3):
        ax.add_patch(Wedge((0, 0), 0.5, 90 + 120 * k, 210 + 120 * k, fc=cols[k], alpha=0.35, ec="white", lw=1.0))
    ax.add_patch(Circle((0, 0), 0.5, fill=False, lw=0.6, color="0.3"))
    ang = np.deg2rad([90, 210, 330])
    tri = 0.5 * np.c_[np.cos(ang), np.sin(ang)]
    for k in range(3):
        ax.plot([0, tri[k, 0]], [0, tri[k, 1]], color="white", lw=1.0, zorder=2)
    ax.plot([tri[1, 0], tri[2, 0], tri[0, 0]], [tri[1, 1], tri[2, 1], tri[0, 1]], color="0.35", lw=0.6,
            ls=(0, (3, 2)), zorder=3)
    ax.plot([tri[0, 0], tri[1, 0]], [tri[0, 1], tri[1, 1]], color="black", lw=0.9, zorder=4)
    ax.plot(tri[:, 0], tri[:, 1], ls="", marker="o", ms=2.5, color="black", zorder=5)
    mid = (tri[0] + tri[1]) / 2
    ax.text(mid[0] + 0.09, mid[1] - 0.052, r"$\sqrt{3}\,D/2$", fontsize=8, rotation=60, ha="center", va="center",
            rotation_mode="anchor")
    ax.annotate("", xy=(0.5, -0.62), xytext=(-0.5, -0.62),
                arrowprops={"arrowstyle": "<->", "lw": 0.5, "color": "0.3", "shrinkA": 0, "shrinkB": 0})
    ax.text(0, -0.66, "$D$", ha="center", va="top", fontsize=8)
    ax.set_aspect("equal"); ax.set_xlim(-0.56, 0.56); ax.set_ylim(-0.74, 0.54); ax.axis("off")
    fig.savefig(OUT / "disk.pdf", bbox_inches="tight", pad_inches=0.01, metadata={"CreationDate": None})
    plt.close(fig)


def fig_hgraph(p=3, q=4):
    """the vertex layout of H(p, q): bands mark the row blocks (P_i and the A_ij) and column blocks (Q_j and the B_ij)"""
    fig, ax = plt.subplots(figsize=(2.55, 2.1))
    dx, dy, top = 1.0, 0.85, 0.75
    xa, xb, yb, yc = -0.24, 0.22, 0.26, -0.27
    w = 0.13
    for i in range(p):
        y = -i * dy
        ax.add_patch(Rectangle((-0.68, y - w / 2), (q - 1) * dx + xa + 0.86, w, fc=STYLE["P"][0], alpha=0.16, lw=0))
        ax.scatter([-0.55], [y], s=24, c=STYLE["P"][0], marker="o", zorder=3)
        ax.text(-0.72, y, f"$P_{i}$", ha="right", va="center", fontsize=7.5)
    for j in range(q):
        x = j * dx
        ax.add_patch(Rectangle((x + xb - w / 2, -(p - 1) * dy + yb - 0.12), w, top + (p - 1) * dy - yb + 0.25,
                               fc=STYLE["Q"][0], alpha=0.16, lw=0))
        ax.scatter([x + xb], [top], s=20, c=STYLE["Q"][0], marker="s", zorder=3)
        ax.text(x + xb, top + 0.15, f"$Q_{j}$", ha="center", va="bottom", fontsize=7.5)
        for i in range(p):
            y = -i * dy
            for orb, (px, py) in (("A", (x + xa, y)), ("B", (x + xb, y + yb)), ("C", (x, y + yc))):
                col, mk = STYLE[orb]
                ax.scatter([px], [py], s=13, c=col, marker=mk, zorder=3, linewidths=0.3, edgecolors="white")
    for txt, (px, py), off in (("$A_{00}$", (xa, 0), (-9, 4)), ("$B_{00}$", (xb, yb), (4, -1)),
                               ("$C_{00}$", (0, yc), (4, -6))):
        ax.annotate(txt, (px, py), xytext=off, textcoords="offset points", fontsize=6.5, color="0.15", zorder=6)
    ax.set_xlim(-1.05, (q - 1) * dx + 0.45); ax.set_ylim(-(p - 1) * dy - 0.45, top + 0.42)
    ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(OUT / "hgraph.pdf", bbox_inches="tight", pad_inches=0.01, metadata={"CreationDate": None})
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    sets = load()
    fig_disk()
    fig_sets(sets)
    fig_hgraph()
    print("wrote", sorted(p.name for p in OUT.glob("*.pdf")))


if __name__ == "__main__":
    main()
