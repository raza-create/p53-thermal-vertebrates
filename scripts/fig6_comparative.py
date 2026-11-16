#!/usr/bin/env python3
"""Figure 3 (comparative): per-residue heat sensitivity + conservation tracks (A), identity by position set (B), surface vs buried (C)."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import fisher_exact
O = Path("results/corrected/figs_animals")
INK, INK2, MUTED, GRID, AXIS = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5, "axes.edgecolor": AXIS, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False})
REG = [("N-terminal strands (94-163)", "N-terminal strands", "#2a78d6"), ("L2 (167-194)", "L2 loop", "#eb6834"),
       ("beta-core (195-236)", "\u03b2-core", "#1baf7a"), ("L3 (237-250)", "L3 loop", "#eda100"),
       ("C-terminal S10/H2 (251-292)", "C-terminal S10/H2", "#e87ba4")]
COL = {r: c for r, _, c in REG}; COL["other"] = AXIS
m = pd.read_csv("results/corrected/per_residue_heat_conservation_rsa.csv")
ident = pd.read_csv(O / "table_species_identity.csv", index_col=0)
fig = plt.figure(figsize=(7.1, 6.6))
outer = fig.add_gridspec(2, 1, height_ratios=[1.25, 1.0], hspace=0.42, left=0.16, right=0.985, top=0.90, bottom=0.08)
top = outer[0].subgridspec(2, 1, height_ratios=[3.2, 0.75], hspace=0.06)
a = fig.add_subplot(top[0]); t = fig.add_subplot(top[1], sharex=a)
a.bar(m.resid, m.dQ350, width=0.8, color=[COL[r] for r in m.region], linewidth=0); a.axhline(0, color=AXIS, lw=0.8)
a.set(ylabel="Q lost, 248\u2192350 K", xlim=(89, 291), ylim=(-0.42, 1.08), yticks=[-0.25, 0, 0.25, 0.5, 0.75, 1.0])
a.grid(axis="y", color=GRID, lw=0.6); a.set_axisbelow(True)
a.tick_params(axis="x", bottom=False, labelbottom=False); a.spines["bottom"].set_visible(False)
a.set_title("A   Per-residue heat sensitivity and conservation (33 vertebrates)", pad=20)
a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in REG], labels=[l for _, l, _ in REG],
         ncol=5, loc="lower left", bbox_to_anchor=(0, 1.0), fontsize=7, handlelength=1.0, columnspacing=1.0, borderaxespad=0.2)
q = m.set_index("resid").dQ350
for r, txt, dx, dy in [(178, "H178", -14, 6), (186, "D186", 0, 4), (223, "P223", -16, 4)]:
    if r in q.index: a.annotate(txt, (r, q[r]), xytext=(dx, dy), textcoords="offset points", ha="center", va="bottom", fontsize=6.5, color=INK2)
if all(r in q.index for r in (227, 228, 229)):
    ytop = q[[227, 228, 229]].max(); a.plot([226.6, 229.4], [ytop + 0.04] * 2, color=INK2, lw=0.7)
    a.text(228, ytop + 0.06, "S227\u2013C229", ha="center", va="bottom", fontsize=6.5, color=INK2)
t.scatter(m.resid[m.invariant], np.ones(m.invariant.sum()), marker="|", s=34, color=INK, lw=0.9)
t.scatter(m.resid[m.divergent], np.zeros(m.divergent.sum()), marker="|", s=34, color=MUTED, lw=0.9)
t.set(yticks=[0, 1], yticklabels=[f"Divergent ({m.divergent.sum()})", f"Invariant ({m.invariant.sum()})"], ylim=(-0.7, 1.7), xlabel="Residue (UniProt P04637)")
t.spines["left"].set_visible(False); t.tick_params(axis="y", length=0)
bottom = outer[1].subgridspec(1, 2, wspace=0.55, width_ratios=[1.05, 1.0])
b = fig.add_subplot(bottom[0]); yy = np.arange(len(ident))[::-1]
for y_ in yy: b.axhline(y_, color=GRID, lw=0.6, zorder=0)
for lab, mk, fc, ec in [("L2 + L3 loops", "o", INK, INK), ("Whole DBD", "o", "none", INK), ("Divergent positions", "^", MUTED, MUTED)]:
    b.scatter(ident[lab], yy, marker=mk, s=26, facecolors=fc, edgecolors=ec, lw=1.1, label=lab, zorder=3)
b.set(yticks=yy, yticklabels=ident.index, xlabel="Identity to human (%)", xlim=(0, 104), ylim=(-0.6, len(ident) - 0.4))
b.set_title("B   Identity by position set"); b.legend(fontsize=6.5, loc="upper left", handletextpad=0.3, borderaxespad=0.2)
c = fig.add_subplot(bottom[1]); s = m.rsa >= 0.20
cnt = {"Surface\n(RSA \u2265 0.2)": (m[s].divergent.sum(), m[s].invariant.sum(), s.sum()), "Buried": (m[~s].divergent.sum(), m[~s].invariant.sum(), (~s).sum())}
x = np.arange(2); w = 0.36
for j, (nm, col, idx) in enumerate([("Divergent", MUTED, 0), ("Invariant", INK2, 1)]):
    for xi, k_ in zip(x, cnt):
        v = 100 * cnt[k_][idx] / cnt[k_][2]; xx = xi + (j - 0.5) * (w + 0.03)
        c.bar(xx, v, width=w, color=col, label=nm if xi == 0 else None)
        c.text(xx, v + 1.2, f"{cnt[k_][idx]}/{cnt[k_][2]}", ha="center", va="bottom", fontsize=6.3, color=INK2)
o1, p1 = fisher_exact([[cnt[k_][0], cnt[k_][2] - cnt[k_][0]] for k_ in cnt]); o2, p2 = fisher_exact([[cnt[k_][1], cnt[k_][2] - cnt[k_][1]] for k_ in cnt])
c.set(xticks=x, xticklabels=list(cnt), ylabel="Positions (%)", ylim=(0, 85), yticks=[0, 20, 40, 60])
c.grid(axis="y", color=GRID, lw=0.6); c.set_axisbelow(True); c.set_title("C   Surface vs buried positions")
c.legend(fontsize=6.5, loc="upper left", borderaxespad=0.2, handlelength=1.2)
c.text(0.99, 0.99, f"Divergent: OR {o1:.2f}, p = {p1:.3f}\nInvariant: OR {o2:.2f}, p = {p2:.3f}", transform=c.transAxes, ha="right", va="top", fontsize=6.3, color=INK2)
fig.savefig(O / "fig6_comparative.png", dpi=300, facecolor="white"); fig.savefig(O / "fig6_comparative.pdf", facecolor="white")
print(f"Wrote {O}/fig6_comparative.(png|pdf); Fisher divergent OR={o1:.2f} p={p1:.4f}, invariant OR={o2:.2f} p={p2:.4f}")
