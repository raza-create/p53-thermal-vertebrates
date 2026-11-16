#!/usr/bin/env python3
"""Main corrected figures for the Animals manuscript. Out: results/corrected/figs_animals/"""
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import fisher_exact
O = Path("results/corrected/figs_animals"); O.mkdir(parents=True, exist_ok=True)
INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.edgecolor": AXIS, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False})
REG = [("N-terminal strands (94-163)", "N-terminal strands", "#2a78d6"), ("L2 (167-194)", "L2 loop", "#eb6834"),
       ("beta-core (195-236)", "β-core", "#1baf7a"), ("L3 (237-250)", "L3 loop", "#eda100"),
       ("C-terminal S10/H2 (251-292)", "C-terminal S10/H2", "#e87ba4")]
COL = {r: c for r, _, c in REG}; COL["other"] = AXIS; TS = [248, 273, 310, 350, 400, 450]
d = pd.read_csv("data/features_corrected/descriptors_corrected.csv")
m = pd.read_csv("results/corrected/per_residue_heat_conservation_rsa.csv")
# ---------- Figure 3: thermal response ----------
fig, ax = plt.subplots(1, 2, figsize=(7.1, 2.9))
g = d.groupby("temperature").Q_bhe45.agg(["mean", "std"]); a = ax[0]
a.grid(axis="y", color=GRID, lw=0.6); a.set_axisbelow(True)
a.axvline(310, color=MUTED, lw=1, ls="--"); a.text(314, 0.47, "human body\ntemperature", color=INK2, fontsize=7)
a.errorbar(g.index, g["mean"], yerr=g["std"], color="#2a78d6", lw=2, marker="o", ms=5, capsize=2, mec=SURF, mew=1)
for T in (248, 310, 350, 450): a.annotate(f"{g.loc[T, 'mean']:.3f}", (T, g.loc[T, "mean"]), xytext=(6, 4), textcoords="offset points", fontsize=7, color=INK2)
a.set(xlabel="Simulation temperature (K)", ylabel="Native-contact fraction Q", ylim=(0.45, 1.0), xticks=TS); a.set_title("A   Global native contacts")
b = ax[1]; reg = m.groupby("region")[[str(T) for T in TS]].mean(); lost = 100 * (1 - reg.div(reg["248"], axis=0))
b.grid(axis="y", color=GRID, lw=0.6); b.set_axisbelow(True); b.axvline(310, color=MUTED, lw=1, ls="--")
ends = []
for r, lab, c in REG:
    b.plot(TS, lost.loc[r].values, color=c, lw=2, marker="o", ms=4, mec=SURF, mew=0.8, label=lab); ends.append([lost.loc[r, "450"], lab])
ends.sort(); y = [e[0] for e in ends]
for i in range(1, len(y)): y[i] = max(y[i], y[i - 1] + 6.5)
for (v, lab), yy in zip(ends, y): b.text(458, yy, f"{lab} ({v:.0f}%)", fontsize=6.5, color=INK2, va="center")
b.set(xlabel="Simulation temperature (K)", ylabel="Native contacts lost vs 248 K (%)", xlim=(240, 545), xticks=TS); b.set_title("B   By structural region")
b.legend(loc="upper left", fontsize=6.5, handlelength=1.5)
for x_ in ax: x_.set_xticks(TS); x_.set_xticklabels([str(T) for T in TS], rotation=45, ha="right")
fig.tight_layout(); fig.savefig(O / "fig3_thermal_response.png", dpi=300, facecolor="white"); fig.savefig(O / "fig3_thermal_response.pdf"); plt.close(fig)
# ---------- per-species identity ----------
def read_fa(p):
    s, k = {}, None
    for l in open(p):
        l = l.strip()
        if not l: continue
        if l.startswith(">"): k = l[1:].split("|")[1]; s[k] = []
        else: s[k].append(l)
    return {a_: "".join(b_) for a_, b_ in s.items()}
aln = read_fa("species/vert_p53_aln.fasta"); H = aln["P04637"]; p2c, k = {}, 0
for c, a_ in enumerate(H):
    if a_ != "-": k += 1; p2c[k] = c
SETS = {"L2 + L3 loops": m[m.region.isin(["L2 (167-194)", "L3 (237-250)"])].resid, "Whole DBD": m.resid, "Divergent positions": m[m.divergent].resid}
SP = [("P04637", "Human"), ("P02340", "Mouse"), ("Q29537", "Dog"), ("P41685", "Cat"), ("P10360", "Chicken"), ("P07193", "Xenopus (ecto.)"), ("P79734", "Zebrafish (ecto.)")]
ident = pd.DataFrame([{lab: 100 * np.mean([aln[acc][p2c[p]] == H[p2c[p]] for p in pos if aln[acc][p2c[p]] != "-"]) for lab, pos in SETS.items()} for acc, _ in SP], index=[n for _, n in SP])
# ---------- Figure 6: comparative ----------
fig = plt.figure(figsize=(7.1, 5.6)); gs = fig.add_gridspec(3, 2, height_ratios=[2.3, 0.5, 2.6], hspace=0.35, wspace=0.45)
a = fig.add_subplot(gs[0, :]); t = fig.add_subplot(gs[1, :], sharex=a)
a.bar(m.resid, m.dQ350, width=0.8, color=[COL[r] for r in m.region], linewidth=0)
for _, r in m.nlargest(6, "dQ350").iterrows(): a.annotate(f"{r.human}{r.resid}", (r.resid, r.dQ350), xytext=(0, 2), textcoords="offset points", ha="center", fontsize=6.5, color=INK2)
a.set(ylabel="Q lost, 248→350 K", xlim=(89, 291)); a.grid(axis="y", color=GRID, lw=0.6); a.set_axisbelow(True); a.tick_params(labelbottom=False)
a.set_title("A   Per-residue heat sensitivity (zinc-free DBD) and conservation across 33 vertebrates")
a.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c) for _, _, c in REG], labels=[l for _, l, _ in REG], ncol=5, fontsize=6.5, loc="upper left", bbox_to_anchor=(0, 1.0), handlelength=1)
t.scatter(m.resid[m.invariant], np.ones(m.invariant.sum()), marker="|", s=40, color=INK, lw=1)
t.scatter(m.resid[m.divergent], np.zeros(m.divergent.sum()), marker="|", s=40, color=MUTED, lw=1)
t.set(yticks=[0, 1], yticklabels=[f"Endo/ecto divergent (n={m.divergent.sum()})", f"Invariant (n={m.invariant.sum()})"], ylim=(-0.6, 1.6), xlabel="Residue (UniProt P04637)")
t.spines["left"].set_visible(False); t.tick_params(axis="y", length=0)
b = fig.add_subplot(gs[2, 0]); yy = np.arange(len(ident))[::-1]
for lab, mk, fc, ec in [("L2 + L3 loops", "o", INK, INK), ("Whole DBD", "o", "none", INK), ("Divergent positions", "^", MUTED, MUTED)]:
    b.scatter(ident[lab], yy, marker=mk, s=30, facecolors=fc, edgecolors=ec, lw=1.2, label=lab, zorder=3)
for y_ in yy: b.axhline(y_, color=GRID, lw=0.6, zorder=0)
b.set(yticks=yy, yticklabels=ident.index, xlabel="Identity to human (%)", xlim=(0, 105)); b.set_title("B   Identity by position set")
b.legend(fontsize=6.5, loc="lower left", handletextpad=0.3)
c = fig.add_subplot(gs[2, 1]); mm = m.copy(); mm["surface"] = mm.rsa >= 0.20
cnt = {}
for lab, sel in [("Surface\n(RSA ≥ 0.2)", mm.surface), ("Buried", ~mm.surface)]:
    cnt[lab] = (mm[sel].divergent.sum(), mm[sel].invariant.sum(), sel.sum())
x = np.arange(2); w = 0.36
for j, (nm, col, idx) in enumerate([("Divergent", MUTED, 0), ("Invariant", INK2, 1)]):
    v = [100 * cnt[k_][idx] / cnt[k_][2] for k_ in cnt]
    c.bar(x + (j - 0.5) * (w + 0.03), v, width=w, color=col, label=nm)
    for xi, k_ in zip(x, cnt): c.text(xi + (j - 0.5) * (w + 0.03), 100 * cnt[k_][idx] / cnt[k_][2] + 1.5, f"{cnt[k_][idx]}/{cnt[k_][2]}", ha="center", fontsize=6.5, color=INK2)
o1, p1 = fisher_exact([[cnt[k_][0], cnt[k_][2] - cnt[k_][0]] for k_ in cnt]); o2, p2 = fisher_exact([[cnt[k_][1], cnt[k_][2] - cnt[k_][1]] for k_ in cnt])
c.set(xticks=x, xticklabels=list(cnt), ylabel="Positions (%)", ylim=(0, 70)); c.grid(axis="y", color=GRID, lw=0.6); c.set_axisbelow(True)
c.legend(fontsize=6.5, loc="upper right"); c.set_title("C   Surface vs buried positions")
c.text(0.02, 0.97, f"Divergent: OR {o1:.2f}, p = {p1:.3f}\nInvariant: OR {o2:.2f}, p = {p2:.3f}", transform=c.transAxes, fontsize=6.5, color=INK2, va="top")
fig.savefig(O / "fig6_comparative.png", dpi=300, facecolor="white", bbox_inches="tight"); fig.savefig(O / "fig6_comparative.pdf", bbox_inches="tight"); plt.close(fig)
top = m.nlargest(20, "dQ350"); ot, pt = fisher_exact([[ (top.rsa >= 0.2).sum(), 20 - (top.rsa >= 0.2).sum()], [(m.rsa >= 0.2).sum() - (top.rsa >= 0.2).sum(), len(m) - 20 - ((m.rsa >= 0.2).sum() - (top.rsa >= 0.2).sum())]])
print(f"top-20 heat-sensitive surface-exposed: {(top.rsa >= 0.2).sum()}/20 vs rest; Fisher OR={ot:.2f} p={pt:.3g}")
print("Per-species identity (%):"); print(ident.round(1).to_string())
ident.round(1).to_csv(O / "table_species_identity.csv"); print(f"Wrote {O}/fig3_thermal_response.(png|pdf), fig6_comparative.(png|pdf), table_species_identity.csv")
