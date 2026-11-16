#!/usr/bin/env python3
"""Per-residue heat sensitivity (corrected BHE Q) vs conservation across 33 vertebrate p53 orthologues."""
import csv, numpy as np, pandas as pd
from collections import Counter
from scipy.stats import spearmanr, mannwhitneyu
OUT = "results/corrected"; TS = [248, 273, 310, 350, 400, 450]
pr = pd.read_csv("data/features_corrected/per_residue_Q_bhe45.csv"); pr.columns = ["resid"] + [int(c) for c in pr.columns[1:]]
n0 = len(pr); pr = pr[pr[248] > 0].copy(); print(f"residues with native contacts: {len(pr)}/{n0}")
for T in (310, 350, 400): pr[f"dQ{T}"] = pr[248] - pr[T]
REG = [("N-terminal strands (94-163)", 94, 163), ("L2 (167-194)", 167, 194), ("beta-core (195-236)", 195, 236),
       ("L3 (237-250)", 237, 250), ("C-terminal S10/H2 (251-292)", 251, 292)]
pr["region"] = pr.resid.map(lambda n: next((r for r, a, b in REG if a <= n <= b), "other"))
reg = pr.groupby("region")[TS].mean(); lost = (100 * (1 - reg.div(reg[248], axis=0))).round(1)
print("\nPer-region mean residue Q (BHE):"); print(reg.round(3).to_string()); print("\n% lost vs 248 K:"); print(lost.to_string())
def read_fa(p):
    s, n = {}, None
    for l in open(p):
        l = l.strip()
        if not l: continue
        if l.startswith(">"): n = l[1:].split("|")[1]; s[n] = []
        else: s[n].append(l)
    return {k: "".join(v) for k, v in s.items()}
aln = read_fa("species/vert_p53_aln.fasta"); aln.pop("Q29480", None)   # donkey: partial DBD
cls = {}
with open("species/uniprot_vert_p53.tsv") as f:
    r = csv.DictReader(f, delimiter="\t"); lc = [c for c in r.fieldnames if "ineage" in c][0]
    for row in r: cls[row["Entry"]] = next((c for c in ["Mammalia", "Aves", "Amphibia", "Actinopteri"] if c in row[lc]), "Other")
endo = [a for a in aln if cls[a] in ("Mammalia", "Aves")]; ecto = [a for a in aln if cls[a] in ("Actinopteri", "Amphibia")]
print(f"\northologues: {len(aln)} (endotherms {len(endo)}, ectotherms {len(ecto)})")
H = aln["P04637"]; p2c, k = {}, 0
for c, a in enumerate(H):
    if a != "-": k += 1; p2c[k] = c
rows = []
for pos in pr.resid:
    c = p2c[pos]; allaa = [aln[a][c] for a in aln if aln[a][c] != "-"]
    e1 = [aln[a][c] for a in endo if aln[a][c] != "-"]; e2 = [aln[a][c] for a in ecto if aln[a][c] != "-"]
    fr = np.array(list(Counter(allaa).values())) / len(allaa)
    m1, m2 = Counter(e1).most_common(1)[0], Counter(e2).most_common(1)[0]
    rows.append(dict(resid=pos, human=H[c], entropy=float(-(fr * np.log2(fr)).sum()), invariant=len(set(allaa)) == 1,
                     endo_maj=f"{m1[0]}({m1[1]}/{len(e1)})", ecto_maj=f"{m2[0]}({m2[1]}/{len(e2)})", divergent=m1[0] != m2[0]))
m = pr.merge(pd.DataFrame(rows), on="resid"); m.to_csv(f"{OUT}/per_residue_heat_conservation.csv", index=False)
print(f"invariant positions: {m.invariant.sum()}/{len(m)}; endotherm-vs-ectotherm divergent: {m.divergent.sum()}/{len(m)}")
print("\nby region: n, invariant, divergent, mean dQ350")
print(m.groupby("region").agg(n=("resid", "size"), invariant=("invariant", "sum"), divergent=("divergent", "sum"), mean_dQ350=("dQ350", "mean")).round(3).to_string())
for T in (310, 350, 400):
    rho, p = spearmanr(m[f"dQ{T}"], m.entropy)
    u1 = mannwhitneyu(m.loc[m.divergent, f"dQ{T}"], m.loc[~m.divergent, f"dQ{T}"]); u2 = mannwhitneyu(m.loc[m.invariant, f"dQ{T}"], m.loc[~m.invariant, f"dQ{T}"])
    print(f"\ndQ{T}: Spearman(heat loss, entropy) rho={rho:.3f} p={p:.3g} | divergent median {m.loc[m.divergent, f'dQ{T}'].median():.3f} vs other {m.loc[~m.divergent, f'dQ{T}'].median():.3f}, MWU p={u1.pvalue:.3g} | invariant median {m.loc[m.invariant, f'dQ{T}'].median():.3f} vs variable {m.loc[~m.invariant, f'dQ{T}'].median():.3f}, MWU p={u2.pvalue:.3g}")
print("\nTop 20 heat-sensitive residues (loss 248->350 K):")
print(m.sort_values("dQ350", ascending=False).head(20)[["resid", "human", "region", 248, 310, 350, "dQ350", "entropy", "invariant", "endo_maj", "ecto_maj", "divergent"]].round(3).to_string(index=False))
d = pd.read_csv("data/features_corrected/descriptors_corrected.csv")
print(f"\nQ_bhe45 vs RMSD across frames: r = {d.Q_bhe45.corr(d.rmsd):.3f}; within 310 K: r = {d[d.temperature==310].Q_bhe45.corr(d[d.temperature==310].rmsd):.3f}")
print(f"\nWrote {OUT}/per_residue_heat_conservation.csv")
