#!/usr/bin/env python3
"""Supplementary Figures S1-S2 and Tables S1-S4 for the Animals manuscript. Out: results/corrected/supplementary/"""
import csv, numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, matplotlib.image as mpimg
from pathlib import Path
O = Path("results/corrected/supplementary"); O.mkdir(parents=True, exist_ok=True); C = Path("results/corrected")
INK2, GRID, AXIS = "#52514e", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": AXIS,
                     "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 9, "axes.titleweight": "bold", "axes.titlelocation": "left", "legend.frameon": False})
TS = [248, 273, 310, 350, 400, 450]; RAMP = ["#86b6ef", "#4f93e5", "#2a78d6", "#1f63b8", "#174f94", "#0d366b"]
d = pd.read_csv("data/features_corrected/descriptors_corrected.csv")
# Figure S1
P = [("Q_bhe45", "Native-contact fraction Q"), ("rmsd", "Cα RMSD (Å)"), ("rg", "Radius of gyration (Å)"),
     ("end_to_end", "End-to-end distance (Å)"), ("hbonds_backbone", "Backbone H-bonds"), ("asphericity", "Asphericity (Å²)")]
fig, ax = plt.subplots(3, 2, figsize=(7.1, 7.4), sharex=True)
for i, (c, lab) in enumerate(P):
    a = ax.flat[i]; a.grid(axis="y", color=GRID, lw=0.6); a.set_axisbelow(True)
    for T, col in zip(TS, RAMP):
        s = d[d.temperature == T].sort_values("time_ns"); a.plot(s.time_ns, s[c].rolling(10, center=True, min_periods=1).mean(), color=col, lw=1.4, label=f"{T} K")
    a.set_ylabel(lab); a.set_title(f"{'ABCDEF'[i]}", loc="left")
for a in ax[-1]: a.set_xlabel("Time (ns)")
h, l = ax.flat[0].get_legend_handles_labels(); fig.legend(h, l, ncol=6, loc="upper center", bbox_to_anchor=(0.5, 1.0), fontsize=7.5)
fig.tight_layout(rect=(0, 0, 1, 0.965)); fig.savefig(O / "FigS1_timeseries.png", dpi=300, facecolor="white"); fig.savefig(O / "FigS1_timeseries.pdf"); plt.close(fig)
# Figure S2
fig = plt.figure(figsize=(7.1, 5.0)); gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.75])
a1 = fig.add_subplot(gs[0]); a1.imshow(mpimg.imread(C / "fig4_corrected.png")); a1.axis("off")
a2 = fig.add_subplot(gs[1]); a2.imshow(mpimg.imread(C / "fig3c_shap_corrected.png")); a2.axis("off"); a2.set_title("C   SHAP importance (share of mean |SHAP|)", fontsize=9)
fig.tight_layout(); fig.savefig(O / "FigS2_rf_shap.png", dpi=300, facecolor="white"); fig.savefig(O / "FigS2_rf_shap.pdf"); plt.close(fig)
# Table S1
cols = {"Q_bhe45": "Q (BHE, 4.5 A)", "Q_8A": "Q (8 A sensitivity)", "rmsd": "CA RMSD (A)", "rg": "Rg (A)", "end_to_end": "End-to-end (A)",
        "hbonds_backbone": "Backbone H-bonds", "hbonds": "Total H-bonds", "asphericity": "Asphericity (A^2)"}
g = d.groupby("temperature")[list(cols)].agg(["mean", "std"])
S1 = pd.DataFrame({cols[c]: [f"{g.loc[T, (c, 'mean')]:.3f} ± {g.loc[T, (c, 'std')]:.3f}" for T in TS] for c in cols}, index=pd.Index(TS, name="Temperature (K)"))
S1.loc["Pearson r with T (frames)"] = [f"{d[c].corr(d.temperature):.3f}" for c in cols]
# Table S2
S2a, S2b, S2c = pd.read_csv(C / "table5_loto.csv"), pd.read_csv(C / "table4_shap.csv"), pd.read_csv(C / "tableS2_shap_per_fold.csv")
# Table S3
m = pd.read_csv(C / "per_residue_heat_conservation_rsa.csv")
S3 = m.rename(columns={"resid": "Residue (UniProt P04637)", "human": "Human aa", "region": "Region", **{str(T): f"Q_res {T} K" for T in TS},
                       "dQ310": "Q lost 248-310 K", "dQ350": "Q lost 248-350 K", "dQ400": "Q lost 248-400 K", "entropy": "Shannon entropy (bits)",
                       "invariant": "Invariant (33 spp.)", "endo_maj": "Endotherm majority aa (n/N)", "ecto_maj": "Ectotherm majority aa (n/N)",
                       "divergent": "Endo/ecto divergent", "rsa": "RSA (2XWR chain A)", "surface": "Surface (RSA>=0.2)"}).drop(columns=["loops"], errors="ignore").round(4)
# Table S4
def read_fa(p):
    s, k = {}, None
    for l in open(p):
        l = l.strip()
        if not l: continue
        if l.startswith(">"): k = l[1:].split("|")[1]; s[k] = []
        else: s[k].append(l)
    return {a: "".join(b) for a, b in s.items()}
aln = read_fa("species/vert_p53_aln.fasta"); H = aln["P04637"]; p2c, k = {}, 0
for c_, a_ in enumerate(H):
    if a_ != "-": k += 1; p2c[k] = c_
rows = []
with open("species/uniprot_vert_p53.tsv") as f:
    r = csv.DictReader(f, delimiter="\t"); lc = [c for c in r.fieldnames if "ineage" in c][0]
    for row in r:
        acc = row["Entry"]; s = aln[acc]; cols_ = [p2c[i] for i in range(94, 293)]; ng = [c for c in cols_ if s[c] != "-"]
        cls = next((c for c in ["Mammalia", "Aves", "Amphibia", "Actinopteri"] if c in row[lc]), "Other")
        rows.append({"UniProt accession": acc, "Species": row["Organism"], "Taxon ID": row["Organism (ID)"], "Class": cls,
                     "Thermal group": "Endotherm" if cls in ("Mammalia", "Aves") else "Ectotherm", "Length (aa)": int(row["Length"]),
                     "DBD coverage (94-292)": round(len(ng) / len(cols_), 3), "DBD identity to human": round(sum(s[c] == H[c] for c in ng) / len(ng), 3),
                     "Included": "No (DBD coverage < 0.9)" if acc == "Q29480" else "Yes"})
S4 = pd.DataFrame(rows).sort_values(["Thermal group", "Class", "DBD identity to human"], ascending=[True, False, False])
tabs = {"S1_descriptors": S1.reset_index(), "S2a_LOTO": S2a, "S2b_SHAP_global": S2b, "S2c_SHAP_per_fold": S2c, "S3_per_residue": S3, "S4_orthologues": S4}
try:
    with pd.ExcelWriter(O / "Supplementary_Tables_S1-S4.xlsx") as w:
        for n, t in tabs.items(): t.to_excel(w, sheet_name=n, index=False)
    print(f"Wrote {O}/Supplementary_Tables_S1-S4.xlsx")
except ImportError as e:
    print(f"openpyxl missing ({e}); writing CSVs instead")
for n, t in tabs.items(): t.to_csv(O / f"Table_{n}.csv", index=False)
print("S1:"); print(S1.to_string()); print(f"\nS4: {len(S4)} entries, included {(S4.Included == 'Yes').sum()}; classes: {S4[S4.Included=='Yes'].Class.value_counts().to_dict()}")
print(f"Wrote {O}/FigS1_timeseries.(png|pdf), FigS2_rf_shap.(png|pdf), Table_*.csv")
