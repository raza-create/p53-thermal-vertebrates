#!/usr/bin/env python3
"""Correlation, RF + LOTO, SHAP on corrected descriptors (manuscript settings). Outputs: results/corrected/"""
import numpy as np, pandas as pd, shap, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import pearsonr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_absolute_error
OUT = Path("results/corrected"); OUT.mkdir(parents=True, exist_ok=True)
df = pd.read_csv("data/features_corrected/descriptors_corrected.csv")
FEAT = {"Q_bhe45": "Native contacts", "rmsd": "RMSD", "rg": "Radius of gyration",
        "end_to_end": "End-to-end distance", "hbonds_backbone": "Backbone H-bonds", "asphericity": "Asphericity"}
RF = dict(n_estimators=200, max_depth=15, min_samples_split=5, random_state=42, n_jobs=-1)
y = df.temperature.values; TS = sorted(df.temperature.unique())
print(f"frames: {len(df)} ({len(df)//len(TS)} per temperature)")
# Table 2: correlations
t2 = pd.DataFrame([(FEAT[c], pearsonr(df[c], y)[0]) for c in FEAT], columns=["Descriptor", "Pearson r"]).sort_values("Pearson r", key=abs, ascending=False)
t2.to_csv(OUT / "table2_correlations.csv", index=False); print("\nTable 2 (r with temperature):"); print(t2.round(3).to_string(index=False))
# Table 3: Q per temperature
q = df.groupby("temperature").Q_bhe45.agg(["mean", "std"]); q["pct_retained_vs248"] = 100 * q["mean"] / q.loc[248, "mean"]; q["pct_lost"] = 100 - q.pct_retained_vs248
q.round(3).to_csv(OUT / "table3_Q.csv"); print("\nTable 3 (Q, BHE 4.5 A):"); print(q.round(3).to_string())
a, b = df.loc[df.temperature == 248, "Q_bhe45"], df.loc[df.temperature == 450, "Q_bhe45"]
sp = np.sqrt(((len(a)-1)*a.var() + (len(b)-1)*b.var()) / (len(a)+len(b)-2)); print(f"\nCohen's d (Q, 248 vs 450 K) = {(a.mean()-b.mean())/sp:.2f}")
def run(feats, tag, figs):
    X = df[list(feats)]
    m = RandomForestRegressor(**RF).fit(X, y); p = m.predict(X)
    print(f"\n[{tag}] in-sample R2 = {r2_score(y, p):.3f}, MAE = {mean_absolute_error(y, p):.1f} K")
    loto, fold = [], {}
    for T in TS:
        tr = y != T; mm = RandomForestRegressor(**RF).fit(X[tr], y[tr]); pp = mm.predict(X[~tr]); fold[T] = mm
        loto.append(dict(heldout_K=T, MAE_K=np.abs(pp - T).mean(), pred_mean_K=pp.mean(), pred_sd_K=pp.std()))
    loto = pd.DataFrame(loto); print(f"[{tag}] LOTO:"); print(loto.round(1).to_string(index=False)); print(f"[{tag}] mean LOTO MAE = {loto.MAE_K.mean():.1f} K")
    samp = df.groupby("temperature", group_keys=False).apply(lambda g: g.sample(833, random_state=42))
    sv = shap.TreeExplainer(m).shap_values(samp[list(feats)]); imp = np.abs(sv).mean(0)
    rng = np.random.default_rng(42); boots = np.array([np.abs(sv[rng.integers(0, len(sv), len(sv))]).mean(0) for _ in range(1000)])
    shp = pd.DataFrame(dict(Descriptor=[feats[c] for c in feats], mean_abs_SHAP_K=imp, CI_low=np.percentile(boots, 2.5, 0),
                            CI_high=np.percentile(boots, 97.5, 0), share_pct=100 * imp / imp.sum())).sort_values("share_pct", ascending=False)
    print(f"[{tag}] SHAP (n={len(samp)}):"); print(shp.round(2).to_string(index=False))
    per = []
    for T, mm in fold.items():
        s2 = df[df.temperature != T].groupby("temperature", group_keys=False).apply(lambda g: g.sample(100, random_state=42))
        v = np.abs(shap.TreeExplainer(mm).shap_values(s2[list(feats)])).mean(0); per.append(dict(heldout_K=T, **{feats[c]: 100*v[i]/v.sum() for i, c in enumerate(feats)}))
    per = pd.DataFrame(per); print(f"[{tag}] per-fold SHAP share (%):"); print(per.round(1).to_string(index=False))
    if figs:
        loto.to_csv(OUT / "table5_loto.csv", index=False); shp.to_csv(OUT / "table4_shap.csv", index=False); per.to_csv(OUT / "tableS2_shap_per_fold.csv", index=False)
        fig, ax = plt.subplots(1, 2, figsize=(10, 4.4))
        idx = np.random.default_rng(0).choice(len(y), 3000, replace=False)
        ax[0].scatter(y[idx] + np.random.default_rng(1).normal(0, 3, len(idx)), p[idx], s=4, alpha=.35, c="#2b6cb0")
        ax[0].set_title(f"A  In-sample (R$^2$ = {r2_score(y, p):.3f}, MAE = {mean_absolute_error(y, p):.1f} K)", loc="left", fontsize=10)
        ax[1].errorbar(loto.heldout_K, loto.pred_mean_K, yerr=loto.pred_sd_K, fmt="o", c="#c05621", capsize=3)
        ax[1].set_title(f"B  Leave-one-temperature-out (mean MAE = {loto.MAE_K.mean():.1f} K)", loc="left", fontsize=10)
        for x_ in ax: x_.plot([230, 470], [230, 470], "k--", lw=.8); x_.set_xlabel("Simulation temperature (K)"); x_.set_ylabel("Predicted temperature (K)"); x_.set_xlim(230, 470); x_.set_ylim(230, 470)
        fig.tight_layout(); fig.savefig(OUT / "fig4_corrected.png", dpi=300); fig.savefig(OUT / "fig4_corrected.pdf"); plt.close(fig)
        fig, ax = plt.subplots(figsize=(5.5, 3.2)); s_ = shp.iloc[::-1]
        ax.barh(s_.Descriptor, s_.share_pct, color="#2b6cb0"); [ax.text(v + .8, i, f"{v:.1f}%", va="center", fontsize=8) for i, v in enumerate(s_.share_pct)]
        ax.set_xlabel("Share of mean |SHAP| (%)"); ax.set_xlim(0, max(s_.share_pct) * 1.18); fig.tight_layout()
        fig.savefig(OUT / "fig3c_shap_corrected.png", dpi=300); fig.savefig(OUT / "fig3c_shap_corrected.pdf"); plt.close(fig)
run(FEAT, "main", True)
run({**{k: v for k, v in FEAT.items() if k != "hbonds_backbone"}, "hbonds": "Total H-bonds"}, "sensitivity: total H-bonds", False)
print(f"\nWrote tables and figures to {OUT}/")
