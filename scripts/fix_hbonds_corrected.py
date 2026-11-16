#!/usr/bin/env python3
"""Recompute intramolecular H-bonds with explicit selections (guess_hydrogens misses OPLS-AA amide H, q = +0.30)."""
import numpy as np, pandas as pd, MDAnalysis as mda
from pathlib import Path
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
D = Path("data/structures/long_simulations"); F = Path("data/features_corrected/descriptors_corrected.csv")
TS = [248, 273, 310, 350, 400, 450]; STRIDE = 10
H_SEL = "protein and name H* and bonded (protein and (name N* or name O*))"
A_SEL = "protein and (name O* or (resname HIS* and name ND1 NE2))"
df = pd.read_csv(F)
if "hbonds_guess" not in df.columns: df = df.rename(columns={"hbonds": "hbonds_guess"})
df = df.drop(columns=[c for c in ["hbonds", "hbonds_backbone"] if c in df.columns])
out = []
for T in TS:
    u = mda.Universe(str(D / f"md_100ns_{T}K.tpr"), str(D / f"md_100ns_{T}K.xtc"))
    if T == 248:
        bbH = u.select_atoms("protein and name H")
        g = HBA(universe=u).guess_hydrogens("protein")
        print(f"backbone amide H: n={bbH.n_atoms}, charge min/max = {bbH.charges.min():.3f}/{bbH.charges.max():.3f}")
        print(f"hydrogens guessed by guess_hydrogens(): {len(u.select_atoms(g))}; explicit polar H: {len(u.select_atoms(H_SEL))}; acceptors: {len(u.select_atoms(A_SEL))}")
    hba = HBA(universe=u, hydrogens_sel=H_SEL, acceptors_sel=A_SEL, d_a_cutoff=3.5, d_h_a_angle_cutoff=120.0)
    hba.run(step=STRIDE)
    hb = hba.results.hbonds; frames = np.arange(0, len(u.trajectory), STRIDE)
    fr, don, acc = hb[:, 0].astype(int), hb[:, 1].astype(int), hb[:, 3].astype(int)
    bb = (u.atoms.names[don] == "N") & (u.atoms.names[acc] == "O")
    tot = pd.Series(fr).value_counts().reindex(frames, fill_value=0).values
    bbc = pd.Series(fr[bb]).value_counts().reindex(frames, fill_value=0).values
    out.append(pd.DataFrame(dict(temperature=T, frame=frames, hbonds=tot, hbonds_backbone=bbc)))
    print(f"{T} K: mean H-bonds {tot.mean():.1f} (backbone {bbc.mean():.1f})")
df = df.merge(pd.concat(out), on=["temperature", "frame"], how="left")
assert df.hbonds.notna().all(), "frame mismatch during merge"
df.to_csv(F, index=False)
cols = ["hbonds_guess", "hbonds", "hbonds_backbone"]
print("\nPer-temperature mean (SD):"); print(df.groupby("temperature")[cols].agg(["mean", "std"]).round(1).to_string())
print("\nPearson r with temperature:"); print(df[cols].corrwith(df.temperature).round(3).to_string())
print(f"\nUpdated {F}")
