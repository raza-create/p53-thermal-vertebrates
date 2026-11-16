#!/usr/bin/env python3
"""Corrected descriptors for the six OPLS-AA (zinc-free) temperature trajectories.
PBC: protein made whole every frame. Common reference: frame 0 of the 248 K run.
Q: Best-Hummer-Eaton (heavy atoms, |i-j|>3, r0<4.5 A, beta=5/A, lambda=1.8); r0<8 A kept as sensitivity.
H-bonds: D-A < 3.5 A and D-H-A > 120 deg. Usage: python3 scripts/extract_features_corrected.py [stride]
"""
import sys, time, numpy as np, pandas as pd, MDAnalysis as mda
from pathlib import Path
from MDAnalysis import transformations as trans
from MDAnalysis.analysis import rms
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from MDAnalysis.lib.distances import capped_distance, calc_bonds
D = Path("data/structures/long_simulations"); OUT = Path("data/features_corrected"); OUT.mkdir(parents=True, exist_ok=True)
TS = [248, 273, 310, 350, 400, 450]; STRIDE = int(sys.argv[1]) if len(sys.argv) > 1 else 10
BETA, LAM, OFFSET = 5.0, 1.8, 90   # topology resid 1 = UniProt 91 (checked below)

def load(T):
    u = mda.Universe(str(D / f"md_100ns_{T}K.tpr"), str(D / f"md_100ns_{T}K.xtc"))
    prot = u.select_atoms("protein"); u.trajectory.add_transformations(trans.unwrap(prot))
    return u, prot, prot.select_atoms("name CA"), prot.select_atoms("not name H*")

def pairs_for(heavy, pos, cut):
    p, d0 = capped_distance(pos, pos, max_cutoff=cut)
    ri = heavy.resindices
    k = (p[:, 0] < p[:, 1]) & (np.abs(ri[p[:, 0]] - ri[p[:, 1]]) > 3)
    return p[k], d0[k]

ru, rprot, rca, rheavy = load(248); ru.trajectory[0]
ref_ca, ref_heavy = rca.positions.copy(), rheavy.positions.copy()
seq = "".join(r.resname[:3] for r in rprot.residues[:6])
print(f"first residues: {seq} (expect TRP PRO LEU SER SER SER = UniProt 91-96)")
P45, D45 = pairs_for(rheavy, ref_heavy, 4.5); P80, D80 = pairs_for(rheavy, ref_heavy, 8.0)
nres = rprot.n_residues; ri = rheavy.resindices - rprot.residues[0].resindex
print(f"native contacts: {len(P45)} (4.5 A), {len(P80)} (8 A); residues {nres}")
rows, perres = [], {}
for T in TS:
    t0 = time.time(); u, prot, ca, heavy = load(T)
    assert heavy.n_atoms == rheavy.n_atoms, "heavy-atom count differs from reference"
    nzn = len(u.select_atoms("resname ZN or name ZN"))
    hba = HBA(universe=u, d_a_cutoff=3.5, d_h_a_angle_cutoff=120.0)
    hba.hydrogens_sel = hba.guess_hydrogens("protein"); hba.acceptors_sel = hba.guess_acceptors("protein")
    hba.run(step=STRIDE); hb = hba.count_by_time()
    qsum = np.zeros(nres); qcnt = np.bincount(ri[P45[:, 0]], minlength=nres) + np.bincount(ri[P45[:, 1]], minlength=nres)
    for k, ts in enumerate(u.trajectory[::STRIDE]):
        x = heavy.positions
        q45 = 1.0 / (1.0 + np.exp(BETA * (calc_bonds(x[P45[:, 0]], x[P45[:, 1]]) - LAM * D45)))
        q80 = 1.0 / (1.0 + np.exp(BETA * (calc_bonds(x[P80[:, 0]], x[P80[:, 1]]) - LAM * D80)))
        qsum += np.bincount(ri[P45[:, 0]], q45, nres) + np.bincount(ri[P45[:, 1]], q45, nres)
        xp = prot.positions - prot.positions.mean(0); ev = np.sort(np.linalg.eigvalsh(xp.T @ xp / len(xp)))
        rows.append(dict(temperature=T, frame=ts.frame, time_ns=ts.time / 1000,
                         rmsd=rms.rmsd(ca.positions, ref_ca, center=True, superposition=True),
                         rg=prot.radius_of_gyration(), end_to_end=np.linalg.norm(ca.positions[-1] - ca.positions[0]),
                         hbonds=int(hb[k]), Q_bhe45=q45.mean(), Q_8A=q80.mean(), asphericity=ev[2] - 0.5 * (ev[0] + ev[1])))
    perres[T] = qsum / np.maximum(qcnt * (k + 1), 1)
    print(f"{T} K: {k+1} frames, Zn atoms = {nzn}, {time.time()-t0:.0f}s")
df = pd.DataFrame(rows); df.to_csv(OUT / "descriptors_corrected.csv", index=False)
pr = pd.DataFrame(perres); pr.insert(0, "resid_uniprot", np.arange(nres) + 1 + OFFSET); pr.to_csv(OUT / "per_residue_Q_bhe45.csv", index=False)
old = pd.read_csv("data/features_100ns/all_features_100ns.csv").groupby("temperature").native_contacts.mean()
cols = ["Q_bhe45", "Q_8A", "rmsd", "rg", "end_to_end", "hbonds", "asphericity"]
s = df.groupby("temperature")[cols].agg(["mean", "std"]).round(3); s[("old_Q_csv", "mean")] = old.round(3)
print("\nPer-temperature mean (SD):"); print(s.to_string())
print("\nPearson r with temperature (frames, descriptive only):")
print(df[cols].corrwith(df.temperature).round(3).to_string())
print(f"\nWrote {OUT}/descriptors_corrected.csv and per_residue_Q_bhe45.csv")
