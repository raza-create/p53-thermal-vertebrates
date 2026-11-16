#!/usr/bin/env python3
"""Where do endotherm/ectotherm-divergent positions sit? Loops vs rest, surface vs buried (2XWR chain A), per-species identity."""
import numpy as np, pandas as pd
from scipy.stats import fisher_exact
from Bio.PDB import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
OUT = "results/corrected"
m = pd.read_csv(f"{OUT}/per_residue_heat_conservation.csv")
MAXASA = dict(A=129, R=274, N=195, D=193, C=167, Q=225, E=223, G=104, H=224, I=197, L=201, K=236, M=224, F=240, P=159, S=155, T=172, W=285, Y=263, V=174)  # Tien et al. 2013, theoretical
THREE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I", LEU="L", LYS="K", MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
s = MMCIFParser(QUIET=True).get_structure("2xwr", "2xwr.cif"); ch = s[0]["A"]
for r in [r for r in ch if r.id[0] != " " or r.get_resname() not in THREE]: ch.detach_child(r.id)
ShrakeRupley().compute(ch, level="R")
rsa = {r.id[1]: r.sasa / MAXASA[THREE[r.get_resname()]] for r in ch}
m["rsa"] = m.resid.map(rsa); print(f"RSA available for {m.rsa.notna().sum()}/{len(m)} residues (2XWR chain A, author numbering)")
mm = m.dropna(subset=["rsa"]).copy(); mm["surface"] = mm.rsa >= 0.20; mm["loops"] = mm.region.isin(["L2 (167-194)", "L3 (237-250)"])
def fx(a, b, lab):
    t = pd.crosstab(mm[a], mm[b]).reindex(index=[True, False], columns=[True, False], fill_value=0)
    o, p = fisher_exact(t.values); print(f"{lab}: {t.loc[True, True]}/{t.loc[True].sum()} vs {t.loc[False, True]}/{t.loc[False].sum()}  OR={o:.2f} p={p:.3g}"); return p
print("\nFisher tests (fraction divergent | fraction invariant):")
fx("loops", "divergent", "divergent in L2+L3 vs rest"); fx("surface", "divergent", "divergent at surface (RSA>=0.2) vs buried")
fx("loops", "invariant", "invariant in L2+L3 vs rest"); fx("surface", "invariant", "invariant at surface vs buried")
print("\nPer-region divergence vs rest (BH-adjusted):")
rows = []
for reg in sorted(mm.region.unique()):
    x = mm.region == reg; t = [[(x & mm.divergent).sum(), (x & ~mm.divergent).sum()], [(~x & mm.divergent).sum(), (~x & ~mm.divergent).sum()]]
    o, p = fisher_exact(t); rows.append(dict(region=reg, n=x.sum(), divergent=t[0][0], OR=o, p=p))
r = pd.DataFrame(rows).sort_values("p"); n = len(r); r["p_BH"] = np.minimum.accumulate((r.p * n / np.arange(1, n + 1))[::-1])[::-1].clip(upper=1)
print(r.round(3).to_string(index=False))
print(f"\nheat-sensitive top-20: surface {mm.nlargest(20, 'dQ350').surface.sum()}/20; overall surface {mm.surface.sum()}/{len(mm)}")
def read_fa(p):
    d, k = {}, None
    for l in open(p):
        l = l.strip()
        if not l: continue
        if l.startswith(">"): k = l[1:].split("|")[1]; d[k] = []
        else: d[k].append(l)
    return {a: "".join(b) for a, b in d.items()}
aln = read_fa("species/vert_p53_aln.fasta"); H = aln["P04637"]; p2c, k = {}, 0
for c, a in enumerate(H):
    if a != "-": k += 1; p2c[k] = c
sets = {"top20_heat": m.nlargest(20, "dQ350").resid, "L2+L3": m[m.region.isin(["L2 (167-194)", "L3 (237-250)"])].resid,
        "divergent": m[m.divergent].resid, "DBD_all": m.resid}
SP = {"P04637": "Human", "P02340": "Mouse", "Q29537": "Dog", "P41685": "Cat", "P10360": "Chicken", "P07193": "Xenopus", "P79734": "Zebrafish"}
print("\nPer-species identity to human (%):")
out = []
for acc, name in SP.items():
    row = dict(species=name)
    for lab, pos in sets.items():
        cols = [p2c[p] for p in pos]; nn = [c for c in cols if aln[acc][c] != "-"]
        row[lab] = round(100 * sum(aln[acc][c] == H[c] for c in nn) / max(1, len(nn)), 1)
    out.append(row)
print(pd.DataFrame(out).to_string(index=False))
mm.to_csv(f"{OUT}/per_residue_heat_conservation_rsa.csv", index=False); print(f"\nWrote {OUT}/per_residue_heat_conservation_rsa.csv")
