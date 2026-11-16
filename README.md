# Loop-first thermal unfolding of the p53 DNA-binding domain and conservation of its heat-sensitive loops across vertebrates

Supplementary data and analysis code for the manuscript of the same title (Mir, Abbas, Abbas & Lee; submitted to *Animals*, MDPI).

## Contents
| Folder | What it holds |
|---|---|
| `simulation_inputs/` | GROMACS inputs: minimisation/equilibration `.mdp`, topology (`topol.top`, `posre.itp`), equilibrated start (`npt.gro`), 1-ns pre-production and 100-ns production `.mdp` (regenerated from the run `.tpr` files with `gmx dump -om`) for 248, 273, 310, 350, 400 and 450 K |
| `data/descriptors/` | Per-frame structural descriptors (`descriptors_corrected.csv`; 1001 frames per temperature) and per-residue native-contact fraction (`per_residue_Q_bhe45.csv`) |
| `data/orthologues/` | 34 reviewed vertebrate p53 sequences from UniProtKB, MAFFT v7.505 alignment, extracted DBDs and key positions |
| `results/tables/` | Per-region, per-residue, correlation, LOTO and SHAP tables |
| `results/figures/` | Main figures (PNG/PDF) |
| `results/supplementary/` | Supplementary Figures S1-S2 and `Supplementary_Tables_S1-S4.xlsx` (also as CSV) |
| `scripts/` | Analysis scripts (run from the project root in the order below) |

## System
Human p53 DBD (UniProt P04637 residues 91-289; PDB 2XWR chain A), zinc-free, OPLS-AA/L force field, SPC/E water, 3 Cl⁻, GROMACS 2023.3. One 100-ns trajectory per temperature after a 1-ns pre-production run at that temperature.

## Analysis order
1. `extract_features_corrected.py` – PBC-corrected descriptors; Best–Hummer–Eaton native contacts (heavy atoms, |i−j|>3, r0<4.5 Å, β=5 Å⁻¹, λ=1.8)
2. `fix_hbonds_corrected.py` – hydrogen bonds with explicit donor/acceptor selections (backbone and total)
3. `rf_shap_corrected.py` – correlations, random forest with leave-one-temperature-out validation, SHAP
4. `dbd_extract.py`, `heat_conservation.py`, `divergence_structure.py` – orthologue alignment mapping, conservation, solvent accessibility and statistics
5. `figs_animals.py`, `fig6_comparative.py`, `fig_structure_pymol.py` + `fig_structure_panel.py`, `supp_materials.py` – figures and supplementary tables

Software: Python 3.12, MDAnalysis 2.10.0, scikit-learn 1.7.2, shap 0.50.0, SciPy 1.17.0, NumPy 2.3.5, pandas 2.3.3, Biopython 1.86, Matplotlib, PyMOL 3.2.

## Trajectories
The six 100-ns trajectories (~1.3 GB each) exceed GitHub's file-size limit and will be deposited on Zenodo; the DOI will be added here.

## License
Code: MIT. Data, tables and figures: CC BY 4.0. See `LICENSE`.

## Contact
Syed Raza Abbas (GitHub: raza-create); corresponding authors: Zeeshan Abbas, Seung Won Lee.
