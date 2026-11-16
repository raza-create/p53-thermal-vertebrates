"""Assemble Figure 2: colorbar in its own row, whitespace cropped, sparse ticks, legend row."""
import numpy as np, matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
O = "results/corrected/figs_animals"; INK, INK2 = "#0b0b0b", "#52514e"
def crop(path, pad=20, thr=0.97):
    im = mpimg.imread(path)[..., :3]; ys, xs = np.where((im < thr).any(axis=2))
    return im[max(ys.min()-pad, 0):min(ys.max()+pad, im.shape[0]), max(xs.min()-pad, 0):min(xs.max()+pad, im.shape[1])]
fig = plt.figure(figsize=(7.1, 7.6))
gs = fig.add_gridspec(5, 2, height_ratios=[1.0, 0.035, 0.10, 1.0, 0.06], hspace=0.08, wspace=0.04, left=0.02, right=0.98, top=0.95, bottom=0.03)
titles = {"heat": "A   Heat sensitivity (native contacts lost, 248\u2192350 K)", "cons": "B   Conservation across 33 vertebrates"}
for tag, r in {"heat": 0, "cons": 3}.items():
    for j, side in enumerate(["front", "back"]):
        ax = fig.add_subplot(gs[r, j]); ax.imshow(crop(f"{O}/_struct_{tag}_{side}.png")); ax.axis("off")
        if j == 0: ax.set_title(titles[tag], loc="left", fontsize=9, fontweight="bold", color=INK, pad=6)
        else: ax.set_title("rotated 180\u00b0", loc="right", fontsize=7.5, color=INK2, pad=6)
cax = fig.add_subplot(gs[1, :].subgridspec(1, 3, width_ratios=[1, 1.4, 1])[0, 1])
cb = fig.colorbar(plt.cm.ScalarMappable(cmap=LinearSegmentedColormap.from_list("wr", ["#ffffff", "#ff0000"]), norm=plt.Normalize(0, 0.9)),
                  cax=cax, orientation="horizontal", ticks=[0, 0.3, 0.6, 0.9])
cb.set_label("Native-contact fraction lost", fontsize=7, color=INK2, labelpad=2); cb.ax.tick_params(labelsize=7, colors=INK2, length=2); cb.outline.set_linewidth(0.5)
lax = fig.add_subplot(gs[4, :]); lax.axis("off")
lax.legend(handles=[Patch(color="#2a78d6", label="Invariant"), Patch(color="#eb6834", label="Endotherm/ectotherm divergent"),
                    Patch(color="#e5e5e5", label="Other variable"), Patch(color="#666666", label="Zn\u00b2\u207a (crystal structure)")],
           loc="center", ncol=4, fontsize=7, frameon=False, handlelength=1.4, columnspacing=1.2)
fig.savefig(f"{O}/fig5_structure.png", dpi=300, facecolor="white"); fig.savefig(f"{O}/fig5_structure.pdf", facecolor="white")
print(f"Wrote {O}/fig5_structure.(png|pdf)")
