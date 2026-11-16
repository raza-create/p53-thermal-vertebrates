"""PyMOL renders for Figure 2. Fixes: flat lighting (white = no loss reads as white, not grey),
per-residue colours (no blending along strands), larger labels on a white background."""
import pandas as pd
from pymol import cmd
m = pd.read_csv("results/corrected/per_residue_heat_conservation_rsa.csv")
O = "results/corrected/figs_animals"
cmd.load("2xwr.cif", "p53"); cmd.remove("not (chain A and (polymer.protein or resn ZN))"); cmd.remove("hydro")
cmd.bg_color("white")
for k, v in [("ray_opaque_background", 1), ("antialias", 2), ("ray_shadow", 0), ("depth_cue", 0), ("ray_trace_fog", 0),
             ("ambient", 0.55), ("direct", 0.45), ("reflect", 0.15), ("specular", 0), ("light_count", 2),
             ("ray_trace_mode", 1), ("ray_trace_gain", 0.05), ("cartoon_discrete_colors", 1),
             ("cartoon_fancy_helices", 1), ("cartoon_highlight_color", -1), ("stick_radius", 0.22)]:
    cmd.set(k, v)
for k, v in [("label_size", -3.2), ("label_color", "black"), ("label_font_id", 7), ("float_labels", 1),
             ("label_bg_color", "white"), ("label_bg_transparency", 0.15), ("label_position", (0, 0, 8))]:
    cmd.set(k, v)
cmd.hide("everything"); cmd.show("cartoon", "polymer"); cmd.show("spheres", "resn ZN"); cmd.set("sphere_scale", 0.7, "resn ZN")
cmd.alter("all", "b=0.0")
for r in m.itertuples(): cmd.alter(f"polymer and resi {r.resid}", f"b={r.dQ350}")
KEY = [175, 176, 179, 186, 229, 248, 273, 282]; aa = {r.resid: r.human for r in m.itertuples()}
cmd.show("sticks", "polymer and resi " + "+".join(map(str, KEY)) + " and not name N+C+O")
cmd.orient("polymer"); view = cmd.get_view()
def shoot(tag):
    cmd.set_view(view); cmd.png(f"{O}/_struct_{tag}_front.png", width=2000, height=1650, dpi=300, ray=1)
    cmd.turn("y", 180); cmd.png(f"{O}/_struct_{tag}_back.png", width=2000, height=1650, dpi=300, ray=1)
cmd.spectrum("b", "white_red", "polymer", minimum=0.0, maximum=0.9); cmd.color("grey40", "resn ZN")
for r in KEY: cmd.label(f"polymer and resi {r} and name CA", f'"{aa.get(r, "?")}{r}"')
shoot("heat")
cmd.label("all", "''")
cmd.set_color("inv", [0x2a/255, 0x78/255, 0xd6/255]); cmd.set_color("div", [0xeb/255, 0x68/255, 0x34/255]); cmd.set_color("oth", [0xe5/255, 0xe5/255, 0xe5/255])
cmd.color("oth", "polymer")
cmd.color("inv", "polymer and resi " + "+".join(str(r) for r in m[m.invariant].resid))
cmd.color("div", "polymer and resi " + "+".join(str(r) for r in m[m.divergent].resid))
cmd.color("grey40", "resn ZN")
shoot("cons"); print("PyMOL renders written")
