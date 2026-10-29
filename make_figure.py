"""Redraw fig_depth_covariate.png from the committed tables under data/ (no release files needed)."""
import numpy as np, pandas as pd, matplotlib as mpl, matplotlib.pyplot as plt
mpl.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
                     "axes.spines.top": False, "axes.spines.right": False})
cells = pd.read_csv("data/per_cell_depth.csv.gz"); ctd = pd.read_csv("data/per_celltype_depth.csv"); pb = pd.read_csv("data/pseudobulk_log1pCP10K.csv.gz")
stages = ["E8.5", "E9.5", "E10.5"]; col = {"E8.5": "#9aa5b1", "E9.5": "#4c72b0", "E10.5": "#c44e52"}
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.3), gridspec_kw=dict(width_ratios=[1, 1.15, 1.25]))
ax = axes[0]
for i, s in enumerate(stages):
    v = cells.loc[cells.stage == s, "umi"].values
    parts = ax.violinplot(np.log10(v), positions=[i], showextrema=False, widths=0.8)
    for b in parts["bodies"]: b.set_facecolor(col[s]); b.set_alpha(0.75); b.set_edgecolor("none")
    med = np.median(v); ax.plot([i - 0.3, i + 0.3], [np.log10(med)] * 2, color="k", lw=1.2)
    ax.text(i, np.log10(med) + 0.05, f"{med/1e3:.1f}k", ha="center", va="bottom", fontsize=7)
ax.set_xticks(range(3)); ax.set_xticklabels(stages); ax.set_yticks([3, 3.5, 4, 4.5]); ax.set_yticklabels(["1k", "3k", "10k", "30k"])
ax.set_ylabel("UMI per cell (recovered counts)"); ax.set_title("E10.5 cells are sequenced ~40% shallower", loc="left"); ax.margins(0.06)
ax = axes[1]; d = ctd.sort_values("ratio_E10.5/E9.5"); y = np.arange(len(d))
ax.hlines(y, 1, d["ratio_E10.5/E9.5"], color="#bbb", lw=1); ax.plot(d["ratio_E10.5/E9.5"], y, "o", color=col["E10.5"], ms=5)
ax.axvline(1, color="k", lw=0.8); ax.axvline(0.607, color=col["E10.5"], lw=0.8, ls="--")
ax.set_yticks(y); ax.set_yticklabels(d["ct"]); ax.set_xlabel("median UMI ratio, E10.5 / E9.5, same label"); ax.set_xlim(0.4, 1.15)
ax.set_title("Drop is within cell type, not composition", loc="left"); ax.text(0.607, len(d) - 0.4, "stage\nmedian 0.61", color=col["E10.5"], fontsize=7, ha="right", va="top")
ax = axes[2]; m = pb.expressed.values; lvl = 0.5 * (pb["pb_E9.5"] + pb["pb_E10.5"]).values
d_meas = (pb["pb_E10.5"] - pb["pb_E9.5"]).values; d_neu = (pb["pb_E10.5"] - pb["pb_E9.5_thinned"]).values
bins = np.quantile(lvl[m], np.linspace(0, 1, 21)); mid = 0.5 * (bins[1:] + bins[:-1]); idx = np.clip(np.digitize(lvl[m], bins) - 1, 0, 19)
ax.scatter(lvl[m], d_meas[m], s=2, color="#d9d9d9", lw=0)
for dd, c, lab in [(d_meas, col["E10.5"], "measured step  ρ = −0.41"), (d_neu, "#4c72b0", "E9.5 thinned to E10.5 depth before the step  ρ = −0.23")]:
    med = np.array([np.median(dd[m][idx == k]) for k in range(20)]); ax.plot(mid, med, "-o", color=c, ms=3.5, lw=1.4, label=lab)
ax.axhline(0, color="k", lw=0.8); ax.set_xlabel("mean expression, log1p CP10K (20 quantile bins)"); ax.set_ylabel("pseudobulk Δ, E10.5 − E9.5")
ax.set_title("Half the level-correlation is depth", loc="left"); ax.set_ylim(-0.6, 0.35); ax.set_xlim(-0.1, 3.6)
ax.legend(frameon=False, loc="lower left", handlelength=1.5, handletextpad=0.4); ax.text(3.5, 0.32, "grey: single genes, measured", ha="right", va="top", fontsize=6, color="#777")
for a, l in zip(axes, "abc"): a.text(-0.12, 1.06, l, transform=a.transAxes, fontweight="bold", fontsize=10)
fig.tight_layout(w_pad=1.2); fig.savefig("fig_depth_covariate.png", dpi=200); print("fig_depth_covariate.png")
