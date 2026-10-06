"""Recompute Figures 12–13 and Table F2 from included aggregate measurements.

This starts after graph filtering and persistent homology. Counts, lifetime
sums and survival counts are inputs; all normalizations and paired summaries
are calculated here. Run: python comparisons.py
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "output"
METHOD = "topological_thresholding"
LABELS = {
    METHOD: "Topological thresholding", "min_df2": "Min. frequency 2",
    "min_df5": "Min. frequency 5", "eyeballed_l10_u150k": "Frequency 10–150,000",
    "weight_ge_2": "Edge weight ≥ 2", "weight_ge_5": "Edge weight ≥ 5",
    "disparity": "Disparity", "lans": "LANS", "marginal_likelihood": "Marginal likelihood",
    "hypergeometric": "Hypergeometric", "noise_corrected_variance": "Noise-corrected",
}
ALTERNATIVES = list(LABELS)[1:]
METRICS = ["H1_finite_per_1k_edges", "H2_finite_per_1k_triangles",
           "H1_total_finite_persistence_per_1k_edges", "H2_total_finite_persistence_per_1k_triangles"]
TITLES = ["H1 feature density", "H2 feature density", "H1 total persistence", "H2 total persistence"]
COLORS = {METHOD: "#FF6F61", "marginal_likelihood": "#2A9D8F", "noise_corrected_variance": "#264653"}


def calculate_metrics(networks):
    """Keep zeros with positive denominators; leave zero denominators undefined."""
    result = networks.copy()
    for k, denominator in [(1, "edges"), (2, "triangles")]:
        den = result[denominator].where(result[denominator] > 0)
        result[f"H{k}_finite_per_1k_{denominator}"] = 1000 * result[f"H{k}_finite"] / den
        result[f"H{k}_total_finite_persistence_per_1k_{denominator}"] = (
            1000 * result[f"H{k}_total_finite_persistence_years"] / den)
    return result


def summarize(measures):
    """Pair within fields, using positive pairs for ratios and defined pairs for wins."""
    rows = []
    for metric in METRICS:
        values = measures.pivot(index="field_code", columns="method", values=metric)
        for alternative in ALTERNATIVES:
            pairs = values[[METHOD, alternative]].replace([np.inf, -np.inf], np.nan).dropna()
            ours, theirs = pairs[METHOD], pairs[alternative]
            positive = (ours > 0) & (theirs > 0)
            logs = np.log2(ours[positive] / theirs[positive])
            ratios = 2 ** logs.quantile([0.5, 0.25, 0.75]).to_numpy()
            tied = np.abs(ours - theirs) <= 1e-12 + 1e-9 * np.abs(theirs)
            rows.append(dict(metric=metric, alternative=alternative,
                             median_ratio=ratios[0], q25=ratios[1], q75=ratios[2],
                             wins=int(((ours > theirs) & ~tied).sum()), ties=int(tied.sum()),
                             losses=int(((ours < theirs) & ~tied).sum()),
                             n_R=int(positive.sum()), n_P=len(pairs)))
    return pd.DataFrame(rows)


def plot_applied_math(measures, survival):
    """Figure 12: normalized aggregates and finite-lifetime survival profiles."""
    applied = measures[measures.field_code == "102"].set_index("method").loc[list(LABELS)]
    fig, axes = plt.subplots(3, 2, figsize=(10, 12), layout="constrained")
    for panel, (ax, metric, title) in enumerate(zip(axes.flat, METRICS, TITLES)):
        ax.hlines(range(11), 0, applied[metric], color="#CCCCCC", lw=0.7)
        for y, method in enumerate(LABELS):
            marker = "D" if method == METHOD else "s" if method == "noise_corrected_variance" else "o"
            ax.scatter(applied.loc[method, metric], y, marker=marker,
                       s=28 if method == METHOD else 23, color=COLORS.get(method, "#CCCCCC"), zorder=3)
        ax.set_yticks(range(11), list(LABELS.values()) if panel % 2 == 0 else [""] * 11, fontsize=9)
        ax.invert_yaxis()
        units = "Features" if panel < 2 else "Years"
        ax.set_xlabel(f"{units} per 1,000 " + ("edges" if panel % 2 == 0 else "triangles"))
        ax.set_title(f"({chr(97 + panel)}) {title}", fontsize=11, loc="left")
        ax.set_axisbelow(True)
        ax.grid(axis="x", alpha=0.2)
    for k, ax in enumerate(axes[2], 1):
        # Plot muted alternatives first, leaving the three highlighted curves visible.
        for method in [m for m in LABELS if m not in COLORS] + list(COLORS):
            curve = survival[(survival.method == method) & (survival.dimension == k)]
            curve = curve.sort_values("min_lifetime_years")
            density = 1000 * curve.finite_count / curve.denominator.where(curve.denominator > 0)
            ax.step(curve.min_lifetime_years, density, where="pre",
                    color=COLORS.get(method, "#CCCCCC"),
                    lw=1.9 if method == METHOD else 1.25 if method in COLORS else 0.65,
                    ls="--" if method == "noise_corrected_variance" else "-",
                    label=LABELS[method] if method in COLORS else None)
        ax.set(xlim=(0, 40), ylim=(0, None), xlabel="Minimum finite lifetime (years)",
               ylabel=f"Finite H{k} features per 1,000 " + ("edges" if k == 1 else "triangles"))
        ax.set_title(f"({chr(100 + k)}) H{k} lifetime profile", fontsize=11, loc="left")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Applied Mathematics comparison · Figure 12\n"
                 "Finite features and lifetimes in years; eight other filters shown in gray", fontsize=14)
    fig.savefig(OUTPUT / "applied_math_comparison.png", dpi=170)
    plt.close(fig)


def plot_across_fields(summary):
    """Figure 13: medians; Table F2's intervals and pair counts are in the CSV."""
    medians = summary.pivot(index="alternative", columns="metric", values="median_ratio")
    medians = medians.loc[ALTERNATIVES, METRICS]
    counts = summary.pivot(index="alternative", columns="metric", values="n_R").loc[ALTERNATIVES, METRICS]
    fig, ax = plt.subplots(figsize=(11, 7), layout="constrained")
    cmap = LinearSegmentedColormap.from_list("comparison", ["#264653", "#FFFFFF", "#FF6F61"])
    plot = ax.imshow(np.log2(medians), cmap=cmap, vmin=-4, vmax=4, aspect="auto")
    ax.set_yticks(range(10), [LABELS[m] for m in ALTERNATIVES])
    ax.set_xticks(range(4), ["H1 feature\ndensity", "H2 feature\ndensity",
                             "H1 normalized\ntotal persistence", "H2 normalized\ntotal persistence"])
    ax.tick_params(top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=9)
    for i in range(10):
        for j in range(4):
            ratio = medians.iloc[i, j]
            ax.text(j, i, f"{ratio:.2f}×\n(n={counts.iloc[i, j]})", ha="center", va="center",
                    fontsize=10, color="#2B2B2B")
    colorbar = fig.colorbar(plot, ax=ax, orientation="horizontal", shrink=0.75,
                            pad=0.08, ticks=[-4, -2, 0, 2, 4])
    colorbar.ax.set_xticklabels(["1/16×", "1/4×", "1× = equal", "4×", "16×"])
    colorbar.set_label("Median topological-thresholding / alternative ratio")
    fig.suptitle("Cross-field comparison · Figure 13\n"
                 "Positive pairs determine ratios; Table F2 CSV also includes defined zero outcomes", fontsize=13)
    fig.savefig(OUTPUT / "cross_field_comparison.png", dpi=170)
    plt.close(fig)


def main():
    OUTPUT.mkdir(exist_ok=True)
    networks = pd.read_csv(HERE / "data/comparison_networks.csv", dtype={"field_code": str})
    if networks.duplicated(["field_code", "method"]).any():
        raise ValueError("Each field and method must have one measurement row.")
    measures = calculate_metrics(networks)
    summary = summarize(measures)
    measures[measures.field_code == "102"].to_csv(OUTPUT / "applied_math_comparison.csv", index=False)
    summary.to_csv(OUTPUT / "table_F2.csv", index=False)
    survival = pd.read_csv(HERE / "data/applied_math_survival.csv")
    plot_applied_math(measures, survival)
    plot_across_fields(summary)
    usable = measures.loc[(measures.method == METHOD) & measures[METRICS].notna().any(axis=1), "field_code"]
    print(f"Recomputed {len(measures)} network measurements from {networks.field_code.nunique()} fields.")
    print(f"{usable.nunique()} fields contribute defined comparisons; saved all {len(summary)} Table F2 rows.")
    print("Saved two CSV tables and the Figure 12 and Figure 13 PNGs in output/.")


if __name__ == "__main__":
    main()
