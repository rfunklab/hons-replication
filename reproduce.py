"""Recompute Applied Mathematics selections and variance from saved diagrams.

Figures 8, 9 and 11: images -> local-change scores -> constrained selections;
the variance calculation is a separate assessment of the selected networks.
"""
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import numpy as np
import pandas as pd

import hons

HERE = Path(__file__).resolve().parent
DATA, OUTPUT = HERE / "data", HERE / "output"
HIGHLIGHTED = [("star", 25, 25), ("plus", 50, 25), ("dot", 75, 50)]


def calculate():
    grid = pd.read_csv(DATA / "applied_math_grid.csv").sort_values(
        ["upper_index", "lower_index"]).reset_index(drop=True)
    lower = grid.groupby("lower_index").lower.first().to_numpy()
    upper = grid.groupby("upper_index").upper.first().to_numpy()
    keys = list(zip(grid.upper_index, grid.lower_index))
    images, h1, h2 = {}, {}, {}
    with np.load(DATA / "applied_math_diagrams.npz", allow_pickle=False) as data:
        if set(data.files) != {f"{i}_{j}" for i, j in keys}:
            raise ValueError("The diagram archive must cover every grid cell.")
        for n, cell in enumerate(keys, 1):
            values = data[f"{cell[0]}_{cell[1]}"]
            dim, birth, death = values.T
            valid = (np.isin(dim, [1, 2]) & np.isfinite(birth) &
                     (birth >= 0) & (birth <= 1) & (death > birth) &
                     ((np.isfinite(death) & (death <= 1)) | np.isposinf(death)))
            if not valid.all():
                raise ValueError(f"Invalid persistence interval at {cell}.")
            diagram = pd.DataFrame(values, columns=["dimension", "birth", "death"])
            images[cell] = hons.rho(diagram)
            h1[cell], h2[cell] = int((dim == 1).sum()), int((dim == 2).sum())
            if n % 9 == 0:
                print(f"Persistence images: {n}/{len(keys)}", flush=True)
    scores = hons.objective(images, lower, upper)
    grid["H1"], grid["H2"] = [h1[k] for k in keys], [h2[k] for k in keys]
    grid["score"] = [scores.get(k, np.nan) for k in keys]
    grid["selectable"] = [k in scores for k in keys]

    # Equation (higher-order variance): mean over all 144 images, n=400 pixels.
    vectors = np.stack([images[k] for k in keys])
    difference = vectors[:, :400] - vectors[:, 400:]
    grid["variance"] = ((difference - difference.mean(axis=0)) ** 2).sum(axis=1) / 399
    eligible = grid[grid.selectable]
    cutoff1, cutoff2 = np.percentile(eligible.H1, 50), np.percentile(eligible.H2, 25)
    grid["primary_feasible"] = grid.selectable & (grid.H1 >= cutoff1) & (grid.H2 >= cutoff2)
    grid["variance_rank_selectable"] = grid.loc[grid.selectable, "variance"].rank(method="min")
    grid["variance_rank_primary"] = grid.loc[grid.primary_feasible, "variance"].rank(method="min")

    settings = [("highlighted", name, p, q) for name, p, q in HIGHLIGHTED]
    settings += [("sweep", "", p, q) for p, q in product(range(1, 92, 10), repeat=2)]
    selections = []
    for kind, name, p, q in settings:
        cell, ties = hons.select(scores, h1, h2, p, q)
        c1, c2 = np.percentile(eligible.H1, p), np.percentile(eligible.H2, q)
        selections.append(dict(
            kind=kind, configuration=name, percentile_H1=p, percentile_H2=q,
            minimum_H1=c1, minimum_H2=c2,
            feasible_cells=int(((eligible.H1 >= c1) & (eligible.H2 >= c2)).sum()),
            upper_index=cell[0] if cell else np.nan, lower_index=cell[1] if cell else np.nan,
            lower=lower[cell[1]] if cell else np.nan, upper=upper[cell[0]] if cell else np.nan,
            score=scores[cell] if cell else np.nan, number_of_ties=len(ties)))
    return grid, pd.DataFrame(selections)


def plot(grid, selections):
    """Draw Figures 8, 9 and 11 on the manuscript's selectable grid."""
    plt.rcParams.update({"font.size": 11, "axes.labelsize": 11,
                         "xtick.labelsize": 10, "ytick.labelsize": 10})
    palette = LinearSegmentedColormap.from_list("charcoal", ["#eff2f2", "#264653"])
    lower = grid.groupby("lower_index").lower.first().to_numpy()
    upper = grid.groupby("upper_index").upper.first().to_numpy()
    highlighted = selections[selections.kind == "highlighted"]

    def heatmap(ax, metric, title):
        values = grid.H1 + grid.H2 if metric == "features" else grid[metric]
        matrix = values.to_numpy().reshape(len(upper), len(lower))[:, 1:]
        norm = LogNorm(vmin=matrix.min(), vmax=matrix.max()) if metric == "variance" else None
        shown = ax.imshow(matrix, origin="lower", aspect="auto", cmap=palette,
                          norm=norm, interpolation="nearest")
        ax.set(title=title, xlabel="Lower bound ℓ (articles)", ylabel="Upper bound u (articles)")
        ax.set_box_aspect(1)
        ax.set_xticks(range(len(lower) - 1),
                      [f"{x:.1f}" if i in {0, 2, 5, 7} else "" for i, x in enumerate(lower[1:])])
        ax.set_yticks(range(len(upper)),
                      ["No upper\ncutoff" if i == len(upper) - 1 else
                       f"{x:.1f}" if i % 2 == 0 else "" for i, x in enumerate(upper)])
        ax.tick_params(length=2)
        colorbar = ax.figure.colorbar(shown, ax=ax, fraction=.046, pad=.035)
        colorbar.ax.tick_params(labelsize=9)

    def markers(ax):
        for row, marker, size in zip(highlighted.itertuples(), ["*", "P", "o"], [150, 85, 55]):
            ax.scatter(row.lower_index - 1, row.upper_index, marker=marker, s=size,
                       color="#2b2b2b", edgecolors="white", linewidth=.8, zorder=5)

    # The first lower-bound column supplies neighbors but is not selectable.
    marker_key = "H₁/H₂ constraint percentiles: star 25/25 · plus 50/25 · dot 75/50"
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4), layout="constrained")
    for ax, metric, title in zip(axes, ["score", "features"],
                                ["Figure 8a · Local-change score", "Figure 8b · H₁ + H₂ features"]):
        heatmap(ax, metric, title)
        markers(ax)
    fig.supxlabel(marker_key, fontsize=11)
    fig.savefig(OUTPUT / "applied_math_selections.png", dpi=170)
    plt.close(fig)

    sweep = selections[selections.kind == "sweep"].dropna(subset=["score"]).copy()
    sweep["strictness"] = sweep.percentile_H1 + sweep.percentile_H2
    bubbles = sweep.groupby(["upper_index", "lower_index"]).agg(
        count=("strictness", "size"), strictness=("strictness", "mean")).reset_index()
    lo, hi = sweep.strictness.min(), sweep.strictness.max()
    shade = (bubbles.strictness - lo) / (hi - lo)
    fig, axes = plt.subplots(1, 2, figsize=(11, 5.4), layout="constrained")
    for ax, metric, title in zip(axes, ["score", "features"],
                                ["Figure 9a · Local-change score", "Figure 9b · H₁ + H₂ features"]):
        heatmap(ax, metric, title)
        ax.scatter(bubbles.lower_index - 1, bubbles.upper_index, s=bubbles["count"] * 17,
                   c=shade, cmap="gray_r", vmin=-.4, vmax=1.15,
                   edgecolors="#2b2b2b", linewidths=.6, zorder=5)
    for row in bubbles.nlargest(3, "count").itertuples():
        axes[0].annotate(str(row.count), (row.lower_index - 1, row.upper_index),
                         xytext=(11, 0), textcoords="offset points", va="center", fontsize=10)
    fig.supxlabel(f"{len(sweep)} feasible settings · {len(bubbles)} selected cells\n"
                  "Bubble area = number of selections; darker bubbles = stricter constraints", fontsize=11)
    fig.savefig(OUTPUT / "applied_math_sweep.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 6.4), layout="constrained")
    heatmap(ax, "variance", "Figure 11 · Higher-order variance")
    markers(ax)
    fig.supxlabel(marker_key, fontsize=9)
    fig.savefig(OUTPUT / "applied_math_variance.png", dpi=170)
    plt.close(fig)


def main():
    OUTPUT.mkdir(exist_ok=True)
    grid, selections = calculate()
    grid.to_csv(OUTPUT / "applied_math_grid.csv", index=False)
    selections.to_csv(OUTPUT / "applied_math_selections.csv", index=False)
    plot(grid, selections)
    primary = selections[selections.configuration == "plus"].iloc[0]
    selected = grid[(grid.upper_index == primary.upper_index) & (grid.lower_index == primary.lower_index)].iloc[0]
    print(selections[selections.kind == "highlighted"].to_string(index=False))
    print(f"Primary variance rank: {selected.variance_rank_primary:.0f} of "
          f"{grid.primary_feasible.sum()} feasible; {selected.variance_rank_selectable:.0f} of "
          f"{grid.selectable.sum()} selectable.")


if __name__ == "__main__":
    main()
