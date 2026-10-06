"""Recompute Applied Mathematics selections and variance from saved diagrams.

Figures 8, 9 and 11: images -> local-change scores -> constrained selections;
the variance calculation is a separate assessment of the selected networks.
"""
from itertools import product
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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
    """One overview of the three reproduced analyses; axes show grid indices."""
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), layout="constrained")
    highlighted = selections[selections.kind == "highlighted"]
    for ax, metric, title in zip(axes.flat, ["score", "features", "variance", "sweep"],
                                ["Local-change score · Figure 8a", "H1 + H2 features · Figure 8b",
                                 "Higher-order variance · Figure 11", "Constraint sweep · Figure 9"]):
        if metric == "sweep":
            sweep = selections[selections.kind == "sweep"].dropna(subset=["score"])
            counts = sweep.groupby(["upper_index", "lower_index"]).size()
            ax.scatter(counts.index.get_level_values(1), counts.index.get_level_values(0),
                       s=counts.to_numpy() * 16, color="#db7866", alpha=.8)
            for (i, j), count in counts.nlargest(3).items():
                ax.annotate(str(count), (j, i), ha="center", va="center", fontsize=9)
            ax.text(.02, .97, f"{len(sweep)} feasible settings; {len(counts)} selected cells\nBubble area = selection count",
                    transform=ax.transAxes, va="top", fontsize=9)
        else:
            values = grid.H1 + grid.H2 if metric == "features" else grid[metric]
            matrix = values.to_numpy().reshape(16, 9)
            matrix = np.ma.masked_invalid(matrix)
            if metric != "features":
                matrix[:, 0] = np.ma.masked
            shown = ax.imshow(matrix, origin="lower", aspect="auto", cmap="viridis",
                              norm=LogNorm(vmin=float(matrix[matrix > 0].min()), vmax=float(matrix.max())))
            colorbar = fig.colorbar(shown, ax=ax, shrink=.8)
            colorbar.ax.tick_params(labelsize=8)
            for row, marker in zip(highlighted.itertuples(), ["*", "+", "o"]):
                ax.scatter(row.lower_index, row.upper_index, marker=marker, s=90,
                           color="white", linewidth=1.5)
        ax.set(title=title, xlabel="Lower frequency threshold (articles)",
               ylabel="Upper frequency threshold (articles)", xlim=(-.5, 8.5), ylim=(-.5, 15.5))
        lower = grid.groupby("lower_index").lower.first()
        upper = grid.groupby("upper_index").upper.first()
        ax.set_xticks(range(9), [f"{x:.1f}" for x in lower], rotation=45)
        ax.set_yticks([0, 3, 6, 9, 12, 15], [f"{upper[i]:,.0f}" for i in [0, 3, 6, 9, 12, 15]])
    fig.suptitle("Applied Mathematics: selections from the saved persistence diagrams\n"
                 "Markers: star = 25/25; plus = 50/25; dot = 75/50 constraints", fontsize=13)
    fig.savefig(OUTPUT / "applied_math.png", dpi=170)
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
