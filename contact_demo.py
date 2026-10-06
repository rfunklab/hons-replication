"""Threshold an observed workplace network by total contact duration.

Each row is one 20-second contact. Aggregate unordered pairs to contact minutes
and their first timestamp; all observed participants remain nodes. HONS selects
the duration cutoff from the candidate grid.
"""
from pathlib import Path
from time import perf_counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

import hons

root = Path(__file__).resolve().parent
output = root / "output"
output.mkdir(exist_ok=True)
raw = pd.read_csv(root / "data/workplace_contacts.tsv", sep=r"\s+", names=["t", "i", "j"])
raw[["i", "j"]] = np.sort(raw[["i", "j"]].to_numpy(), axis=1)
pairs = raw.groupby(["i", "j"]).agg(intervals=("t", "size"), first_seen=("t", "min"))
pairs["minutes"] = pairs.intervals / 3
network = nx.from_pandas_edgelist(pairs.reset_index(), "i", "j", edge_attr=True)
people, counts = np.unique(raw[["i", "j"]].to_numpy(), return_counts=True)
network.add_nodes_from((int(p), {"contact_intervals": int(n)}) for p, n in zip(people, counts))
window = (int(raw.t.min()), int(raw.t.max()))
network.graph.update(source="SocioPatterns workplace contacts", observation_window=window)

# Candidate lower cutoffs span 0–5 minutes in 20-second steps.
lower = np.arange(0, 301, 20) / 60
upper = [float(pairs.minutes.max())]
constraints = (50, 25)
started = perf_counter()
selected, details = hons.threshold(
    network, edge_attribute="minutes", lower=lower, upper=upper,
    filtration="first_seen", filtration_range=window,
    constraints=constraints, return_details=True,
)
elapsed = perf_counter() - started
grid = details["grid"].copy()
q1 = np.percentile(grid.loc[grid.selectable, "H1"], constraints[0])
q2 = np.percentile(grid.loc[grid.selectable, "H2"], constraints[1])
grid["H1_min"], grid["H2_min"] = q1, q2
grid["feasible"] = grid.selectable & (grid.H1 >= q1) & (grid.H2 >= q2)
grid["selected"] = [(i, j) == details["cell"] for i, j in zip(grid.upper_index, grid.lower_index)]
grid.to_csv(output / "contact_grid.csv", index=False)

plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), gridspec_kw={"height_ratios": [1.45, 1]})
positions = nx.spring_layout(network, seed=7)
xy = np.array(list(positions.values()))
for ax, graph, title in zip(axes[0], [network, selected], ["All observed contacts", "HONS-selected contacts"]):
    if graph is None:
        ax.text(.5, .5, "No cutoff satisfies both constraints", ha="center", transform=ax.transAxes)
        ax.set_title(title)
    else:
        isolates = set(nx.isolates(graph))
        nx.draw_networkx_edges(graph, positions, ax=ax, width=.6, alpha=.35, edge_color="#507b89")
        nx.draw_networkx_nodes(graph, positions, ax=ax, node_size=24, linewidths=.3,
                              edgecolors="#50646e", node_color=["#e2e7e9" if n in isolates else "#e0805c" for n in graph])
        ax.set_title(f"{title}\n{len(graph)} participants · {graph.number_of_edges()} edges · {len(isolates)} isolates")
    ax.set_xlim(xy[:, 0].min() - .1, xy[:, 0].max() + .1)
    ax.set_ylim(xy[:, 1].min() - .1, xy[:, 1].max() + .1)
    ax.set_aspect("equal")
    ax.axis("off")

score_ax, count_ax = axes[1]
score_ax.plot(grid.lower, grid.score, "o-", color="#55608e", ms=4, label="HONS score")
score_ax.scatter(grid.loc[grid.feasible, "lower"], grid.loc[grid.feasible, "score"],
                 s=45, facecolors="none", edgecolors="#19766b", label="Meets both constraints", zorder=3)
score_ax.set(title="Score across the predefined grid", ylabel="Score (image change per minute)")
for dimension, minimum, color in [("H1", q1, "#19766b"), ("H2", q2, "#b06438")]:
    count_ax.plot(grid.lower, grid[dimension], "o-", ms=4, color=color, label=dimension)
    count_ax.axhline(minimum, color=color, ls=":", label=f"{dimension} constraint: {minimum:g}")
count_ax.set(title="Positive-lifetime feature counts", ylabel="Number of features")
for ax in axes[1]:
    if selected is not None:
        ax.axvline(selected.graph["hons"]["lower"], color="#ac364f", ls="--", label="Selected cutoff")
    ax.set(xlabel="Minimum total contact duration (minutes)", xlim=(-.1, 5.1))
    ax.grid(alpha=.2)
    ax.legend(fontsize=8)
seconds = None if selected is None else round(selected.graph["hons"]["lower"] * 60)
result = "No feasible cutoff" if seconds is None else f"Selected minimum duration: {seconds // 60} min {seconds % 60:02d} s"
fig.suptitle(f"Workplace contacts: raw records to HONS\n{result}", fontsize=15)
fig.text(.5, .025, "All participants remain; pale points are isolates. The 0-minute cell supplies a neighbor and is not selectable.",
         ha="center", fontsize=9)
fig.tight_layout(rect=(0, .05, 1, .94), h_pad=2)
fig.savefig(output / "contact_demo.png", dpi=180)
plt.close(fig)
print(f"{len(raw)} raw rows; {len(network)} participants; {network.number_of_edges()} pairs; window={window}")
print(f"HONS runtime: {elapsed:.3f} s; count constraints: H1 >= {q1:g}, H2 >= {q2:g}")
print("Selection:", None if selected is None else selected.graph["hons"])
print(f"Saved {output / 'contact_grid.csv'} and {output / 'contact_demo.png'}")
