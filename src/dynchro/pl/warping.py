import matplotlib.pyplot as plt
from typing import Collection
import numpy as np
from anndata import AnnData
import pandas as pd

def norm(sequence: np.ndarray) -> np.ndarray:
    """
    Normalize a sequence so that each element lies between 0 and 1.
    """
    return (sequence - np.min(sequence)) / (np.max(sequence) - np.min(sequence))

def flatten(values: list) -> list:
    return [item for sublist in values for item in sublist]

def plot_warping(
        warped_datasets: Collection[AnnData],
        extra_datasets: Collection[AnnData] | None = None,
        pseudotime_key: str = "pseudotime",
        dtw_key: str = "dtw",
        dimred_key: str = "X_umap",
        ax = None,
):
# fix legends as well
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))

    # check only two warped datasets are provided
    assert  len(warped_datasets) <= 2, "Only two warped datasets are supported for plotting."

    colors = ["#ff681c", "#1b2944", "#3abbba", "#315b25"]

    for dataset, color in zip(warped_datasets + extra_datasets if extra_datasets else [], colors):
        dimred = dataset.obsm[dimred_key]
        pseudotime = norm(dataset.obs[pseudotime_key])
        ax.scatter(pseudotime, dimred[:, 1], c=color, edgecolors="black", linewidths=1.1, s=60)

    warping1 = flatten(warped_datasets[0].obs[f"{dtw_key}_path"])
    warping2 = flatten(warped_datasets[1].obs[f"{dtw_key}_path"])

    ds1_x = norm(warped_datasets[0].obs[pseudotime_key])[warping2]
    ds1_y = warped_datasets[0].obsm[dimred_key][warping2, 1]
    ds2_x = norm(warped_datasets[1].obs[pseudotime_key])[warping1]
    ds2_y = warped_datasets[1].obsm[dimred_key][warping1, 1]

    for x1, y1, x2, y2 in zip(ds1_x, ds1_y, ds2_x, ds2_y):
        ax.plot([x1, x2], [y1, y2], color='black', alpha = 0.3)

    ax.legend()

    ax.spines[['top', 'right', 'bottom', 'left']].set_visible(False)
    x0 = ax.get_xlim()[0]
    y0 = ax.get_ylim()[0]

    ax.arrow(x0, y0 - 2, 0.2, 0, head_width=0.3, head_length=0.01, fc='k', ec='k')

    ax.set_xlabel("Pseudotime")
    ax.xaxis.set_label_coords(0.11, 0.02)

    ax.set_xticks([])
    ax.set_yticks([])

    return ax
