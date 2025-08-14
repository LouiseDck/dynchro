from anndata import AnnData
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from typing import Collection


def plot_avg_cost_path(
        dataset: AnnData,
        pseudotime_key: str = "pseudotime",
        dtw_key: str = "dtw",
        trendline: bool = True,
        ax: Axes | None = None
) -> Axes:
    if ax is None:
        fig, ax = plt.subplots()

    average_cost = dataset.obs.sort_values(pseudotime_key)[f"{dtw_key}_avg_cost"].values
    xvals = dataset.obs.sort_values(pseudotime_key)[pseudotime_key].values

    # Plot the average cost path
    ax.scatter(xvals, average_cost, label="Average Cost Path", color = "orange", alpha=0.5)

    b, a = np.polyfit(xvals[:-1], average_cost[:-1], 1)
    ax.plot(xvals, b * xvals + a, color='black', linestyle='--', label='Trend Line (Cost)')

    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Average Cost")
    ax.set_title("Average Cost on Path")
    ax.legend()

    return ax

