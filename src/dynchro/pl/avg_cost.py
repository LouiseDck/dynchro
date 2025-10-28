from anndata import AnnData
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import numpy as np
from typing import Collection
from moepy import lowess


def plot_avg_cost_path(
        dataset: AnnData,
        pseudotime_key: str = "pseudotime",
        dtw_key: str = "dtw",
        trendline: bool = True,
        color = "#ffc000",
        ax: Axes | None = None
) -> Axes:
    if ax is None:
        fig, ax = plt.subplots()

    average_cost = dataset.obs.sort_values(pseudotime_key)[f"{dtw_key}_avg_cost"].values
    xvals = dataset.obs.sort_values(pseudotime_key)[pseudotime_key].values



    if trendline:
        # Fit and plot a LOWESS trend line

        quant_reg_func = lowess.calc_quant_reg_betas
        quantile_model = lowess.Lowess(reg_func=quant_reg_func)

        quantile_model.fit(xvals, average_cost, frac=0.3)
        x_pred = np.linspace(xvals.min(), xvals.max(), 100)
        y_pred = quantile_model.predict(x_pred)

        df_quantiles = lowess.quantile_model(xvals, average_cost, frac=0.3, num_fits=100)

        ax.plot(x_pred, y_pred, color='blue', linestyle='-', label='Trend Line (LOWESS)')
        ax.fill_between(df_quantiles.index, df_quantiles[0.1], df_quantiles[0.9], color='r', edgecolor='k', alpha=0.25, label='10-90% Prediction Interval')

        b, a = np.polyfit(xvals[:-1], average_cost[:-1], 1)
        ax.plot(xvals, b * xvals + a, color='black', linestyle='--', label='Trend Line (Cost)')

    # Plot the average cost path
    ax.scatter(xvals, average_cost, label="Average Cost Path", color = color, alpha=1, edgecolors='k')

    ax.set_xlabel("Pseudotime")
    ax.set_ylabel("Average Cost")
    ax.set_title("Average Cost on Path")
    ax.legend()

    return ax

