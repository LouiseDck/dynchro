from collections.abc import Iterable
from anndata import AnnData
import matplotlib.pyplot as plt

def plot_kde_comparison(
    dataset: AnnData,
    pseudotime_keys: Iterable[str] = "pseudotime"
):

    colors = ["black", "blue", "red", "green", "orange", "purple"]
    kde_vals = [dataset.uns.get(f"{pseudotime_key}_kde") for pseudotime_key in pseudotime_keys]
    # x_values, y_values = get_kde_eval(dataset, pseudotime_key, bandwidth, n_points)
    fig = plt.figure(figsize=(6, 3))
    ax = fig.add_subplot(111)
    for kde_val, color in zip(kde_vals, colors):
        ax.plot(kde_val["x"], kde_val["y"], label='KDE', color=color)
    ax.set_title(f'KDE of {pseudotime_keys}')

    plt.xlabel('Pseudotime')
    plt.ylabel('Density')
    plt.legend()
    plt.show()
