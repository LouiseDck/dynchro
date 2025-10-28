import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from anndata import AnnData

def plot_cost_matrix(
    cost: np.ndarray,
    reference: AnnData | None = None,
    query: AnnData | None = None,
    dtw_key: str = "dtw",
    plot_paths: bool = False,
    cmap = "viridis",
    path_color = "red",
    ax = None
) -> None:
    """
    Plot the cost matrix of the DTW algorithm.
    """
    if ax is None:
        ax = plt.gca()
    sns.heatmap(cost, cmap=cmap, ax=ax)
    ax.set_title(f"{dtw_key} Cost Matrix")
    if plot_paths and reference is not None and query is not None:
        path1 = flatten(reference.obs[f"{dtw_key}_path"])
        path2 = flatten(query.obs[f"{dtw_key}_path"])
        ax.plot(path1, path2, color=path_color, linewidth=1)
    ax.set_xlabel("Query")
    ax.set_ylabel("Reference")
    # plt.show()


def plot_distance_matrix(
    distance: np.ndarray,
    reference: AnnData | None = None,
    query: AnnData | None = None,
    dtw_key: str = "dtw",
    plot_paths: bool = True,
    ax = None
) -> None:
    """
    Plot the distance matrix of the DTW algorithm.
    """
    if ax is None:
        ax = plt.gca()
    sns.heatmap(distance, cmap="viridis", ax=ax)
    ax.set_title(f"{dtw_key} Distance Matrix\nTotal distance: {reference.uns[f"{dtw_key}_distance"]:.2f}")
    if plot_paths and reference is not None and query is not None:
        path1 = flatten(reference.obs[f"{dtw_key}_path"])
        path2 = flatten(query.obs[f"{dtw_key}_path"])
        ax.plot(path1, path2, color='red', linewidth=1)
    ax.set_xlabel("Query")
    ax.set_ylabel("Reference")
    # plt.show()

def flatten(values : list) -> list:
    """
    Flatten a list of lists into a single list.
    """
    flat_list = []

    for sublist in values:
        if isinstance(sublist, list):
            # If the item is a list, extend the result with its contents
            flat_list.extend(sublist)
        else:
            # If the item is not a list, append it directly if it is not None
            if sublist is not None:
                flat_list.append(sublist)

    return flat_list

def plot_dtw_matrices(
    reference: AnnData,
    query: AnnData,
    dtw_key: str = "dtw",
    axes = None
) -> None:
    """
    Plot the DTW matrices of the DTW algorithm.
    """
    cost = reference.obsm[f"{dtw_key}_cost"]
    distance = reference.obsm[f"{dtw_key}_D"]

    if axes is None:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    plot_cost_matrix(cost, reference, query, dtw_key, ax=axes[0])
    plot_distance_matrix(distance, reference, query, dtw_key, ax=axes[1])
