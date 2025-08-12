import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats import wasserstein_distance
from warnings import warn

from anndata import AnnData

def get_kde_eval(
    dataset: AnnData,
    pseudotime_key: str = "pseudotime",
    n_points: int = 100,
    mode: str = "copy"
) -> tuple[np.ndarray, np.ndarray] | AnnData | None:
    vector = dataset.obs[pseudotime_key].values

    # check that vector is normalized between 0 and 1
    if np.min(vector) < 0 or np.max(vector) > 1:
        warn(f"Vector {pseudotime_key} is not normalized between 0 and 1. Normalizing it now and storing as 'norm_{pseudotime_key}'.")
        vector = (vector - np.min(vector)) / (np.max(vector) - np.min(vector))
        dataset.obs[f"norm_{pseudotime_key}"] = vector

    kde = scipy.stats.gaussian_kde(vector)
    x_values = np.linspace(0, 1, n_points)
    y_values = kde(x_values)

    if mode == "only_results":
        return x_values, y_values

    dataset.uns[f"{pseudotime_key}_kde"] = {
        "x": x_values,
        "y": y_values
    }

    if mode == "copy":
        return dataset

    return None

# def get_wasserstein_distance(
#     adata : AnnData | None = None,
#     pseudotime_keys: tuple[str, str] | None = None,
#     distr1: tuple[np.ndarray, np.ndarray] | None = None,
#     distr2: tuple[np.ndarray, np.ndarray] | None = None
# ) -> float:

#     if adata is not None and pseudotime_keys is not None:
#         if len(pseudotime_keys) != 2:
#             raise ValueError("pseudotime_keys must be a tuple of two strings.")
#         if pseudotime_keys[0] not in adata.uns or pseudotime_keys[1] not in adata.uns:
#             raise ValueError(f"pseudotime_keys {pseudotime_keys} not found in adata.uns.")

#         kde_res1 = adata.uns[pseudotime_keys[0]]
#         kde_res2 = adata.uns[pseudotime_keys[1]]

#         distr1 = (kde_res1["x"], kde_res1["y"])
#         distr2 = (kde_res2["x"], kde_res2["y"])

#     return wasserstein_distance(distr1[0], distr2[0], distr1[1], distr2[1])


def get_wasserstein_distance(
    adata : AnnData | None = None,
    pseudotime_keys: tuple[str, str] | None = None,
    distr1: tuple[np.ndarray, np.ndarray] | None = None,
    distr2: tuple[np.ndarray, np.ndarray] | None = None,
    cutoff: float | None = None
) -> float:

    if adata is not None and pseudotime_keys is not None:
        if len(pseudotime_keys) != 2:
            raise ValueError("pseudotime_keys must be a tuple of two strings.")
        if pseudotime_keys[0] not in adata.uns or pseudotime_keys[1] not in adata.uns:
            raise ValueError(f"pseudotime_keys {pseudotime_keys} not found in adata.uns.")

        kde_res1 = adata.uns[pseudotime_keys[0]]
        kde_res2 = adata.uns[pseudotime_keys[1]]

        distr1 = (kde_res1["x"], kde_res1["y"])
        distr2 = (kde_res2["x"], kde_res2["y"])

    if cutoff is not None:
        distr1_before = get_cutoff_kde(distr1[0], distr1[1], (0.0, cutoff))
        distr2_before = get_cutoff_kde(distr2[0], distr2[1], (0.0, cutoff))

        before = wasserstein_distance(distr1_before[0], distr2_before[0], distr1_before[1], distr2_before[1])

        distr1_after = get_cutoff_kde(distr1[0], distr1[1], (cutoff, 1.0))
        distr2_after = get_cutoff_kde(distr2[0], distr2[1], (cutoff, 1.0))

        after = wasserstein_distance(distr1_after[0], distr2_after[0], distr1_after[1], distr2_after[1])

        return before, after

    else:
        return wasserstein_distance(distr1[0], distr2[0], distr1[1], distr2[1])


def get_cutoff_kde(
    xs: np.ndarray,
    ys: np.ndarray,
    between: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    """
    Get the cutoff point for the KDE plot.

    Parameters
    ----------
    xs : np.ndarray
        The x values of the KDE.
    ys : np.ndarray
        The y values of the KDE.
    cutoff : float
        The cutoff value to apply on the x values.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        The x and y values of the KDE after applying the cutoff.
    """
    mask = (xs >= between[0]) & (xs <= between[1])

    return xs[mask], ys[mask]
