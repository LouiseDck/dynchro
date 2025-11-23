import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats import wasserstein_distance
from warnings import warn

from anndata import AnnData

def get_kde_eval(
    dataset: AnnData,
    pseudotime_key: str = "pseudotime",
    n_points: int = 100,
    bandwidth: float = 0.1,
    mode: str = "copy"
) -> tuple[np.ndarray, np.ndarray] | AnnData | None:
    vector = dataset.obs[pseudotime_key].values

    # check that vector is normalized between 0 and 1
    if np.min(vector) < 0 or np.max(vector) > 1:
        warn(f"Vector {pseudotime_key} is not normalized between 0 and 1. Normalizing it now and storing as 'norm_{pseudotime_key}'.")
        vector = (vector - np.min(vector)) / (np.max(vector) - np.min(vector))
        dataset.obs[f"norm_{pseudotime_key}"] = vector

    kde = gaussian_kde(vector, bw_method=bandwidth)
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
