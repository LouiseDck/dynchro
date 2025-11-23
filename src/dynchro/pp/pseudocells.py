from typing import List

import numpy as np
import pandas as pd
import scipy.stats
from anndata._core.anndata import AnnData

from .kernels import apply_gaussian_kernel, get_gaussian_transition_matrix


# Main function to calculate the pseudocells.
# Issue -> cannot put it in layers, might need to be a separate object as the cells themselves change?
# Might be put into varm but that is hacky as well
# Or into .uns as a separate anndata
# varp is not an option as it requires (var, var) as dimensions
#
# Current solution -> put it into .varm
# Should actually do mudata :/
def calculate_pseudocells(adata: AnnData, amount: int, lineage_label: str, bw: float = 0.05) -> AnnData:
    """Main function used to calculate pseudocells. It will generate the amount of pseudocells specified by the amount parameter for the specified lineage in the anndata object.
    The pseudocells will be stored in adata.varm[f"pseudocells_{amount}_{lineage_label}"] and the pseudotime of the pseudocells will be stored in adata.uns[f"pseudocells_{amount}_{lineage_label}_pseudotime"].

    Args:
        adata (AnnData): The anndata object.
        amount (int): The amount of pseudocells to generate for this lineage.
        lineage_label (str): The lineage for which to generate pseudocells.

    Returns:
        AnnData: The anndata object, updated with the pseudocells and pseudotime, stored in adata.varm[f"pseudocells_{amount}_{lineage_label}"] and adata.uns[f"pseudocells_{amount}_{lineage_label}_pseudotime"] respectively.
    """
    trunc_anndata = adata[adata.obs[lineage_label]]
    pseudocell_pseudotimes = interpolate_uniform(trunc_anndata, amount)
    trunc_anndata, pseudocells = smooth_pseudocells(trunc_anndata, pseudocell_pseudotimes, amount, bw)

    pseudocell_names = [f"pseudocell_{lineage_label}_{i}" for i in range(amount)]
    pseudocell_pseudotimes = {k: v for k, v in zip(pseudocell_names, pseudocell_pseudotimes)}

    # Save information into anndata
    trunc_anndata.varm[f"pseudocells_{amount}_{lineage_label}"] = pseudocells.T
    trunc_anndata.uns[f"pseudocells_{amount}_{lineage_label}_pseudotime"] = pd.Series(pseudocell_pseudotimes)
    if "pseudocells" not in trunc_anndata.uns:
        trunc_anndata.uns["pseudocells"] = []
    trunc_anndata.uns["pseudocells"] += [f"pseudocells_{amount}_{lineage_label}"]
    return trunc_anndata


def interpolate_uniform(adata, amount: int) -> List[int]:
    # interpolate pseudotime uniformly between min and max pseudotime
    min_pt = adata.obs.pseudotime.min()
    max_pt = adata.obs.pseudotime.max()
    return list(np.linspace(min_pt, max_pt, amount))


def interpolate_pseudocells(adata: AnnData, amount: int) -> List[int]:
    # interpolate pseudotime based on the distribution of real cells in the pseudotime
    return np.percentile(sorted(adata.obs.pseudotime.values), np.linspace(0.0, 100, amount))


def smooth_pseudocells(adata: AnnData, pseudocell_pseudotimes: List[int], amount: int, bw: float) -> np.ndarray:
    distances = calculate_pseudotime_distance(pseudocell_pseudotimes, sorted(adata.obs.pseudotime.values))
    adata.uns[f"transition_matrix_{amount}"] = np.array(
        [get_gaussian_transition_matrix(distances[i,]) for i in range(distances.shape[0])]
    )
    return adata, np.array([apply_gaussian_kernel(distances[i,], adata.X, bw) for i in range(distances.shape[0])])


def calculate_pseudotime_distance(pseudocell_pseudotimes: List[int], trajectory_pseudotime: np.ndarray) -> np.ndarray:
    distances = []
    for cell in pseudocell_pseudotimes:
        distance = trajectory_pseudotime - cell
        distances.append(distance)
    return np.array(distances)


def map_back(adata: AnnData, pseudocell_adata: AnnData, transition_matrix, warped_pt_label: str, treshold: float = 0.001) -> AnnData:
    mapped_warped_pseudotimes = []
    for i in range(transition_matrix.shape[1]):
        res = transition_matrix[:, i] / sum(transition_matrix[:, i])

        warped_pseudotimes = pseudocell_adata.obs[res > treshold][warped_pt_label].values
        contribution_amounts = res[res > treshold]

        result = np.average(warped_pseudotimes, weights=contribution_amounts)
        mapped_warped_pseudotimes.append(result)

    adata.obs[f"mapped_{warped_pt_label}"] = mapped_warped_pseudotimes
    return adata
