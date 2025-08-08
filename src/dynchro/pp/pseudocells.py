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
def calculate_pseudocells(adata: AnnData, amount: int, lineage_label: str) -> AnnData:
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
    pseudocell_pseudotimes = interpolate_pseudocells(trunc_anndata, amount)
    trunc_anndata, pseudocells = smooth_pseudocells(trunc_anndata, pseudocell_pseudotimes, amount)

    pseudocell_names = [f"pseudocell_{lineage_label}_{i}" for i in range(amount)]
    pseudocell_pseudotimes = {k: v for k, v in zip(pseudocell_names, pseudocell_pseudotimes)}

    # Save information into anndata
    adata.varm[f"pseudocells_{amount}_{lineage_label}"] = pseudocells.T
    adata.uns[f"pseudocells_{amount}_{lineage_label}_pseudotime"] = pd.Series(pseudocell_pseudotimes)
    if "pseudocells" not in adata.uns:
        adata.uns["pseudocells"] = []
    adata.uns["pseudocells"] += [f"pseudocells_{amount}_{lineage_label}"]
    return adata


def interpolate_pseudocells(adata: AnnData, amount: int) -> List[int]:
    # interpolate pseudotime based on the distribution of real cells in the pseudotime
    return np.percentile(sorted(adata.obs.pseudotime.values), np.linspace(0.0, 100, amount))


def smooth_pseudocells(adata: AnnData, pseudocell_pseudotimes: List[int], amount: int) -> np.ndarray:
    distances = calculate_pseudotime_distance(pseudocell_pseudotimes, sorted(adata.obs.pseudotime.values))
    adata.uns[f"transition_matrix_{amount}"] = np.array(
        [get_gaussian_transition_matrix(distances[i,], adata.X) for i in range(distances.shape[0])]
    )
    return adata, np.array([apply_gaussian_kernel(distances[i,], adata.X) for i in range(distances.shape[0])])


def calculate_pseudotime_distance(pseudocell_pseudotimes: List[int], trajectory_pseudotime: np.ndarray) -> np.ndarray:
    distances = []
    for cell in pseudocell_pseudotimes:
        distance = trajectory_pseudotime - cell
        distances.append(distance)
    return np.array(distances)


# Merges pseudocells that appear in multiple lineages, so that pseuodcells lie minimum 1 percentile apart.
# Might need to be tested for common pseudocells in:
# - 3 lineages: DONE & works
# - at the end of a lineage
# - in the middle of a lineage
# TODO: let it work for arbitrary overlapping sections between the lineages?
def merge_pseudocells_lineages(adata: AnnData, lineages: str, pseudocell_amount: int) -> AnnData:
    """Merges pseudocells from the different lineages, so that there are less overlapping pseudocells.
    Tries to ensure that the individual lineage dynamics stay preserved.
    Important: will only work if there is 1 overlapping section between the lineages.

    Args:
        adata (AnnData): The anndata object containing the pseudocells.
        lineages (str): The lineages for which to merge the pseudocells.
        pseudocell_amount (int): The amount of pseudocells that were generated for each lineage.

    Returns:
        AnnData: An updated anndata object, containing the merged pseudocells.
    """
    common_cells = adata[adata.obs[lineages].all(axis=1)]
    minpt = min(common_cells.obs.pseudotime)
    maxpt = max(common_cells.obs.pseudotime)

    cell_pseudotime = common_cells.obs.pseudotime

    # construct the labels of the lineages
    pseudocell_lineages = [
        adata.varm[f"pseudocells_{pseudocell_amount}_{lineage_label}"].T for lineage_label in lineages
    ]
    pseudocell_lineages_pseudotimes = [
        adata.uns[f"pseudocells_{pseudocell_amount}_{lineage_label}_pseudotime"]  # .to_numpy()
        for lineage_label in lineages
    ]

    # construct the masks in order to identify the common pseudocells, and use it to
    masks = [(plp >= minpt) & (plp <= maxpt) for plp in pseudocell_lineages_pseudotimes]

    common_pseudocells = np.concatenate([pll[np.where(mask)[0], :] for pll, mask in zip(pseudocell_lineages, masks)])
    pseudotimes = pd.concat([plp[mask] for plp, mask in zip(pseudocell_lineages_pseudotimes, masks)])
    # sort for the iterative merging & replacement
    pseudocells = common_pseudocells[np.argsort(pseudotimes), :]
    pseudotimes_sorted = pseudotimes[np.argsort(pseudotimes)]

    # Determine which pseudocells to keep of the min/max range
    rest_pseudocells, rest_pseudotimes = iterative_keep_pseudocells(pseudocells, cell_pseudotime, pseudotimes_sorted)

    # replace the pseudocells in the anndata object with the kept pseudocells, must be done for each lineage
    for pll, plp, lin in zip(pseudocell_lineages, pseudocell_lineages_pseudotimes, lineages):
        leftindex = np.searchsorted(plp, minpt)
        rightindex = np.searchsorted(plp, maxpt)

        adata.varm[f"pseudocells_{pseudocell_amount}_{lin}"] = np.concatenate(
            [pll[:leftindex, :], rest_pseudocells, pll[rightindex:, :]]
        ).T
        adata.uns[f"pseudocells_{pseudocell_amount}_{lin}_pseudotime"] = np.concatenate(
            [plp[:leftindex], rest_pseudotimes, plp[rightindex:]]
        )

    return adata


def iterative_keep_pseudocells(pseudocells, cell_pseudotime, pseudotimes_sorted):
    i = 0
    cur = pseudocells[i]
    nxt = pseudocells[i + 1]

    cur = i
    nxt = i + 1

    keeping = [True for _ in range(len(pseudocells))]

    while i < len(pseudocells) - 2:
        cur_score = scipy.stats.percentileofscore(cell_pseudotime, pseudotimes_sorted[cur])
        nxt_score = scipy.stats.percentileofscore(cell_pseudotime, pseudotimes_sorted[nxt])

        if nxt_score - cur_score < 1:
            keeping[i] = False
            nxt = i + 2
        else:
            cur = i + 1
            nxt = i + 2
        i += 1

    rest_pseudocells = pseudocells[keeping, :]
    rest_pseudotimes = pseudotimes_sorted[keeping]

    return rest_pseudocells, rest_pseudotimes
