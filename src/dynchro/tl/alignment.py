from typing import List, Tuple

import numpy as np
import pandas as pd
import scipy.optimize
import scipy.spatial
from anndata import AnnData

from .dtw import dtw, traceback
from ..pl.matrices import plot_dtw_matrices

def get_matching_lineages(
    dataset1: AnnData,
    dataset2: AnnData,
    dtw_key: str = "dtw",
    pseudotime_key: str = "pseudotime",
    lineage_labels_key: str = "lineage_labels",
    plot: bool = False,
    mode: str = "only_results"
):
    common_vars = list(set(dataset1.var_names).intersection(set(dataset2.var_names)))

    lineages1 = [
        dataset1[dataset1.obs[label]][:, common_vars]
        for label in dataset1.uns[lineage_labels_key]
    ]

    lineages2 = [
        dataset2[dataset2.obs[label]][:, common_vars]
        for label in dataset2.uns[lineage_labels_key]
    ]

    # sort according to pseudotime_key
    lineages1 = [lin[lin.obs[pseudotime_key].argsort()] for lin in lineages1]
    lineages2 = [lin[lin.obs[pseudotime_key].argsort()] for lin in lineages2]

    distances = np.zeros((len(lineages1), len(lineages2)))

    for i, (label1, lineage1) in enumerate(zip(dataset1.uns[lineage_labels_key], lineages1)):
        for j, (label2, lineage2) in enumerate(zip(dataset2.uns[lineage_labels_key], lineages2)):
            if mode == "only_results":
                total_dist, cost, distances_matrix = dtw(
                    lineage1.X, lineage2.X, distance="euclidean", mode=mode
                )
                path1, path2 = traceback(D = distances_matrix)
                norm_distance = total_dist / (cost.shape[0] * cost.shape[1])# sum(len(x) for x in path1.values())
            elif mode == "copy":
                ll_dtw_key = f"{dtw_key}_{label1}_{label2}"

                lineage1, lineage2 = dtw(
                    lineage1, lineage2, distance="euclidean", mode=mode, dtw_key=ll_dtw_key
                )
                total_dist = lineage1.uns[f"{ll_dtw_key}_distance"]
                lineage1, lineage2 = traceback(lineage1, lineage2, dtw_key=ll_dtw_key, mode=mode, diag_mult=2.0)
                norm_distance = total_dist / (len(lineage1) * len(lineage2))#sum((len(x) for x in lineage1.obs[f"{ll_dtw_key}_path"].values))

            distances[i, j] = norm_distance

            if plot:
                plot_dtw_matrices(lineage1, lineage2, dtw_key=f"{dtw_key}_{label1}_{label2}")

    row_index, column_index = scipy.optimize.linear_sum_assignment(np.array(distances))

    return row_index, column_index, distances, lineages1, lineages2


def align_lineages(lineage1, lineage2, normalize=False, plot=False, mode = "only_results", dtw_key = "dtw"):
    if mode == "only_results":
        total_dist, cost, distances = dtw(lineage1, lineage2, distance="euclidean", mode="only_results")
    elif mode == "copy":
        lineage1, lineage2 = dtw(lineage1, lineage2, distance="euclidean", mode="copy", dtw_key=dtw_key)
        total_dist = lineage1.uns[f"{dtw_key}_distance"]
        cost = lineage1.obsm[f"{dtw_key}_cost"]
        distances = lineage1.obsm[f"{dtw_key}_D"]

    if normalize:
        total_dist = total_dist / (lineage1.shape[0] * lineage2.shape[0])

    if plot:
        print(dtw_key)

        lineage1, lineage2 = traceback(lineage1, lineage2, dtw_key=dtw_key, mode=mode)
        plot_dtw_matrices(lineage1, lineage2, dtw_key=dtw_key)

    return total_dist



# def get_matching_lineages(trajectories: List[AnnData], config=None) -> List[Tuple[str]]:
#     """Get the matching lineages between two trajectories.

#     Args:
#         trajectories (List[AnnData]): The two trajectories to compare, both AnnData objects.
#         config (Dict[str, str], optional): A configuration dictionary, defining when to use pseudocells. Defaults to None.

#     Returns:
#         List[Tuple[str]]: A List of matching lineages labels, where the first element of the tuple is the lineage label of the first trajectory, and the second element is the lineage label of the second trajectory.
#     """

#     # Only select the genes / variables that are present in all trajectories
#     common_vars = list(set.intersection(*(set(trajectory.var_names) for trajectory in trajectories)))

#     # Get the lineages for each trajectory
#     lineages = [
#         [
#             get_counts_common_vars(trajectory, config, "compare_trajectories_pseudocells", linlabel, common_vars)
#             for linlabel in trajectory.uns["lineage_labels"]
#         ]
#         for trajectory in trajectories
#     ]

#     # Calculate the dtw distance between all the lineages
#     distances = [
#         [align_lineages(lineage1, lineage2, normalize=False, traceback=False)[2] for lineage1 in lineages[0]]
#         for lineage2 in lineages[1]
#     ]

#     # Find the optimal matching between the lineages
#     row_index, column_index = scipy.optimize.linear_sum_assignment(np.array(distances))
#     matching = [
#         (trajectories[0].uns["lineage_labels"][i], trajectories[1].uns["lineage_labels"][i])
#         for i, j in zip(column_index, row_index)
#     ]

#     matching2 = [
#         (
#             trajectories[0],
#             trajectories[1],
#             trajectories[0].uns["lineage_labels"][i],
#             trajectories[1].uns["lineage_labels"][i],
#             f"{trajectories[0].uns['id']}_{trajectories[0].uns['lineage_labels'][i]}_{trajectories[1].uns['id']}_{trajectories[1].uns['lineage_labels'][i]}",
#         )
#         for i, j in zip(column_index, row_index)
#     ]

#     return row_index, column_index, distances  # matching, matching2


def align_trajectories(matching_lineages: List[Tuple[str]], pseudocells: bool = False) -> List[AnnData]:
    """Aligns two trajectories based on the matching lineages.

    Requires that the first two trajectories in each matching_lineages spot are the same anndata each time.

    Args:
        matching_lineages (List[Tuple[str]]): matching lineages, obtained using get_matching_lineages.
        trajectories (List[AnnData]): the trajectories to align.
        pseudocells (bool, optional): whether or not to use pseudocells. Defaults to False.

    Returns:
        List[AnnData]: The aligned trajectories.
    """
    for trajectory1, trajectory2, label1, label2, label in matching_lineages:
        common_vars = list(set(trajectory1.var_names).intersection(trajectory2.var_names))
        lineages = [
            get_common_vars(trajectory, label, common_vars, pseudocells=pseudocells, get_x=False)
            for trajectory, label in zip([trajectory1, trajectory2], [label1, label2])
        ]
        pseudotimes = [
            get_pseudotime(trajectory, label, pseudocells=pseudocells)
            for trajectory, label in zip([trajectory1, trajectory2], [label1, label2])
        ]

        if pseudocells:
            path1, path2, total_dist, cost, distances = align_lineages(lineages[0], lineages[1])
        else:
            lineages[0] = lineages[0][np.argsort(pseudotimes[0]), :]
            lineages[1] = lineages[1][np.argsort(pseudotimes[1]), :]
            path1, path2, total_dist, cost, distances = align_lineages(lineages[0].X, lineages[1].X)

        trajectory1.uns["id"]
        trajectory2.uns["id"]

        trajectory1.uns[f"alignment_costs_{label}"] = cost
        trajectory1.uns[f"alignment_path_{label}"] = path1
        trajectory1.uns[f"alignment_total_dist_{label}"] = total_dist
        trajectory1.uns[f"alignment_dist_{label}"] = distances
        trajectory2.uns[f"alignment_costs_{label}"] = cost
        trajectory2.uns[f"alignment_path_{label}"] = path2
        trajectory2.uns[f"alignment_total_dist_{label}"] = total_dist
        trajectory2.uns[f"alignment_dist_{label}"] = distances

        # TODO: save all this in the anndatas, or return this information somehow -> use an ID (alignment_id -> corrected pseudotimes in obs)
        # TODO: use alignment_id to store path & distance matrix in varm of uns
        index1, wpt1, index2, wpt2 = warp_pseudotime(
            lineages[0], lineages[1], path1, path2, pseudotimes[0], pseudotimes[1], pseudocells=pseudocells
        )

        add_warped_pseudotime(trajectory1, index1, wpt1, label)
        add_warped_pseudotime(trajectory2, index2, wpt2, label)

    trajectories = [matching_lineages[0][0], matching_lineages[0][1]]
    lineages = [
        [matching_lineages[i][2] for i in range(len(matching_lineages))],
        [matching_lineages[i][3] for i in range(len(matching_lineages))],
    ]
    labels = [matching_lineages[i][4] for i in range(len(matching_lineages))]

    merged = [merge_lineages(trajectory, lineages[i], labels=labels) for i, trajectory in enumerate(trajectories)]

    # to_merge = [[matching_lineages[i][j] for i in range(len(matching_lineages))] for j in range(2,4)]
    # to_merge = [(x, y[4]) for x, y in zip(to_merge, matching_lineages)]
    # # labels should be labels for the alignment path etc
    # # TODO: maybe merge the costs & paths
    # merged = [merge_lineages(trajectory, rest[0], labels = rest[1]) for trajectory, rest in zip([trajectory1, trajectory2], to_merge)]

    return merged


#######################################
# ACCESSORS FOR PSEUDOCELLS AND CELLS #
#######################################


def get_counts_common_vars(trajectory, config, step, lineage, common_vars, get_x=True):
    if config[step]:
        # pseudocells
        f"{step}_label"
        return trajectory[:, common_vars].varm[f"{config[f'{step}_label']}_{lineage}"].T.toarray().squeeze()

    else:
        if get_x:
            return trajectory[trajectory.obs[lineage]][:, common_vars].X
        else:
            return trajectory[trajectory.obs[lineage]][:, common_vars]


# Get certain vars from a lineage
def get_common_vars(trajectory1, label, vars, pseudocells=False, get_x=True):
    if pseudocells:
        label = f"pseudocells_100_{label}"
        return trajectory1[:, vars].varm[label].T.toarray().squeeze()
    else:
        if get_x:
            return trajectory1[trajectory1.obs[label]][:, vars].X
        else:
            return trajectory1[trajectory1.obs[label]][:, vars]


def get_pseudotime(adata, label, pseudocells=False):
    if pseudocells:
        label = f"pseudocells_100_{label}"
        return adata.uns[f"{label}_pseudotime"]
    else:
        return adata[adata.obs[label]].obs["pseudotime"]


# align two lineages using dtw
# def align_lineages(lineage1, lineage2, normalize=False, traceback=True):
#     total_dist, cost, distances = dtw(lineage1, lineage2, distance="euclidean")

#     if normalize:
#         total_dist = total_dist / (lineage1.shape[0] * lineage2.shape[0])

#     path1, path2 = [], []
#     if traceback:
#         path1, path2 = traceback(distances[1:, 1:])

#     return path1, path2, total_dist, cost, distances


def warp_pseudotime(cells1, cells2, path1, path2, pseudot1, pseudot2, pseudocells=False):
    changed_pt1 = {}
    changed_pt2 = {}
    for p1, p2 in zip(path1, path2):
        if pseudocells:
            pt1 = pseudot1[p1]
            pt2 = pseudot2[p2]
            # eeek iffy if a multiple pseudocells have the same pseudotime
            changed_pt1 = add_to_dict(changed_pt1, pt1, pt2)
            changed_pt2 = add_to_dict(changed_pt2, pt2, pt1)
        else:
            pt1 = pseudot1.iloc[p1]
            pt2 = pseudot2.iloc[p2]

            changed_pt1 = add_to_dict(changed_pt1, cells1.obs.iloc[p1].name, pt2)
            changed_pt2 = add_to_dict(changed_pt2, cells2.obs.iloc[p2].name, pt1)

    index1 = list(changed_pt1.keys())
    index2 = list(changed_pt2.keys())
    cpt1 = [sum(x) / len(x) for x in changed_pt1.values()]
    cpt2 = [sum(x) / len(x) for x in changed_pt2.values()]

    return index1, cpt1, index2, cpt2


def add_to_dict(d, key, value):
    if key in d:
        d[key].append(value)
        # nr, cur_avg = d[key]
        # cur_avg_update = (cur_avg + value) / (nr + 1)
        # d[key] = (nr+1, cur_avg_update)
    else:
        d[key] = [value]
    return d


# add warped pseudotime to lineage
def add_warped_pseudotime(lineage, new_index, warped_pt, label, pseudocells=False):
    if pseudocells:
        # TODO sort or sth? Or keep the pseudotime as a pandas series/dataframe in .uns?
        lineage.uns[label + "_warped_pseudotime"] = pd.Series(warped_pt, index=new_index)
    else:
        # lineage.obs.iloc[new_index, label + "_warped_pseudotime"] = warped_pt
        lineage.obs[f"{label}_warped_pseudotime"] = pd.NA
        lineage.obs.update(pd.DataFrame({f"{label}_warped_pseudotime": warped_pt}, index=new_index))
        lineage.obs[f"{label}_warped_pseudotime"] = pd.to_numeric(lineage.obs[f"{label}_warped_pseudotime"])

    return lineage


# merge two lineages of one single trajectory
def merge_lineages(trajectory, lineages, labels):
    common_cells_mask = trajectory.obs[lineages].all(axis=1)
    pseudotime_lineage_labels = [f"{label}_warped_pseudotime" for label in labels]
    common_label = f"{'_'.join(labels)}_warped_pseudotime"

    trajectory.obs[common_label] = pd.NA

    trajectory.obs.loc[common_cells_mask, common_label] = (
        trajectory[common_cells_mask].obs[pseudotime_lineage_labels].mean(axis=1)
    )

    trajectory.obs[common_label] = pd.to_numeric(trajectory.obs[common_label])

    # TODO merge the alignment paths as well

    return trajectory
