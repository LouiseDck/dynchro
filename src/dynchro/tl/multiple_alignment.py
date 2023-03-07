import anndata as ad
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy
import scipy.spatial.distance
from sklearn.cluster import AgglomerativeClustering

from .alignment import (
    align_lineages,
    align_trajectories,
    get_counts_common_vars,
    get_matching_lineages,
)

#####################
# MULTI COMPARISONS #
#####################


# compare multiple trajectories & align them one after one
# in a hierarchical mannerk
# returns the normalized distances between all the (pseudocell) trajectories
# TOD) will now match with itself all the time
def multi_compare(trajectories, lineage_labels, config=None, pseudocell_labels=None):
    results = []
    indices = []
    for i, traj1 in enumerate(trajectories):
        results2 = []
        indices.append(traj1.uns["id"])
        for j in range(0, len(trajectories)):
            traj2 = trajectories[j]
            # This is not new right, this should just use normal matching & alignment
            cost = compare_trajectories(traj1, lineage_labels, traj2, lineage_labels, config=config)

            results2.append(cost)
            # results_euclidean2.append(cost2)
        results.append(results2)
        # results_euclidean.append(results_euclidean2)

    distances = pd.DataFrame(data=results, index=indices, columns=indices)

    linkage = scipy.cluster.hierarchy.linkage(scipy.spatial.distance.squareform(results), method="average")

    clustering2 = AgglomerativeClustering(
        n_clusters=3, affinity="precomputed", linkage="complete", compute_full_tree=True
    ).fit(results)

    import itertools

    ii = itertools.count(len(results))
    thing = [{"node_id": next(ii), "left": x[0], "right": x[1]} for x in clustering2.children_]

    # clustering_euclidean = sklearn.cluster.AgglomerativeClustering(n_clusters=3, affinity = 'precomputed', linkage='complete', compute_full_tree=True).fit(results_euclidean)
    # thing_euclidean = [{'node_id': next(ii), 'left': x[0], 'right':x[1]} for x in clustering_euclidean.children_]

    references = []
    for label in set(clustering2.labels_):
        cluster_trajectories = [trajectories[i] for i in np.where(clustering2.labels_ == label)[0]]
        aligned = iterative_align(cluster_trajectories, config=config)
        references.append(aligned)

    return references


# compute a distance matrix between all trajectories
def compare_trajectories(t1, llabels1, t2, llabels2, config=None):
    common_vars = list(set(t1.var_names).intersection(t2.var_names))

    lineages1 = [
        get_counts_common_vars(t1, config, "compare_trajectories_pseudocells", linlabel, common_vars)
        for linlabel in llabels1
    ]
    lineages2 = [
        get_counts_common_vars(t2, config, "compare_trajectories_pseudocells", linlabel, common_vars)
        for linlabel in llabels2
    ]
    # # split in lineages
    # lineages1 = [t1[t1.obs[label]] for label in llabels1]
    # lineages2 = [t2[t2.obs[label]] for label in llabels2]

    results = [
        [align_lineages(lineage1, lineage2, normalize=True)[2] for lineage1 in lineages1] for lineage2 in lineages2
    ]
    row_index, column_index = scipy.optimize.linear_sum_assignment(np.array(results))
    cost = np.sum(np.array(results)[row_index, column_index])

    return cost


def compare_euclidean(t1, llabels1, t2, llabels2, config=None):
    common_vars = list(set(t1.var_names).intersection(t2.var_names))

    lineages1 = [
        get_counts_common_vars(t1, config, "compare_trajectories_pseudocells", linlabel, common_vars)
        for linlabel in llabels1
    ]
    lineages2 = [
        get_counts_common_vars(t2, config, "compare_trajectories_pseudocells", linlabel, common_vars)
        for linlabel in llabels2
    ]

    results = [
        [euclidean_lineages(lineage1, lineage2, normalize=True) for lineage1 in lineages1] for lineage2 in lineages2
    ]
    row_index, column_index = scipy.optimize.linear_sum_assignment(np.array(results))
    cost = np.sum(np.array(results)[row_index, column_index])

    return cost


def euclidean_lineages(lineage1, lineage2, normalize=True):
    total_dist = np.sum(
        np.linalg.norm(lineage1 - lineage2, axis=-1)
    )  # scipy.spatial.distance.euclidean(lineage1, lineage2)

    if normalize:
        total_dist = total_dist / (lineage1.shape[0] * lineage2.shape[0])

    return total_dist


# iterative add trajectories to each other
# reference is now the first trajectory
def iterative_align(trajectories, config):
    orig = trajectories[0]
    orig_label = orig.uns["id"]
    orig.obs["orig"] = orig_label
    for trajectory in trajectories[1:]:
        # TODO: save the lineage labels in .uns

        matching_lineages = get_matching_lineages(
            orig, trajectory, ["lineage1", "lineage2"], ["lineage1", "lineage2"], config=config
        )
        align_trajectories(
            matching_lineages, [orig, trajectory], [["lineage1", "lineage2"], ["lineage1", "lineage2"]], plot=False
        )

        trajectory.obs["orig"] = trajectory.uns["id"]
        # do sth so that orig = merged trajectories
        orig = merge_two_trajectories(orig, trajectory)
    return orig


# Is this reasonable?
# todo, merge pseudocells?
def merge_two_trajectories(traj1, traj2):
    result = ad.concat([traj1, traj2], merge="same", uns_merge="first")
    varm1 = traj1.varm
    varm2 = traj2.varm

    new_varm = varm1

    for key, value in varm2.items():
        if key not in varm1:
            new_varm[key] = value
        else:
            new_varm[key] = np.concatenate((varm1[key], varm2[key]), axis=1)

    result.varm = new_varm

    return result


# TODO
def merge_uns(dicts):
    # get labels from uns

    return 0
