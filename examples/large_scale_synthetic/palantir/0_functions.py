import anndata as ad
import numpy as np
import pandas as pd

import dynchro


def write_with_pc(adata, name):
    for pckey in adata.uns["pseudocells"]:
        if isinstance(adata.uns[f"{pckey}_pseudotime"], pd.core.series.Series):
            pseudocells = adata.uns[f"{pckey}_pseudotime"]
            pseudocells_index = pseudocells.index.values
            pseudocells_values = pseudocells.values

            adata.uns[f"{pckey}_pseudotime_index"] = pseudocells_index
            adata.uns[f"{pckey}_pseudotime_values"] = pseudocells_values

            adata.uns[f"{pckey}_pseudotime"] = "removed"

    adata.write_h5ad(f"{name}.h5ad")


def read_with_pc(file):
    adata = ad.read_h5ad(file)

    for pckey in adata.uns["pseudocells"]:
        if f"{pckey}_pseudotime" in adata.uns and adata.uns[f"{pckey}_pseudotime"] == "removed":
            pseudocells_index = adata.uns[f"{pckey}_pseudotime_index"]
            pseudocells_values = adata.uns[f"{pckey}_pseudotime_values"]

            pseudocells = pd.Series(pseudocells_values, index=pseudocells_index)
            adata.uns[f"{pckey}_pseudotime"] = pseudocells

    return adata


def get_config(trajectories, sort=True):
    config = {
        "compare_lineages_pseudocells": False,
        "compare_lineages_pseudocells_label": "pseudocells_50",
        "compare_trajectories_pseudocells": False,
        "compare_trajectories_pseudocells_label": "pseudocells_50",
        "align_pseudocells": False,
    }

    common_vars = list(set.intersection(*(set(trajectory.var_names) for trajectory in trajectories)))

    lineages = [
        [
            dynchro.tl.get_counts_common_vars(
                trajectory, config, "compare_trajectories_pseudocells", linlabel, common_vars, get_x=False
            )
            for linlabel in trajectory.uns["lineage_labels"]
        ]
        for trajectory in trajectories
    ]

    if sort:
        lineages = [[lin[np.argsort(lin.obs.pseudotime)] for lin in lineages1] for lineages1 in lineages]

    return lineages


import matplotlib.pyplot as plt
import seaborn as sns


def plot_cost_distances(cost, distances, path1, path2):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.set_size_inches(12, 4)
    sns.heatmap(cost, ax=ax1)
    distances = distances[1:, 1:]
    sns.heatmap(distances, ax=ax2)
    plt.scatter(path2, path1, c="white", s=3)
    plt.show()


class DTWResult:
    def __init__(self, lin1, lin2, cost, distances, path1, path2, dist, dist_norm, dist_norm_path):
        self._lin1 = lin1
        self._lin2 = lin2
        self._cost = cost
        self._distances = distances
        self._path1 = path1
        self._path2 = path2
        self._dist = dist
        self._dist_norm = dist_norm
        self._dist_norm_path = dist_norm_path

    def cost(self):
        return self._cost

    def distances(self):
        return self._distances

    def path1(self):
        return self._path1

    def path2(self):
        return self._path2

    def dist(self):
        return self._dist

    def dist_norm(self):
        return self._dist_norm

    def dist_norm_path(self):
        return self._dist_norm_path

    def get_warping(self):
        return self._path1, self._path2

    def plot_cost_distances(self):
        plot_cost_distances(self._cost, self._distances, self._path1, self._path2)


def dynchro_wrapper(ad1, ad2, distance_measure="correlation"):
    trajectories = [ad1, ad2]
    lineages = get_config(trajectories)

    results = []

    for lineage1, label1 in zip(lineages[0], ad1.uns["lineage_labels"]):
        for lineage2, label2 in zip(lineages[1], ad2.uns["lineage_labels"]):
            # print(lineage2.X.shape, lineage1.X.shape)
            lin1, lin2 = lineage1.X, lineage2.X
            total_dist, cost, distances = dynchro.tl.skip_dtw(lin1, lin2, distance=distance_measure)
            path1, path2 = dynchro.tl.skip_traceback(distances)

            total_dist_norm = total_dist / (lineage1.X.shape[0] + lineage2.X.shape[0])
            total_dist_norm2 = total_dist / len(path1)

            res = DTWResult(
                label1, label2, cost, distances, path1, path2, total_dist, total_dist_norm, total_dist_norm2
            )
            results.append(res)

            # plot_cost_distances(cost, distances, path1, path2)

    return results
