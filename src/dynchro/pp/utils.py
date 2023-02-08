import numpy as np
import scanpy as sc


def find_path_in_tree(tree, paths, start, stops):
    fr = start
    tos = tree[fr]
    paths = [start]

    for to in tos:
        path = [to]
        while to in tree:
            to = tree[to][0]
            path.append(to)

        paths.append(path)

    return paths


def find_paths_in_tree(tree, starts, stops):
    return [find_path_in_tree(tree, [], start, stops) for start in starts]


def label_lineages(anndata, tree, starts, stops):
    # makes a dict from the pandas indexer with "from" as keys, the rest as a list of values
    # TODO: probably does not work yet with multiple tos from a single from
    tree = tree.set_index("from").T.to_dict("list")
    for key, value in tree.items():
        tree[key] = [value[0]]
    find_paths_in_tree(tree, starts, stops)


def label_lineage(anndata, clusterlabels, lineage_clusters, label):
    anndata.obs[label] = np.where(anndata.obs[clusterlabels].isin(lineage_clusters), True, False)
    if "lineage_labels" not in anndata.uns:
        anndata.uns["lineage_labels"] = []

    anndata.uns["lineage_labels"].append(label)

    return anndata


def preprocess_ti(adata, root=None, to_remove=None):
    if root is not None:
        adata.uns["iroot"] = np.flatnonzero(adata.obs["milestones"] == root)[0]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.scale(adata, max_value=10)
    sc.pp.neighbors(adata, random_state=1)
    sc.tl.umap(adata, random_state=1)
    sc.tl.paga(adata, groups="milestones")

    if to_remove is not None:
        adata = adata[adata.obs.leiden != to_remove]
    sc.tl.dpt(adata)

    adata.obs.rename(columns={"dpt_pseudotime": "pseudotime"}, inplace=True)
    # sc.pl.umap(adata, color = ["leiden", "dpt_pseudotime", "milestones"])
    return adata


def preprocess_batch_effect(dataset, remove=True):
    if remove:
        dataset = dataset[dataset.obs.milestones != "sEndC"]
    return dataset
