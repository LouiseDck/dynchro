import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns


def plot_dtw_lineage(t1, t2, mtx_label="dist", lineage_label1="lineage1", lineage_label2="lineage1"):
    # remove infinite values from the distance matrix
    fig = plt.figure()
    sns.heatmap(t1.uns[f"alignment_{mtx_label}_{lineage_label1}"][1:, 1:])

    return fig


# plot cost and dist matrices for each lineage in labels
def plot_dtw_trajectories(t1, t2, labels, mtx_labels=["costs", "dist"], fig=None, axes=None):
    # figure out minimum and maximum values for colorbars heatmaps
    # maybe the minimum value of costs should always be 0?
    min_val = min([t1.uns[f"alignment_{mtx_labels[0]}_{label}"][1:, 1:].min() for label in labels])
    max_val = max([t1.uns[f"alignment_{mtx_labels[0]}_{label}"][1:, 1:].max() for label in labels])
    min_val_dist = min([t1.uns[f"alignment_{mtx_labels[1]}_{label}"][1:, 1:].min() for label in labels])
    max_val_dist = max([t1.uns[f"alignment_{mtx_labels[1]}_{label}"][1:, 1:].max() for label in labels])

    # Plot both lineages
    if fig is None:
        fig, axes = plt.subplots(
            len(mtx_labels),
            len(labels) + 1,
            figsize=(len(labels) * 6, len(mtx_labels) * 5),
            gridspec_kw=dict(width_ratios=[1, 1, 0.05]),
        )
    for i, label in enumerate(labels):
        sns.heatmap(
            t1.uns[f"alignment_{mtx_labels[0]}_{label}"][1:, 1:], ax=axes[0, i], cbar=False, vmin=min_val, vmax=max_val
        ).set_title(label)
        sns.heatmap(
            t1.uns[f"alignment_{mtx_labels[1]}_{label}"][1:, 1:],
            ax=axes[1, i],
            cbar=False,
            vmin=min_val_dist,
            vmax=max_val_dist,
        ).set_title(label)

    fig.colorbar(axes[0][0].collections[0], cax=axes[0][2])
    fig.colorbar(axes[1][0].collections[0], cax=axes[1][2])

    return fig


def plot_dtw_clustermap(
    t1, t2, mtx_label="dist", lineage_label1="lineage1", lineage_label2="lineage1", name="clustermap.png"
):
    # remove infinite values from the distance matrix
    fig = plt.figure()
    rf = sns.clustermap(t1.uns[f"alignment_{mtx_label}_{lineage_label1}"][1:, 1:])
    rf.savefig(name)
    return fig


def plot_scatter(t1, label, offset=0, color="#3abbba", ax=None):
    if "X_umap" not in t1.obsm:
        sc.tl.pca(t1)
        sc.pp.neighbors(t1)
        sc.tl.umap(t1)

    sns.scatterplot(
        x=t1.obs[f"{label}_warped_pseudotime"],
        y=t1.obsm["X_umap"][:, 0] + offset,
        color=color,
        zorder=5,
        s=75,
        linewidth=1,
        edgecolor="black",
        ax=ax,
    )


def plot_scatter_trajectories(t1, labels, offsets=[0, 10, 20], colors=["#3abbba", "#ff681c", "#4B0082"], ax=None):
    if "X_umap" not in t1.obsm:
        sc.tl.pca(t1)
        sc.pp.neighbors(t1)
        sc.tl.umap(t1)

    for color, label, offset in zip(colors, labels, offsets):
        pseudotime = t1.obs[~t1.obs[f"{label}_warped_pseudotime"].isnull()][f"{label}_warped_pseudotime"]
        cells = pseudotime.index
        cells_integerindex = t1.obs.index.get_indexer(cells)
        ys = t1.obsm["X_umap"][cells_integerindex, 0]
        pseudotime = pseudotime.values.astype(float)

        sns.scatterplot(
            x=pseudotime,
            y=ys + offset,
            color=color,
            zorder=5,
            s=75,
            linewidth=1,
            edgecolor="black",
            ax=ax,
        )


def plot_dtw_trajectories_scatter(t1, t2, labels, mtx_labels=["costs", "dist"]):
    trajectories = [t1, t2]

    fig, axes = plt.subplots(
        len(labels) + 1,
        len(mtx_labels) + 1,
        figsize=(len(labels) * 6, (len(mtx_labels) + 1) * 5),
        gridspec_kw=dict(width_ratios=[1, 1, 0.05]),
    )

    fig = plot_dtw_trajectories(t1, t2, labels, mtx_labels, fig, axes[1:,])

    labels += ["dataset1_lineage1_dataset2_lineage1_dataset1_lineage2_dataset2_lineage2"]
    plot_scatter_trajectories(trajectories[0], labels, ax=axes[0, 0])
    plot_scatter_trajectories(trajectories[1], labels, ax=axes[0, 1])

    return fig


def plot_scatter_lines(t1, t2):
    sns.set(rc={"figure.figsize": (13, 7.5)})
    sns.set_style(style="white")

    plot_scatter(t1)
    plot_scatter(t2, offset=-5, color="#ff681c")


def calculate_coords(t1, t2):
    # is this lineage based?
    pass
