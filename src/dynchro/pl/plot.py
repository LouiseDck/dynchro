import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns


def plot_dtw(t1, dataset_label="dataset2"):
    print("Implement a plotting function here.")

    # Get all lineages
    # Print each lineage separately -> dtw distance matrix & then

    for lineage in t1.uns["lineage_labels"]:
        sns.heatmap(t1.uns[f"alignment_costs_{lineage}_{dataset_label}"])
        plt.show()

    return 0


def plot_scatter(t1, offset=0, color="#3abbba"):
    if "X_umap" not in t1.obsm:
        sc.tl.pca(t1)
        sc.pp.neighbors(t1)
        sc.tl.umap(t1)

    sns.scatterplot(
        x=t1.obs["lineage1_lineage2_warped_pseudotime"],
        y=t1.obsm["X_umap"][:, 0] + offset,
        color=color,
        zorder=5,
        s=75,
        linewidth=1,
        edgecolor="black",
    )


def plot_scatter_lines(t1, t2):
    sns.set(rc={"figure.figsize": (13, 7.5)})
    sns.set_style(style="white")

    plot_scatter(t1)
    plot_scatter(t2, offset=-5, color="#ff681c")


def calculate_coords(t1, t2):
    # is this lineage based?
    pass
