import anndata
import matplotlib.pyplot as plt
import numpy as np
import pytest

import dynchro


@pytest.fixture(scope="session")
def dataset1():
    d1 = anndata.read("tests/data/batcheffect_dataseta0.h5ad")
    d1.uns["id"] = "dataset1"
    d1 = dynchro.pp.preprocess_ti(d1, root="sB")
    d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sC", "sEndC"], "lineage1")
    d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sD", "sEndD"], "lineage2")

    d1.X = d1.X.todense()

    # TODO: make this into one function
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage1")
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage2")
    d1 = dynchro.pp.merge_pseudocells_lineages(d1, ["lineage1", "lineage2"], 100)

    yield d1


@pytest.fixture(scope="session")
def dataset1_shuffled():
    d1 = anndata.read("tests/data/batcheffect_dataseta0.h5ad")
    # d1.X[cells, genes] -> only shuffle the genes
    # this only works when shuffling everything

    d1.uns["id"] = "dataset1"
    d1 = dynchro.pp.preprocess_ti(d1, root="sB")
    d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sC", "sEndC"], "lineage1")
    d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sD", "sEndD"], "lineage2")

    d1.X = d1.X.todense()
    np.random.shuffle(d1.X)

    # TODO: make this into one function
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage1")
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage2")
    d1 = dynchro.pp.merge_pseudocells_lineages(d1, ["lineage1", "lineage2"], 100)

    yield d1


@pytest.fixture(scope="session")
def dataset2_shuffled():
    d1 = anndata.read("tests/data/batcheffect99_datasetb0.h5ad")
    # d1.X[cells, genes] -> only shuffle the genes
    # this only works when shuffling everything

    d1.uns["id"] = "dataset2"
    d1 = dynchro.pp.preprocess_ti(d1, root="sB")
    d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sC", "sEndC"], "lineage1")
    d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sD", "sEndD"], "lineage2")

    d1.X = d1.X.todense()
    np.random.shuffle(d1.X)

    # TODO: make this into one function
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage1")
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage2")
    d1 = dynchro.pp.merge_pseudocells_lineages(d1, ["lineage1", "lineage2"], 100)

    yield d1


@pytest.fixture(scope="session")
def dataset2():
    d2 = anndata.read("tests/data/batcheffect99_datasetb0.h5ad")
    d2.uns["id"] = "dataset2"
    d2 = dynchro.pp.preprocess_ti(d2, root="sB")
    d2 = dynchro.pp.label_lineage(d2, "milestones", ["sA", "sB", "SBmid", "sC", "sEndC"], "lineage1")
    d2 = dynchro.pp.label_lineage(d2, "milestones", ["sA", "sB", "SBmid", "sD", "sEndD"], "lineage2")
    d2.X = d2.X.todense()

    d2 = dynchro.pp.calculate_pseudocells(d2, 100, "lineage1")
    d2 = dynchro.pp.calculate_pseudocells(d2, 100, "lineage2")
    d2 = dynchro.pp.merge_pseudocells_lineages(d2, ["lineage1", "lineage2"], 100)

    yield d2


@pytest.fixture(scope="session")
def dataset4():
    d2 = anndata.read("tests/data/batcheffect99_datasetb1.h5ad")
    d2.uns["id"] = "dataset2"
    d2 = dynchro.pp.preprocess_ti(d2, root="sB")
    d2 = dynchro.pp.label_lineage(d2, "milestones", ["sA", "sB", "SBmid", "sC", "sEndC"], "lineage1")
    d2 = dynchro.pp.label_lineage(d2, "milestones", ["sA", "sB", "SBmid", "sD", "sEndD"], "lineage2")
    d2.X = d2.X.todense()

    d2 = dynchro.pp.calculate_pseudocells(d2, 100, "lineage1")
    d2 = dynchro.pp.calculate_pseudocells(d2, 100, "lineage2")
    d2 = dynchro.pp.merge_pseudocells_lineages(d2, ["lineage1", "lineage2"], 100)

    yield d2


@pytest.fixture(scope="session")
def dataset_complex1():
    d1 = anndata.read("tests/data/tree_batcheffect75_dataseta1.h5ad")
    d1.uns["iroot"] = np.flatnonzero(d1[d1.uns["iroot"]])[0]
    d1 = dynchro.pp.preprocess_ti(d1, root="sB")
    lin1 = ["sA", "sB", "sC", "sCmid", "sE", "sEndE"]
    lin2 = ["sA", "sB", "sC", "sCmid", "sD", "sDmid", "sF", "sEndF"]
    lin3 = ["sA", "sB", "sC", "sCmid", "sD", "sDmid", "sG", "sEndG"]

    d1 = dynchro.pp.label_lineage(d1, "milestones", lin1, "lineage1")
    d1 = dynchro.pp.label_lineage(d1, "milestones", lin2, "lineage2")
    d1 = dynchro.pp.label_lineage(d1, "milestones", lin3, "lineage3")

    # TODO: make this into one function
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage1")
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage2")
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage3")
    d1 = dynchro.pp.merge_pseudocells_lineages(d1, ["lineage1", "lineage2", "lineage3"], 100)

    yield d1


@pytest.fixture(scope="session")
def dataset_complex2():
    d1 = anndata.read("tests/data/tree_batcheffect75_datasetb1.h5ad")
    d1.uns["iroot"] = np.flatnonzero(d1[d1.uns["iroot"]])[0]
    d1 = dynchro.pp.preprocess_ti(d1, root="sB")
    lin1 = ["sA", "sB", "sC", "sCmid", "sE", "sEndE"]
    lin4 = ["sA", "sB", "sC", "sCmid", "sD", "sDmid"]

    d1 = dynchro.pp.label_lineage(d1, "milestones", lin1, "lineage1")
    d1 = dynchro.pp.label_lineage(d1, "milestones", lin4, "lineage2")

    # TODO: make this into one function
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage1")
    d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage2")
    d1 = dynchro.pp.merge_pseudocells_lineages(d1, ["lineage1", "lineage2"], 100)

    yield d1


class TestAlignment:
    def test_matching_lineages(self, dataset1, dataset2):
        config = {
            "compare_lineages_pseudocells": False,
            "copmare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": False,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": False,
        }

        matching_lineages = dynchro.tl.get_matching_lineages([dataset1, dataset2], config=config)
        assert matching_lineages == [("lineage1", "lineage2"), ("lineage1", "lineage2")]

    def test_matching_pseudocells(self, dataset1, dataset2):
        config = {
            "compare_lineages_pseudocells": True,
            "copmare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": False,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": False,
        }

        matching_lineages = dynchro.tl.get_matching_lineages([dataset1, dataset2], config=config)
        assert matching_lineages == [("lineage1", "lineage2"), ("lineage1", "lineage2")]

    # @pytest.mark.skip(reason="Don't know how to test this yet, maybe plots?")
    def test_align_trajectories_simple(self, dataset1, dataset2):
        config = {
            "compare_lineages_pseudocells": False,
            "compare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": False,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": False,
        }

        matching_lineages, matching = dynchro.tl.get_matching_lineages(
            [dataset1, dataset2],
            config=config,
        )

        result = dynchro.tl.align_trajectories(
            matching,
            pseudocells=False,
        )
        dynchro.pl.plot_scatter_lines(dataset1, dataset2)
        plt.show()

        # TODO: how to test this?
        assert 1 == 0

    def test_align_trajectories_complex(self, dataset_complex1, dataset_complex2):
        config = {
            "compare_lineages_pseudocells": False,
            "copmare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": False,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": False,
        }

        matching_lineages = dynchro.tl.get_matching_lineages(
            [dataset_complex1, dataset_complex2],
            ["lineage1", "lineage2", "lineage3"],
            ["lineage1", "lineage2"],
            config=config,
        )
        dynchro.tl.align_trajectories(
            matching_lineages,
            [dataset_complex1, dataset_complex2],
            pseudocells=False,
        )

        # TODO: how to test this?
        assert 1 == 0

    def _calculate_correlations(self, cost_matrix, p1, p2):
        import scipy as sp

        pseud_dist = []
        for x in p1:
            pdt = []
            for y in p2:
                pdt.append(abs(x - y))
            pseud_dist.append(pdt)

        correlations = []
        for i in enumerate(cost_matrix):
            corr = sp.stats.pearsonr(i[1], pseud_dist[i[0]])[0]
            correlations.append(corr)
        return correlations

    def test_shuffled(self, dataset1_shuffled, dataset2_shuffled, name="test_shuffled.png"):
        config = {
            "compare_lineages_pseudocells": False,
            "compare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": False,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": False,
        }

        matching_lineages, matching = dynchro.tl.get_matching_lineages(
            [dataset1_shuffled, dataset2_shuffled],
            config=config,
        )

        result = dynchro.tl.align_trajectories(
            matching,
            pseudocells=False,
        )

        # get pseudotime differences
        # get differences in cost matrix
        cost_matrix = result[0].uns["alignment_costs_dataset1_lineage1_dataset2_lineage1"]
        pseudotime_d1_l1 = result[0][result[0].obs.lineage1, :].obs.pseudotime
        pseudotime_d2_l1 = result[1][result[1].obs.lineage1, :].obs.pseudotime

        correlations = self._calculate_correlations(cost_matrix, pseudotime_d1_l1, pseudotime_d2_l1)

        import matplotlib.pyplot as plt

        plt.scatter(range(len(correlations)), correlations)

        _ = 0
        # get correlation between pseudotime and cost matrix

        fig3 = dynchro.pl.plot_dtw_clustermap(
            result[0],
            result[1],
            mtx_label="costs",
            lineage_label1=matching[0][4],
            lineage_label2=matching[0][4],
            name="cmap_shuffled.png",
        )

        fig1 = dynchro.pl.plot_dtw_lineage(
            result[0], result[1], mtx_label="dist", lineage_label1=matching[0][4], lineage_label2=matching[0][4]
        )
        fig1.savefig(f"dist_{name}")

        fig2 = dynchro.pl.plot_dtw_lineage(
            result[0], result[1], mtx_label="costs", lineage_label1=matching[0][4], lineage_label2=matching[0][4]
        )
        fig2.savefig(f"cost_{name}")

    def test_non_shuffled(self, dataset1, dataset4, name="test_non_shuffled.png"):
        config = {
            "compare_lineages_pseudocells": False,
            "compare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": False,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": False,
        }

        matching_lineages, matching = dynchro.tl.get_matching_lineages(
            [dataset1, dataset4],
            config=config,
        )

        result = dynchro.tl.align_trajectories(
            matching,
            pseudocells=False,
        )

        cost_matrix = result[0].uns["alignment_costs_dataset1_lineage1_dataset2_lineage1"]
        pseudotime_d1_l1 = result[0][result[0].obs.lineage1, :].obs.pseudotime
        pseudotime_d2_l1 = result[1][result[1].obs.lineage1, :].obs.pseudotime
        correlations = self._calculate_correlations(cost_matrix, pseudotime_d1_l1, pseudotime_d2_l1)

        import matplotlib.pyplot as plt

        plt.scatter(range(len(correlations)), correlations)

        fig3 = dynchro.pl.plot_dtw_clustermap(
            result[0],
            result[1],
            mtx_label="costs",
            lineage_label1=matching[0][4],
            lineage_label2=matching[0][4],
            name="cmap_non_shuffled.png",
        )
        fig3.savefig(f"clmap_costs_{name}")

        fig05 = dynchro.pl.plot_dtw_trajectories_scatter(result[0], result[1], labels=[matching[0][4], matching[1][4]])
        fig05.savefig(f"trajectories_scatter_{name}")

        fig0 = dynchro.pl.plot_dtw_trajectories(result[0], result[1], labels=[matching[0][4], matching[1][4]])
        fig0.savefig(f"trajectories_{name}")

        fig1 = dynchro.pl.plot_dtw_lineage(
            result[0], result[1], mtx_label="dist", lineage_label1=matching[0][4], lineage_label2=matching[0][4]
        )
        fig1.savefig(f"dist_{name}")

        fig2 = dynchro.pl.plot_dtw_lineage(
            result[0], result[1], mtx_label="costs", lineage_label1=matching[0][4], lineage_label2=matching[0][4]
        )
        fig2.savefig(f"cost_{name}")

        _ = 0

    def test_non_shuffled_pseudocells(self, dataset1, dataset2, name="test_non_shuffled.png"):
        config = {
            "compare_lineages_pseudocells": True,
            "compare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": True,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": True,
        }

        matching_lineages, matching = dynchro.tl.get_matching_lineages(
            [dataset1, dataset2],
            config=config,
        )

        result = dynchro.tl.align_trajectories(
            matching,
            pseudocells=True,
        )

        fig = dynchro.pl.plot_dtw_lineage(
            result[0], result[1], lineage_label1=matching[0][4], lineage_label2=matching[0][4]
        )

        fig.savefig(name)

    def test_shuffled_non_shuffed(self, dataset1, dataset2, dataset1_shuffled, dataset2_shuffled):
        self.test_shuffled(dataset1_shuffled, dataset2_shuffled, name="shuffled.png")
        self.test_non_shuffled(dataset1, dataset2, name="normal")

        self.test_non_shuffled(dataset1, dataset1, name="same.png")
        self.test_shuffled_same()

    def test_shuffled_same(self):
        d1 = anndata.read("tests/data/batcheffect_dataseta0.h5ad")
        # d1.X[cells, genes] -> only shuffle the genes
        # this only works when shuffling everything

        d1.uns["id"] = "dataset1"
        d1 = dynchro.pp.preprocess_ti(d1, root="sB")
        d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sC", "sEndC"], "lineage1")
        d1 = dynchro.pp.label_lineage(d1, "milestones", ["sA", "sB", "SBmid", "sD", "sEndD"], "lineage2")

        d1.X = d1.X.todense()
        np.random.shuffle(d1.X)

        # TODO: make this into one function
        d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage1")
        d1 = dynchro.pp.calculate_pseudocells(d1, 100, "lineage2")
        d1 = dynchro.pp.merge_pseudocells_lineages(d1, ["lineage1", "lineage2"], 100)

        d2 = anndata.read("tests/data/batcheffect_dataseta0.h5ad")
        d2.uns["id"] = "dataset2"

        d2.X = d2.X.todense()
        # np.random.shuffle(d2.X)

        d2 = dynchro.pp.preprocess_ti(d2, root="sB")
        d2 = dynchro.pp.label_lineage(d2, "milestones", ["sA", "sB", "SBmid", "sC", "sEndC"], "lineage1")
        d2 = dynchro.pp.label_lineage(d2, "milestones", ["sA", "sB", "SBmid", "sD", "sEndD"], "lineage2")

        # TODO: make this into one function
        d2 = dynchro.pp.calculate_pseudocells(d2, 100, "lineage1")
        d2 = dynchro.pp.calculate_pseudocells(d2, 100, "lineage2")
        d2 = dynchro.pp.merge_pseudocells_lineages(d2, ["lineage1", "lineage2"], 100)

        self.test_shuffled(d1, d2, name="same_shuffled.png")

