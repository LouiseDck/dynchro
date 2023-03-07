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
            "copmare_lineages_pseudocells_label": "pseudocells_100",
            "compare_trajectories_pseudocells": False,
            "compare_trajectories_pseudocells_label": "pseudocells_100",
            "align_pseudocells": False,
        }

        matching_lineages = dynchro.tl.get_matching_lineages([dataset1, dataset2], config=config)
        aligned = dynchro.tl.align_trajectories(matching_lineages, [dataset1, dataset2], pseudocells=False)

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

    @pytest.mark.skip(reason="not implemented yet")
    def test_align_trajectories_simple_pseudocells(self, dataset1, dataset2):
        assert 1 == 0

    @pytest.mark.skip(reason="not implemented yet")
    def test_align_trajectories_complex_pseudocells(self, dataset_complex1, dataset_complex2):
        assert 1 == 0
