import os

import anndata
import pytest

import dynchro


@pytest.fixture(scope="session")
def dataset1():
    print(os.getcwd())
    yield anndata.read("tests/data/batcheffect_dataseta0.h5ad")


@pytest.fixture(scope="session")
def dataset2():
    yield anndata.read("tests/data/batcheffect99_datasetb0.h5ad")


class TestPreprocessing:
    def test_preprocess_ti(self, dataset1):
        dynchro.pp.preprocess_ti(dataset1, root="sB")

        # check that data is normalised
        # TODO: how to?
        # check that data is log transformed
        assert "log1p" in dataset1.uns, "Counts are not log transformed"
        # check that data is scaled
        # TODO: how to?
        # check neighbors have been computed
        assert "neighbors" in dataset1.uns, "Neighbors have not been computed"
        # check umap has been computed
        assert "X_umap" in dataset1.obsm, "Umap has not been computed"
        # check that paga has been run
        assert "paga" in dataset1.uns, "Paga has not been run"
        # check that root has been set
        assert "iroot" in dataset1.uns, "Root has not been set"
        # check that dpt has been run
        assert "pseudotime" in dataset1.obs, "DPT has not been run, or 'pseudotime' not added to obs"

    def test_preprocess_batch_effect(self, dataset2):
        preprocessed = dynchro.pp.preprocess_batch_effect(dataset2)

        assert "sEndC" not in preprocessed.obs.milestones.unique(), "sEndC is still in the dataset"

        # check that label has been removed

    def test_labelling(self, dataset1):
        dataset1 = dynchro.pp.label_lineage(dataset1, "milestones", ["sB", "sEndC"], "lineage1")
        # check that the right column name has been added

        assert "lineage1" in dataset1.obs, "lineage1 column not added to obs"
        assert "lineage1" in dataset1.uns["lineage_labels"], "lineage1 not added to uns[lineage_labels]"

        # check that the right value has been added
        assert dataset1.obs.loc[
            dataset1.obs["milestones"] == "sB", "lineage1"
        ].all(), "not all sB cells labelled as lineage1"
        assert dataset1.obs.loc[
            dataset1.obs["milestones"] == "sEndC", "lineage1"
        ].all(), "not all sEndC cells labelled as lineage1"


class TestPseudocells:
    def test_spacing(self, dataset1):
        assert 1 == 0

    def test_interpolation(self, dataset1):
        assert 1 == 0

    def test_merge(self, dataset1):
        assert 1 == 0
