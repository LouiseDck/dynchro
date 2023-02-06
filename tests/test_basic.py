import anndata
import pytest

import dynchro


def test_package_has_version():
    dynchro.__version__


@pytest.mark.skip(reason="This decorator should be removed when test passes.")
def test_example():
    assert 1 == 0  # This test is designed to fail.


@pytest.fixture
def dataset1():
    yield anndata.read("dyngen_pooled/data/batcheffect_dataseta0.h5ad")


@pytest.fixture
def dataset2():
    yield anndata.read("dyngen_pooled/data/batcheffect99_datasetb0.h5ad")


@pytest.mark.skip(reason="Not implemented yet")
class TestTwoSimpleDatasets:
    def test_pseudocells(self):
        assert 1 == 0

    def test_merge_pseudocells(self):
        assert 1 == 0

    def test_align_pseudocells(self):
        assert 1 == 0

    def test_align_cells(self):
        assert 1 == 0
