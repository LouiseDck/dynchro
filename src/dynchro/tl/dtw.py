from typing import Iterable
import numpy as np
import scipy.spatial

from anndata import AnnData

from collections import defaultdict

# to test
# should fail
#     reference, query and cost matrix
#     reference and cost matrix
#     query and cost matrix
#     reference and query with different number of features
#     - cost matrix and only_results False?
#     cost matrix and pseudotime_key
#     reference and query with unsorted pseudotime
#     reference and wrong pseudotime_key
#     query and wrong pseudotime_key
# should pass
#     reference and query with sorted, same key pseudotime
#     reference and query with mode "copy" or "inplace"
#     reference and query with mode "only_results"
#     cost_matrix with mode "only_results" and with no mode
#
def dtw(
        reference: AnnData | np.ndarray | None = None,
        query: AnnData | np.ndarray | None = None,
        cost_matrix: np.ndarray | None = None,
        dtw_key: str = "dtw",
        pseudotime_key: str = "pseudotime",
        distance = "euclidean",
        mode: str = "copy"
        ) -> tuple[np.ndarray, np.ndarray, float] | tuple[AnnData, AnnData] | None:
    """
    Computes Dynamic Time Warping (DTW) distance between two datasets.

    Parameters
    ----------
    reference : Union[AnnData, np.ndarray]
        The reference dataset, either as an AnnData object or a numpy array.
        You must provide either `reference` and `query`, or `cost_matrix` but not both.
    query : Union[AnnData, np.ndarray]
        The query dataset, either as an AnnData object or a numpy array.
        You must provide either `reference` and `query`, or `cost_matrix` but not both.
    cost_matrix : Union[np.ndarray, None], optional
        A precomputed cost matrix to use for DTW. If None, it will be computed
        using the specified distance metric, by default None.
    dtw_key : str, optional
        The key under which the DTW results will be stored in the AnnData object, by default "dtw".
    pseudotime_key : str, optional
        The key in `obs` of the AnnData object where pseudotime values are stored, by default "pseudotime".
        If the input is a numpy array, this parameter is ignored.
    distance : str, optional
        The distance metric to use for DTW, by default "euclidean".
        Other options can be "manhattan", "cosine", etc., as supported by `scipy.spatial.distance.cdist`.
    mode : str, optional
        The mode of operation, either "copy", "inplace" or "only_results".
        If "copy", the function returns two AnnData objects with DTW results added.
        If "inplace", it modifies the input AnnData objects in place.
        If "only_results", it returns only the DTW results without modifying the input objects.
    copy : bool, optional
        If True, returns copies of the input AnnData objects with DTW results added.
        If False, modifies the input objects in place, by default True.
    only_results : bool, optional
        If True, returns only the DTW distance, cost matrix, and distance matrix without modifying the (optional)
        AnnData objects. If no
        If False, modifies the AnnData objects to include the DTW results, by default False.

    """

    # check if either reference and query or cost_matrix is provided
    assert (reference is not None and query is not None) or (cost_matrix is not None), \
        "You must provide either `reference` and `query`, or `cost_matrix` but not both."

    assert mode in ["copy", "inplace", "only_results"], \
        "`mode` must be one of 'copy', 'inplace', or 'only_results'."

    if reference is not None and query is not None:
        if isinstance(reference, AnnData) and isinstance(query, AnnData):
            assert reference.shape[1] == query.shape[1], \
                "Reference and query AnnData objects must have the same number of features (columns)."

            # check if pseudotime_key exists in both AnnData objects
            assert pseudotime_key in reference.obs and pseudotime_key in query.obs, \
                f"Both AnnData objects must contain the key '{pseudotime_key}' in their obs."

            # check if pseudotime values are sorted
            assert np.all(np.diff(reference.obs[pseudotime_key].values) >= 0), \
                f"Pseudotime values in reference AnnData object must be sorted in ascending order for key '{pseudotime_key}'."
            assert np.all(np.diff(query.obs[pseudotime_key].values) >= 0), \
                f"Pseudotime values in query AnnData object must be sorted in ascending order for key '{pseudotime_key}'."

            # sort the values in both AnnData objects
            reference_matrix = get_values_sorted(reference, pseudotime_key).X.toarray()
            query_matrix = get_values_sorted(query, pseudotime_key).X.toarray()

        elif isinstance(reference, np.ndarray) and isinstance(query, np.ndarray):
            assert reference.shape[1] == query.shape[1], \
                "Reference and query numpy arrays must have the same number of features (columns)."

        dtw_distance, cost, D = _dtw(reference_matrix, query_matrix, distance)

    else:
        assert isinstance(cost_matrix, np.ndarray), \
            "If `cost_matrix` is provided, it must be a 2D numpy array."
        assert cost_matrix.ndim == 2, "Cost matrix must be a 2D array."

        dtw_distance, cost, D = _dtw(cost_matrix)

    if mode == "only_results":
        return dtw_distance, cost, D

    reference.obsm[f"{dtw_key}_cost"] = cost
    reference.obsm[f"{dtw_key}_D"] = D
    reference.uns[f"{dtw_key}_distance"] = dtw_distance

    query.obsm[f"{dtw_key}_cost"] = cost.T
    query.obsm[f"{dtw_key}_D"] = D.T
    query.uns[f"{dtw_key}_distance"] = dtw_distance

    if mode == "copy":
        return reference, query

    return None

def traceback(
        reference: AnnData | None = None,
        query: AnnData | None = None,
        dtw_key: str = "dtw",
        D: np.ndarray | None = None,
        end_x: int | str | None = None,
        end_y: int | str | None = None,
        mode: str = "copy"
        ) -> tuple[AnnData, AnnData] | tuple[np.ndarray, np.ndarray] | None:
    """
    Traceback the DTW distance matrix to find the optimal path.

    Parameters
    ----------
    """
    if D is not None:
        mode = "only_results"

    assert (reference is not None and query is not None) or (D is not None), \
        "You must provide either `reference` and `query`, or `D` but not both."

    if reference is not None and query is not None:
        assert f"{dtw_key}_D" in reference.obsm and f"{dtw_key}_D" in query.obsm, \
            f"Both AnnData objects must contain the key '{dtw_key}_D' in their obsm."
        assert np.array_equal(reference.obsm[f"{dtw_key}_D"], query.obsm[f"{dtw_key}_D"].T), \
            "The DTW distance matrices in both AnnData objects must be the same."

        D = reference.obsm[f"{dtw_key}_D"]

    if end_x is None:
        end_x = D.shape[0] - 1
    if end_y is None:
        end_y = D.shape[1] - 1

    if isinstance(end_x, str):
        assert end_x == "free", \
            "`end_x` must be either an integer, 'free', or None."
        assert isinstance(end_y, int), \
            "`end_y` must be an integer or None if `end_x` is 'free'."
        end_x = np.argmin(D[:, end_y])

    elif isinstance(end_y, str):
        assert end_y == "free", \
            "`end_y` must be either an integer, 'free', or None."
        assert isinstance(end_x, int), \
            "`end_x` must be an integer or None if `end_y` is 'free'."
        end_y = np.argmin(D[end_x, :])

    assert 0 <= end_x < D.shape[0] and 0 <= end_y < D.shape[1], \
        (
            "`end_x` and `end_y` must be integers within the range of the distance matrix.",
            f"`end_x` must be between 0 and {D.shape[0] - 1}, and `end_y` must be between 0 and {D.shape[1] - 1}."
        )
    path1, path2 = traceback_start(D, end_x, end_y)
    agg_path1, agg_path2 = aggregate_path(path1, path2)

    if mode == "only_results":
        return agg_path1, agg_path2

    reference.obs[f"{dtw_key}_path"] = list(agg_path1.values()) + [None] * (len(reference.obs) - len(agg_path1.values()))
    query.obs[f"{dtw_key}_path"] = list(agg_path2.values()) + [None] * (len(query.obs) - len(agg_path2.values()))

    cells_path1 = assign_cells_from_path(agg_path1, query)
    cells_path2 = assign_cells_from_path(agg_path2, reference)

    reference.obs[f"{dtw_key}_cells_path"] = list(cells_path1.values()) + [None] * (len(reference.obs) - len(cells_path1.values()))
    query.obs[f"{dtw_key}_cells_path"] = list(cells_path2.values()) + [None] * (len(query.obs) - len(cells_path2.values()))

    if mode == "copy":
        return reference, query

    return None

def _dtw(x, y, distance="euclidean"):
    cost = scipy.spatial.distance.cdist(x, y, distance)
    tol = 1e-15
    cost[abs(cost) < tol] = 0

    return _dtw_cost(cost)

def _dtw_cost(cost):
    r, c = cost.shape
    distances = np.zeros((r + 1, c + 1))
    distances[0, 1:], distances[1:, 0] = np.inf, np.inf
    distances[0, 0] = 0

    for i in range(1, r + 1):
        for j in range(1, c + 1):
            options = [distances[i - 1, j - 1], distances[i - 1, j], distances[i, j - 1]]
            distances[i, j] = min(options) + cost[i - 1, j - 1]

    return distances[-1, -1], cost, distances[1:, 1:]


def traceback_start(D, i, j):
    p, q = [i], [j]
    while (i > 0) or (j > 0):
        tb = np.argmin((D[i - 1, j - 1], D[i, j - 1], D[i - 1, j]))
        if tb == 0:
            i -= 1
            j -= 1
        elif tb == 1:
            j -= 1
        elif tb == 2:
            i -= 1
        p.insert(0, i)
        q.insert(0, j)
    return np.array(p), np.array(q)


def _traceback(D):
    i, j = np.array(D.shape) - 1

    return traceback_start(D, i, j)


def traceback_yedge(D):
    i = np.argmin(D[:, -1])
    j = D.shape[1] - 1

    return traceback_start(D, i, j)


def traceback_xedge(D):
    j = np.argmin(D[-1, :])
    i = D.shape[0] - 1

    return traceback_start(D, i, j)


def get_values_sorted(adata: AnnData, pseudotime_key: str = "pseudotime") -> np.ndarray:
    """
    Sorts the matrix of an AnnData object based on the pseudotime values stored in `adata.obs`.

    Parameters
    ----------
    adata : AnnData
        The AnnData object containing the pseudotime values.
    pseudotime_key : str, optional
        The key in `adata.obs` where the pseudotime values are stored, by default "pseudotime".
        The result will be sorted according to these values.
        The values are expected to be numeric.

    Returns
    -------
    np.ndarray
        Sorted pseudotime values according to the specified key.
    """
    if pseudotime_key not in adata.obs:
        raise ValueError(f"pseudotime_key '{pseudotime_key}' not found in AnnData.obs")

    order = np.argsort(adata.obs[pseudotime_key].values)

    return adata[order, :]


def aggregate_path(path1: Iterable[int], path2: Iterable[int]):

    reference_path = defaultdict(list)
    query_path = defaultdict(list)

    for i, j in zip(path1, path2):
        reference_path[i].append(j)
        query_path[j].append(i)

    return reference_path, query_path

def assign_cells_from_path(
        agg_path: dict[int, list[int]],
        data: AnnData
) -> dict[int, list[str]]:
    """
    Assigns cells to their corresponding paths based on the aggregated path.

    Parameters
    ----------
    agg_path : dict[int, list[int]]
        The aggregated path where keys are indices and values are lists of indices.
    data : AnnData
        The AnnData object containing the cell information.

    Returns
    -------
    dict[int, list[str]]
        A dictionary where keys are indices and values are lists of cell names.
    """
    cell_paths = defaultdict(list)

    for i, cells in agg_path.items():
        for cell in cells:
            cell_paths[i].append(data.obs_names[cell])

    return cell_paths
