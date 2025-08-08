from anndata import AnnData

from dynchro.tl.dtw import dtw
from dynchro.tl.dtw import traceback
from dynchro.tl.warp import warp_pseudotime

from dynchro.pl.matrices import plot_dtw_matrices

def dynchronize(
    reference: AnnData | None = None,
    query: AnnData | None = None,
    pseudotime_key: str = "pseudotime",
    dtw_key: str = "dtw",
    mode: str = "copy"
):
    # run dtw, traceback, warp_pseudotime and plot matrices
    reference, query = dtw(reference, query, dtw_key = dtw_key, pseudotime_key = pseudotime_key, mode = mode)
    reference, query = traceback(reference, query, dtw_key = dtw_key, mode = mode)
    query = warp_pseudotime(reference, query, dtw_key = dtw_key, pseudotime_key = pseudotime_key, mode = mode)
    reference = warp_pseudotime(query, reference, dtw_key = dtw_key, pseudotime_key = pseudotime_key, mode = mode)

    plot_dtw_matrices(reference, query, dtw_key)
