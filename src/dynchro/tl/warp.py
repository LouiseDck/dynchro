from collections import defaultdict

from anndata import AnnData

def warp_pseudotime(
        reference: AnnData,
        query: AnnData,
        dtw_key: str = "dtw",
        pseudotime_key: str = "pseudotime",
        mode: str = "copy"
) -> tuple[AnnData, AnnData] | None:

    warpings = query.obs[f"{dtw_key}_cells_path"]

    warped_pseudotime = defaultdict(list)

    for index, cell in zip(warpings.index, warpings):
        # get pseudotime from cell in reference
        pseudotime = reference.obs.loc[cell, pseudotime_key].values[0]
        warped_pseudotime[index].append(pseudotime)

    # average the pseudotime values
    for index, values in warped_pseudotime.items():
        warped_pseudotime[index] = sum(values) / len(values)

    if mode == "only_results":
        return warped_pseudotime

    query.obs[f"{dtw_key}_warped_pseudotime"] = warped_pseudotime

    if mode == "copy":
        return query

    return None
