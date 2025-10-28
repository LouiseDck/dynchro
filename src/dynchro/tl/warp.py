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
        if cell is not None:
            # print(index, cell)
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


def avg_cost_path(
        query: AnnData,
        dtw_key: str = "dtw",
        mode: str = "copy"
) -> tuple[AnnData, AnnData] | None:

    cost_matrix = query.obsm[f"{dtw_key}_cost"]
    warpings = query.obs[f"{dtw_key}_path"]

    avg_cost = defaultdict(list)

    for i, (index, cells) in enumerate(zip(warpings.index, warpings)):
        if cells is not None:
            for cell in cells:

                # get cost from cell in reference
                cost = cost_matrix[i, cell]
                avg_cost[index].append(cost)


    for index, value in avg_cost.items():
        avg_cost[index] = sum(value) / len(value)

    # averaged_cost = [sum(x) / len(x) for x in avg_cost.values()]

    if mode == "only_results":
        return avg_cost

    print(avg_cost)
    query.obs[f"{dtw_key}_avg_cost"] = avg_cost

    if mode == "copy":
        return query

    return None

