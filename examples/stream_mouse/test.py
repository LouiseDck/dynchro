def map_new_data(adata_ref, adata_new, use_radius=False):
    """Map new data to the inferred trajectories

    Parameters
    ----------
    adata_ref: reference AnnData
        Annotated data matrix.
    adata_new: new AnnData
        Annotated data matrix for new data (to be mapped).
    use_radius: `bool`, optional (default: False)
        Only valid when `method = 'mlle'`. If True, when searching for the neighbors for each cell in MLLE space, STREAM uses a fixed radius instead of a fixed number of cells.
    first_pc: `bool`, optional (default: False)
        Only valid when `feature='top_pcs'` If True, the first principal component will be included
    top_pcs_feature: `str`, optional (default: None)
        Choose from {{'var_genes'}}
        Only valid when `feature='top_pcs'`.
        Features used for pricipal component analysis
        If None, all the genes will be used.
        IF 'var_genes', the most variable genes obtained from select_variable_genes() will be used.
    Returns
    -------

    Combined AnnData object. Only the shared obs_keys and var_keys will be preserved
    updates `adata` with the following fields.
    batch: `pandas.Series` (`adata.obs['batch']`,dtype `str`)
        The annotation of each cell. It consists of 'ref' for reference cells and 'new' for the cells to map

    updates `adata_new` with the following fields.(depending on the `feature` or `method`)
    var_genes: `numpy.ndarray` (`adata_new.obsm['var_genes']`)
        Store #observations × #var_genes data matrix used mapping.
    var_genes: `pandas.core.indexes.base.Index` (`adata_new.uns['var_genes']`)
        The selected variable gene names.
    top_pcs: `numpy.ndarray` (`adata_new.obsm['top_pcs']`)
        Store #observations × n_pc data matrix used for subsequent dimension reduction.
    X_dr : `numpy.ndarray` (`adata_new.obsm['X_dr']`)
        A #observations × n_components data matrix after dimension reduction.
    X_mlle : `numpy.ndarray` (`adata_new.obsm['X_mlle']`)
        Store #observations × n_components data matrix after mlle.
    X_umap : `numpy.ndarray` (`adata_new.obsm['X_umap']`)
        Store #observations × n_components data matrix after umap.
    X_pca : `numpy.ndarray` (`adata_new.obsm['X_pca']`)
        Store #observations × n_components data matrix after pca.
    X_vis_umap: `numpy.ndarray` (`adata.obsm['X_vis_umap']`)
        Store #observations × 2 umap data matrix.
    X_vis_tsne: `numpy.ndarray` (`adata.obsm['X_vis_tsne']`)
        Store #observations × 2 tsne data matrix.
    """

    feature = adata_ref.uns["params"]["dimension_reduction"]["feature"]
    method = adata_ref.uns["params"]["dimension_reduction"]["method"]
    assert method in ["mlle", "umap", "pca"], (
        "'%s' is not supported. Please use one of 'mlle','umap', and 'pca' for dimension reduction." % method
    )
    if feature == "var_genes":
        print("Top variable genes are being used for mapping ...")
        if "var_genes" not in adata_ref.uns_keys():
            raise ValueError("variable genes are not selected yet in `adata_ref`")
        for x in adata_ref.uns["var_genes"]:
            if x not in adata_new.var_names:
                raise ValueError("variable gene '%s' does not exist in `adata_new.var_names`" % (x))
        adata_new.uns["var_genes"] = adata_ref.uns["var_genes"].copy()
        adata_new.obsm["var_genes"] = adata_new[:, adata_new.uns["var_genes"]].X.copy()
        input_data = adata_new.obsm["var_genes"]
    # if(feature == 'all'):
    #     print('All genes are being used for mapping ...')
    #     if(not set(adata_ref.var_names) <= set(adata_new.var_names)):
    #         raise ValueError("`adata_new.var_names` does not contain all the genes in `adata_ref.var_names`")
    #     input_data = adata_new[:,adata_ref.var.index].X
    # if(feature == 'top_pcs'):
    #     print('Top principal components are being used for mapping ...')
    #     if('top_pcs' not in adata_ref.uns_keys()):
    #         raise ValueError("top principal components are not selected yet in `adata_ref`")
    #     trans = adata_ref.uns['top_pcs']
    #     top_pcs_feature = adata_ref.uns['params']['select_top_principal_components']['feature']
    #     if(top_pcs_feature == 'var_genes'):
    #         if('var_genes' not in adata_ref.uns_keys()):
    #             raise ValueError("variable genes are not selected yet in `adata_ref`")
    #         adata_new.uns['var_genes'] = adata_ref.uns['var_genes'].copy()
    #         adata_new.obsm['var_genes'] = adata_new[:,adata_new.uns['var_genes']].X.copy()
    #         X_pca = trans.transform(adata_new.obsm['var_genes'])
    #     else:
    #         X_pca = trans.transform(adata_new[:,adata_ref.var.index].X)
    #     n_pc = adata_ref.obsm['top_pcs'].shape[1]
    #     first_pc = adata_ref.uns['params']['select_top_principal_components']['first_pc']
    #     if(first_pc):
    #         adata_new.obsm['top_pcs'] = X_pca[:,0:(n_pc)]
    #     else:
    #         #discard the first Principal Component
    #         adata_new.obsm['top_pcs'] = X_pca[:,1:(n_pc+1)]
    #     input_data = adata_new.obsm['top_pcs']
    adata_new.uns["epg"] = adata_ref.uns["epg"].copy()
    adata_new.uns["flat_tree"] = adata_ref.uns["flat_tree"].copy()

    print("method '%s' is being used for mapping ..." % method)
    if method == "mlle":
        if "trans_mlle" in adata_ref.uns_keys():
            trans = adata_ref.uns["trans_mlle"]
            if use_radius:
                dist_nb = trans.nbrs_.kneighbors(input_data, n_neighbors=trans.n_neighbors, return_distance=True)[0]
                ind = trans.nbrs_.radius_neighbors(input_data, radius=dist_nb.max(), return_distance=False)
                new_X_mlle = np.empty((input_data.shape[0], trans.n_components))
                for i in range(input_data.shape[0]):
                    weights = barycenter_weights_modified(input_data[i], trans.nbrs_._fit_X[ind[i]], reg=trans.reg)
                    new_X_mlle[i] = np.dot(trans.embedding_[ind[i]].T, weights)
                adata_new.obsm["X_mlle_mapping"] = new_X_mlle
            else:
                adata_new.obsm["X_mlle_mapping"] = trans.transform(input_data)
            adata_new.obsm["X_dr"] = adata_new.obsm["X_mlle_mapping"].copy()
        else:
            raise Exception("Please run 'st.dimension_reduction()' using 'mlle' first.")
    # if(method == 'umap'):
    #     if('trans_umap' in adata_ref.uns_keys()):
    #         trans = adata_ref.uns['trans_umap']
    #         adata_new.obsm['X_umap_mapping'] = trans.transform(input_data)
    #         adata_new.obsm['X_dr'] = adata_new.obsm['X_umap_mapping'].copy()
    #     else:
    #         raise Exception("Please run 'st.dimension_reduction()' using 'umap' first.")
    # if(method == 'pca'):
    #     if('trans_pca' in adata_ref.uns_keys()):
    #         trans = adata_ref.uns['trans_pca']
    #         adata_new.obsm['X_pca_mapping'] = trans.transform(input_data)
    #         adata_new.obsm['X_dr'] = adata_new.obsm['X_pca_mapping'].copy()
    #     else:
    #         raise Exception("Please run 'st.dimension_reduction()' using 'pca' first.")

    if "plot_visualization_2D" in adata_ref.uns["params"].keys():
        print("Visualizing new cells on 2D plane ...")
        vis_method = adata_ref.uns["params"]["plot_visualization_2D"]["method"]
        vis_trans = "vis_trans_" + vis_method
        trans = adata_ref.uns[vis_trans]
        adata_new.obsm["X_vis_" + vis_method] = trans.transform(adata_new.obsm["X_dr"])
        if adata_ref.uns["params"]["epg"]["use_vis"]:
            print("Using the manifold from `plot_visualization_2D()` ")
            adata_new.obsm["X_dr"] = adata_new.obsm["X_vis_" + vis_method]

    project_cells_to_epg(adata_new)
    calculate_pseudotime(adata_new)
    adata_combined = adata_ref.concatenate(adata_new, batch_categories=["ref", "new"])
    shared_obs_key = [x for x in adata_new.obs_keys() if x in adata_ref.obs_keys()]
    shared_var_key = [x for x in adata_new.var_keys() if x in adata_ref.var_keys()]
    adata_combined.obs = adata_combined.obs[shared_obs_key + ["batch"]]
    adata_combined.var = adata_combined.var[shared_var_key]
    adata_combined.uns["workdir"] = adata_new.uns["workdir"]
    for key in adata_ref.uns_keys():
        if key in ["var_genes", "epg", "flat_tree", "vis_trans_tsne", "vis_trans_umap"]:
            adata_combined.uns[key] = adata_ref.uns[key]
        if key.split("_")[-1] == "color":
            if key in adata_new.uns_keys():
                adata_combined.uns[key] = {**adata_ref.uns[key], **adata_new.uns[key]}
    return adata_combined
