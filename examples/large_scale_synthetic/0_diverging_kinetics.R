library(tidyverse)
library(dyngen)

# Generate a batch effect in a dataset, using kinetics
# Kinetics are in between the two models, depending on the alpha used
# Alpha closer to 0 -> closer to model1
generate_diverging_kinetics <- function(model1, model2, alpha) {
    print("here")
    model_between <- model1
    fii <- model1$feature_info
    fni <- model1$feature_network

    genes <- c()

    ## TODO the branch specificness does not work :'(
    print("Network")
    for (i in seq(from = 1, to = nrow(fni))) {
        if (grepl("D[0-9]*", fni[i, "from"])) {
            genes <- c(genes, fni[i, ]$from)
            genes <- c(genes, fni[i, ]$to)
            genes <- unique(genes)
            print(genes)
              for (col in colnames(fii)) {
                  if (is.numeric(model1$feature_network[[col]]) && !all(model1$feature_network[[col]] == model2$feature_network[[col]])) {
                      fii[[i, col]] <- model1$feature_network[[i, col]] * (1 - alpha) + model2$feature_network[[i, col]] * alpha
                  }
              }
          }
    }

    # TODO: maybe this can be done better than with a double for loop?
    # Maybe subset the columns first & don't loop over all columns?
    for (i in seq(from = 1, to = nrow(fii))) {
        # if (grepl("D[0-9]*", fii[i, "module_id"])) {
        # print(fii[i, "feature_id"])
        if (fii[i, ]$feature_id %in% genes) {
            print(i)
            for (col in colnames(fii)) {
                if (is.numeric(model1$feature_info[[col]]) && !all(model1$feature_info[[col]] == model2$feature_info[[col]])) {
                    fii[[i, col]] <- model1$feature_info[[i, col]] * (1 - alpha) + model2$feature_info[[i, col]] * alpha
                }
            }
        }
    }


    # for (col in colnames(fii)) {
    #     if (is.numeric(model1$feature_info[[col]]) && !all(model1$feature_info[[col]] == model2$feature_info[[col]])) {
    #         fii[[col]] <- model1$feature_info[[col]] * (1 - alpha) + model2$feature_info[[col]] * alpha
    #     }
    # }

    # for (col in colnames(fni)) {
    #     if (is.numeric(model1$feature_network[[col]]) && !all(model1$feature_network[[col]] == model2$feature_network[[col]])) {
    #         fni[[col]] <- model1$feature_network[[col]] * (1 - alpha) + model2$feature_network[[col]] * alpha
    #     }
    # }

    model_between$feature_info <- fii
    model_between$feature_network <- fni

    model_between <- model_between %>%
        generate_gold_standard() %>%
        generate_cells() %>%
        generate_experiment()

    model_between
}
