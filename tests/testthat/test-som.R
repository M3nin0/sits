test_that("Calculating time-series distances using DTW", {
    set.seed(2903)

    # test 1 (example from `IncDTW` article)
    t1 <- c(3, 4, 5, 6)
    t2 <- c(1, 3, 3, 5, 6)

    expect_equal(dtw2vec(t1, t2), 3)

    # test 2 (using ndvi time-series)
    t1 <- c(0.5001, 0.5060, 0.5079, 0.7141, 0.4262)
    t2 <- c(0.3974, 0.5144, 0.5052, 0.6463, 0.5139)

    # 0.4323 -> Value confirmed using `dtw` and `IncDTW` packages as reference
    expect_equal(dtw2vec(t1, t2), 0.4323)
})


test_that("Creating clustering using Self-organizing Maps with DTW distance", {
    set.seed(2903)

    som_map <- sits_som_map(
        samples_modis_ndvi,
        grid_xdim = 4,
        grid_ydim = 4,
        distance  = "dtw" # custom distance only available in the sits package
    )

    expect_true(all(colnames(som_map$labelled_neurons) %in%
        c("id_neuron", "label_samples", "count", "prior_prob", "post_prob")))

    expect_true(som_map$labelled_neurons[1, ]$prior_prob >= 0)
    expect_true(som_map$labelled_neurons[1, ]$post_prob >= 0)
    expect_true(all(unique(som_map$labelled_neurons$id_neuron) %in% 1:16))

    cleaned_samples <- sits_som_clean_samples(som_map)
    expect_true("eval" %in% names(cleaned_samples))
    expect_true("post_prob" %in% names(cleaned_samples))
    expect_true(all(cleaned_samples$eval %in% c("clean", "analyze", "remove")))

    expect_true(cleaned_samples[1, ]$post_prob > 0)

    cluster_purity <- suppressMessages(sits_som_evaluate_cluster(som_map))

    expect_true(cluster_purity[1, ]$mixture_percentage > 60)
    expect_true(cluster_purity[2, ]$mixture_percentage < 40)
    expect_error(sits_som_clean_samples(samples_modis_ndvi))
    expect_error(sits_som_evaluate_samples(samples_modis_ndvi))
})
