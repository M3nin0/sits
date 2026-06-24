test_that("TWDTW - logistic weight (default)", {
    twdtw_model <- sits_train(samples_modis_ndvi, sits_twdtw())

    expect_equal(sits_bands(twdtw_model), "NDVI")
    expect_true(all(
        sits_labels(twdtw_model) %in% sits_labels(samples_modis_ndvi)
    ))

    point_ndvi <- sits_select(point_mt_6bands, bands = "NDVI")
    point_class <- sits_classify(
        data = point_ndvi,
        ml_model = twdtw_model,
        multicores = 1,
        progress = FALSE
    )

    # predicted classes belong to the training labels
    expect_true(all(
        point_class$predicted[[1]]$class %in% sits_labels(samples_modis_ndvi)
    ))

    # one prediction per yearly interval
    expect_true(nrow(sits_show_prediction(point_class)) == 17)

    # probabilities form a valid simplex (sum to 1 per row)
    probs <- point_class$predicted[[1]][, sits_labels(samples_modis_ndvi)]

    expect_true(all(
        abs(rowSums(probs) - 1) < 1.0e-06
    ))
})

test_that("TWDTW - no weight (plain DTW) and linear weight", {
    point_ndvi <- sits_select(point_mt_6bands, bands = "NDVI")

    model_none <- sits_train(samples_modis_ndvi, sits_twdtw(weight = "none"))
    class_none <- sits_classify(
        point_ndvi, model_none, multicores = 1, progress = FALSE
    )

    # classes must be valid labels
    expect_true(all(
        class_none$predicted[[1]]$class %in% sits_labels(samples_modis_ndvi)
    ))

    model_lin <- sits_train(
        samples = samples_modis_ndvi,
        ml_method = sits_twdtw(
            weight = "linear",
            a = 0.01,
            b = 0
        )
    )

    class_lin <- sits_classify(point_ndvi,
                               model_lin,
                               multicores = 1,
                               progress = FALSE)

    # classes must be valid labels
    expect_true(all(
        class_lin$predicted[[1]]$class %in% sits_labels(samples_modis_ndvi)
    ))
})

test_that("TWDTW - user-supplied patterns", {
    # Pre calculate patterns
    patterns <- sits_patterns(samples_modis_ndvi)

    # Train model with patterns
    twdtw_model <- sits_train(
        samples = samples_modis_ndvi,
        ml_method = sits_twdtw(
            patterns = patterns
        )
    )

    # Classify
    point_ndvi <- sits_select(point_mt_6bands, bands = "NDVI")
    point_class <- sits_classify(
        data = point_ndvi,
        ml_model = twdtw_model,
        multicores = 1,
        progress = FALSE
    )

    # classes must be valid labels
    expect_true(all(
        point_class$predicted[[1]]$class %in% sits_labels(samples_modis_ndvi)
    ))
})

test_that("TWDTW - GAM templates", {
    twdtw_model <- sits_train(
        samples = samples_modis_ndvi,
        ml_method = sits_twdtw(
            pattern_method = "gam"
        )
    )

    # Select NDVI
    point_ndvi <- sits_select(point_mt_6bands, bands = "NDVI")
    point_class <- sits_classify(
        data = point_ndvi,
        ml_model = twdtw_model,
        multicores = 1,
        progress = FALSE
    )

    # Classes must be valid labels
    expect_true(all(
        point_class$predicted[[1]]$class %in% sits_labels(samples_modis_ndvi)
    ))
})
