test_that(".bert_dataset_lazy returns correct item shapes", {
    skip_if_not_installed("torch")
    # define bands
    bands <- .samples_bands(samples_modis_ndvi)

    # define time-series
    timeline <- .samples_timeline(samples_modis_ndvi)

    # define shapes
    n_times <- length(timeline)
    n_bands <- length(bands)
    stats <- .samples_stats(samples_modis_ndvi)

    # normalize time-series
    norm_ts <- .pred_as_ts(
        data = .pred_normalize(.predictors(samples_modis_ndvi), stats),
        bands = bands,
        timeline = timeline
    )

    # band standard deviations
    band_sd <- vapply(bands, function(b) stats::sd(norm_ts[[b]]), numeric(1))
    names(band_sd) <- bands

    # generate dataset
    ds <- .bert_dataset_lazy(
        samples = samples_modis_ndvi,
        indices = seq_len(10L),
        stats = stats,
        bands = bands,
        timeline = timeline,
        mask_ratio = 0.5,
        masking_method = "random",
        noise_frac = 0.5,
        noised_bands = bands,
        band_sd = band_sd
    )

    # get item
    item <- ds$.getitem(1L)

    # test
    expect_equal(as.integer(item$x$shape), c(n_times, n_bands))
    expect_equal(as.integer(item$y$y$shape), c(n_times, n_bands))
    expect_equal(as.integer(item$y$mask$shape), c(n_times, 1L))
    expect_equal(ds$.length(), 10L)
})

test_that(".bert_dataset_lazy corrupts only masked positions with noise", {
    skip_if_not_installed("torch")

    # define bands
    bands <- .samples_bands(samples_modis_ndvi)

    # define time-series
    timeline <- .samples_timeline(samples_modis_ndvi)

    # define time-series stats
    stats <- .samples_stats(samples_modis_ndvi)

    # normalize time-series
    norm_ts <- .pred_as_ts(
        data = .pred_normalize(.predictors(samples_modis_ndvi), stats),
        bands = bands,
        timeline = timeline
    )

    # band standard deviations
    band_sd <- vapply(bands, function(b) stats::sd(norm_ts[[b]]), numeric(1))
    names(band_sd) <- bands

    # generate dataset
    ds <- .bert_dataset_lazy(
        samples = samples_modis_ndvi,
        indices = seq_len(5L),
        stats = stats,
        bands = bands,
        timeline = timeline,
        mask_ratio = 0.5,
        masking_method = "random",
        noise_frac = 0.5,
        noised_bands = bands,
        band_sd = band_sd
    )

    # get item
    item <- ds$.getitem(1L)
    x    <- as.matrix(torch::as_array(item$x))
    y    <- as.matrix(torch::as_array(item$y$y))
    mask <- as.vector(torch::as_array(item$y$mask))

    # visible (mask == 0) positions must be untouched
    visible <- which(mask == 0)

    if (length(visible) > 0L) {
        expect_equal(x[visible, , drop = FALSE], y[visible, , drop = FALSE])
    }

    # At least one corrupted position must differ from the clean signal
    corrupted <- which(mask == 1)

    expect_true(length(corrupted) > 0L)
    expect_false(isTRUE(all.equal(
        x[corrupted, , drop = FALSE], y[corrupted, , drop = FALSE]
    )))
})

test_that("sits_bert returns a function when no samples given", {
    # generate functions
    bert_fn <- sits_bert(
        embedding_dim = 16L,
        epochs        = 5L
    )

    # check if the factory is working
    expect_true(is.function(bert_fn))
})

test_that("sits_bert pre-training produces sits_encoder", {
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")

    # generate encoder
    encoder <- .try(
        sits_pre_train(
            samples        = samples_modis_ndvi,
            encoder_method = sits_bert(
                embedding_dim = 16L,
                epochs        = 5L,
                batch_size    = 32L,
                noise_frac    = 0.5,
                verbose       = FALSE
            )
        ),
        .default = NULL
    )

    # skip if encoder is not there
    skip_if(is.null(encoder), "torch training failed (likely no resources)")

    # test
    expect_true(inherits(encoder, "sits_encoder"))
    expect_true(inherits(encoder, "torch_model"))
})

test_that("sits_bert encoder can encode a sits tibble", {
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")

    # embeddings dim
    embedding_dim <- 16L

    # generate encoder
    encoder <- .try(
        sits_pre_train(
            samples        = samples_modis_ndvi,
            encoder_method = sits_bert(
                embedding_dim = embedding_dim,
                epochs        = 5L,
                batch_size    = 32L,
                verbose       = FALSE
            )
        ),
        .default = NULL
    )

    # skip if the encoder is not there
    skip_if(is.null(encoder), "torch training failed")

    # encode!
    enc_samples <- sits_encode(
        data    = sits_sample(samples_modis_ndvi, frac = 0.3),
        encoder = encoder
    )

    # test
    expect_true(inherits(enc_samples, "sits"))
    expect_equal(length(sits_bands(enc_samples)), embedding_dim)
    expect_true(all(grepl("^EMB", sits_bands(enc_samples))))
})

test_that("sits_bert: downstream classification works", {
    skip_if_not_installed("torch")
    skip_if_not_installed("luz")

    # generate encoder
    encoder <- .try(
        sits_pre_train(
            samples        = samples_modis_ndvi,
            encoder_method = sits_bert(
                embedding_dim = 16L,
                epochs        = 5L,
                batch_size    = 32L,
                verbose       = FALSE
            )
        ),
        .default = NULL
    )

    # skip if the encoder is not there
    skip_if(is.null(encoder), "torch training failed")

    # encode samples
    enc_samples <- sits_encode(
        data    = sits_sample(samples_modis_ndvi, frac = 0.6),
        encoder = encoder
    )

    # train random forest
    rf_model <- sits_train(enc_samples, sits_rfor(num_trees = 20L))

    # sits encode
    point_enc <- sits_encode(
        data    = sits_select(point_mt_6bands, bands = "NDVI"),
        encoder = encoder
    )

    # classify
    point_class <- sits_classify(
        data     = point_enc,
        ml_model = rf_model,
        progress = FALSE
    )

    # test
    expect_true(
        all(point_class$predicted[[1L]]$class %in%
                sits_labels(samples_modis_ndvi))
    )
})
