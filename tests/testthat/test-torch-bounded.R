test_that("bounded module reproduces the unwrapped module", {
    skip_on_cran()
    skip_if_not_installed("torch")

    torch::torch_manual_seed(42L)
    # A module shaped like a sits classifier:
    # (batch, n_times, n_bands) -> (batch, n_labels)
    module <- torch::nn_module(
        initialize = function() {
            self$net <- torch::nn_sequential(
                torch::nn_flatten(),
                torch::nn_linear(24L * 3L, 16L),
                torch::nn_relu(),
                torch::nn_linear(16L, 5L),
                torch::nn_softmax(dim = 2L)
            )
        },
        forward = function(x) self$net(x)
    )()
    module$eval()

    values <- torch::torch_randn(c(500L, 24L, 3L))
    expected <- module(values)

    # Slicing must not change the result, whatever the number of slices:
    # smaller than, equal to, larger than and an exact divisor of the input
    for (batch_size in c(64L, 128L, 250L, 500L, 1000L)) {
        bounded <- .torch_module_bounded(
            module = module, batch_size = batch_size
        )
        bounded$eval()
        output <- bounded(values)
        expect_equal(dim(output), dim(expected))
        expect_true(as.logical(torch::torch_allclose(output, expected)))
    }
})

test_that("bounded module handles blocks without valid pixels", {
    skip_on_cran()
    skip_if_not_installed("torch")

    module <- torch::nn_module(
        initialize = function() self$lin <- torch::nn_linear(3L, 5L),
        forward = function(x) self$lin(x)
    )()
    module$eval()
    bounded <- .torch_module_bounded(module = module, batch_size = 64L)
    bounded$eval()

    # Chunks fully covered by NA pixels reach the model with zero rows
    output <- bounded(torch::torch_randn(c(0L, 3L)))
    expect_equal(dim(output), c(0L, 5L))
})

test_that("bounded module concatenates list outputs", {
    skip_on_cran()
    skip_if_not_installed("torch")

    torch::torch_manual_seed(42L)
    # Encoders may return more than one tensor per forward pass
    module <- torch::nn_module(
        initialize = function() self$lin <- torch::nn_linear(3L, 2L),
        forward = function(x) list(self$lin(x), x * 2.0)
    )()
    module$eval()

    values <- torch::torch_randn(c(300L, 3L))
    expected <- module(values)
    bounded <- .torch_module_bounded(module = module, batch_size = 64L)
    bounded$eval()
    output <- bounded(values)

    expect_length(output, length(expected))
    expect_true(as.logical(torch::torch_allclose(output[[1L]], expected[[1L]])))
    expect_true(as.logical(torch::torch_allclose(output[[2L]], expected[[2L]])))
})

test_that("bounded module keeps the predict method of the model", {
    skip_on_cran()
    skip_if_not_installed("torch")

    torch::torch_manual_seed(42L)
    # Encoders such as sits_ssl_mae define a predict method which runs only
    # part of the model. luz calls it instead of forward, and so must the
    # bounded module
    module <- torch::nn_module(
        initialize = function() {
            self$encoder <- torch::nn_linear(3L, 2L)
            self$decoder <- torch::nn_linear(2L, 3L)
        },
        forward = function(x) self$decoder(self$encoder(x)),
        predict = function(x) self$encoder(x)
    )()
    module$eval()

    values <- torch::torch_randn(c(300L, 3L))
    bounded <- .torch_module_bounded(module = module, batch_size = 64L)
    bounded$eval()

    # predict() must return embeddings, not the reconstruction
    expect_true(as.logical(torch::torch_allclose(
        bounded$predict(values), module$predict(values)
    )))
    expect_true(as.logical(torch::torch_allclose(
        bounded(values), module(values)
    )))
})

test_that("batch_size is validated against the training batch size", {
    # The usual mistake is reusing the batch_size used to train a model
    expect_error(.check_batch_size(64L))
    expect_error(.check_batch_size(0L))
    expect_silent(.check_batch_size(2L^15L))
})

test_that("gpu_memory is deprecated", {
    expect_warning(.check_gpu_memory_deprecated(TRUE), "deprecated")
    expect_silent(.check_gpu_memory_deprecated(FALSE))
})
