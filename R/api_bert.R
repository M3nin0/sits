#' @title SITS-BERT noise-corruption time series dataset
#' @name .bert_dataset_lazy
#'
#' @description
#' Internal Torch dataset for SITS-BERT self-supervised pretraining. It is a
#' noise-corruption variant of \code{\link{.mae_dataset_lazy}}: instead of
#' overwriting a random subset of observations with a fixed \code{mask_value},
#' it adds zero-mean band-relative Gaussian noise to the selected positions.
#' The model is trained to recover the original (clean) values at the corrupted
#' positions, following the masked-observation pretext task of Yuan & Lin
#' (2021).
#'
#' Each item returns a corrupted input tensor \code{x} and a target structure
#' \code{y = list(y, mask)} where \code{y} is the clean (normalized) time
#' series and \code{mask} flags the corrupted positions used for the
#' masked-MSE loss.
#'
#' @param samples       A \code{sits} samples object.
#' @param indices       Integer indices into \code{samples} for this partition.
#' @param stats         Quantile statistics from \code{.samples_stats()} used
#'   for normalization.
#' @param bands         Character vector of band names.
#' @param timeline      Sample timeline (vector of dates).
#' @param mask_ratio    Numeric in (0, 1). Fraction of timesteps to corrupt.
#' @param masking_method Character. Position-selection strategy passed to
#'   \code{\link{.mae_mask_index}} (\code{"random"} or \code{"contiguous"}).
#' @param noise_frac    Numeric. Noise standard deviation as a fraction of each
#'   band's standard deviation (computed in the normalized domain).
#' @param noised_bands  Character vector of bands eligible for noise injection.
#' @param band_sd       Named numeric vector of per-band standard deviations
#'   (in the normalized domain), indexed by band name.
#'
#' @return A Torch dataset compatible with \code{torch::dataloader}.
#'
#' @author Felipe Carlos \email{efelipecarlos@@gmail.com}
#'
#' @references
#' Yuan, Y., & Lin, L. (2021). Self-Supervised Pretraining of Transformers for
#' Satellite Image Time Series Classification. IEEE JSTARS, 14, 474-487.
#'
#' @keywords internal
#' @noRd
.bert_dataset_lazy <- torch::dataset(
    name = ".BertDatasetLazy",

    initialize = function(samples,
                          indices,
                          stats,
                          bands,
                          timeline,
                          mask_ratio,
                          masking_method,
                          noise_frac,
                          noised_bands,
                          band_sd) {
        self$samples <- samples
        self$indices <- indices
        self$stats <- stats
        self$bands <- bands
        self$timeline <- timeline
        self$mask_ratio <- mask_ratio
        self$masking_method <- masking_method
        self$noise_frac <- noise_frac
        self$noised_bands <- noised_bands
        self$band_sd <- band_sd
    },

    .getitem = function(i) {
        idx <- self$indices[i]
        data <- self$samples[idx, , drop = FALSE]

        pred <- .predictors(data)
        pred <- .pred_normalize(pred, self$stats)

        ts <- .pred_as_ts(pred, self$bands, self$timeline)

        # original (clean) normalized data
        .ts(data) <- ts

        # select positions to corrupt
        # (reuses the MAE position picker)
        masked <- .mae_mask_index(
            n_samples = length(idx),
            timeline = self$timeline,
            mask_ratio = self$mask_ratio,
            masking_method = self$masking_method
        )

        # inject band-relative Gaussian noise at
        # the selected positions
        ts_noisy <- ts
        n_idx <- length(masked$masked_idx)

        for (b in self$noised_bands) {
            # standard deviation
            sd_b <- self$noise_frac * self$band_sd[[b]]

            # generate noise
            noise <- stats::rnorm(n_idx, mean = 0, sd = sd_b)

            # save ts and noise
            ts_noisy[masked$masked_idx, b] <-
                ts_noisy[masked$masked_idx, b] + noise
        }

        # define corrupted data
        data_noisy <- data

        # save corrupted time-series
        .ts(data_noisy) <- ts_noisy

        # if there is only one index, define a proper shape
        if (length(idx) == 1) {
            dim <- c(length(self$timeline), length(self$bands))
        }

        # otherwise, the shape is layers x time x bands
        else {
            dim <- c(length(idx), length(self$timeline), length(self$bands))
        }

        # transform inpute data as predictor
        x <- torch::torch_tensor(
            array(
                data = as.matrix(.pred_features(.predictors(data_noisy))),
                dim = dim
            ),
            dtype = torch::torch_float()
        )

        # define reference data as predictors
        y <- torch::torch_tensor(
            array(
                data = as.matrix(.pred_features(.predictors(data))),
                dim = dim
            ),
            dtype = torch::torch_float()
        )

        # define mask dimension
        mask_dim <- NULL

        if (length(idx) == 1) {
            mask_dim <- c(length(self$timeline), 1L)
        } else {
            mask_dim <- c(length(idx), length(self$timeline), 1L)
        }

        # transform mask to torch tense
        mask <- torch::torch_tensor(
            array(masked$mask, dim = mask_dim),
            dtype = torch::torch_float()
        )

        # return!
        list(
            x = x,
            y = list(
                y = y, mask = mask
            )
        )
    },
    .getbatch = function(i) {
        self$.getitem(i)
    },
    .length = function() {
        length(self$indices)
    }
)
