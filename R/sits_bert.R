#' @title Pre-train a SITS-BERT encoder on SITS time-series data
#'
#' @name sits_bert
#'
#' @description
#' \code{sits_bert()} creates a SITS-BERT pretraining factory compatible with
#' \code{\link{sits_pre_train}}. It implements the self-supervised
#' masked-observation pretext task of Yuan & Lin (2021): a random subset of
#' observations in each time series is corrupted with zero-mean, band-relative
#' Gaussian noise, and a transformer encoder is trained to recover the original
#' (clean) values at the corrupted positions. After pretraining the
#' reconstruction decoder is discarded and the pretrained encoder is returned
#' for downstream use via \code{\link{sits_encode}}.
#'
#' The noise standard deviation for each band is \code{noise_frac} times that
#' band's standard deviation, measured in the normalized domain. This makes the
#' corruption magnitude comparable across bands of different dynamic range.
#'
#' The function can be used in two ways:
#' \itemize{
#' \item If \code{samples} is provided, it trains immediately and returns
#' an encoder-ready model object (see Value).
#' \item If \code{samples = NULL}, it returns a training function with
#' signature \code{function(samples)} that can be passed to
#' \code{\link{sits_pre_train}} or called later.
#' }
#'
#' @param samples A \code{sits} samples object. If \code{NULL} (default),
#'   returns a training function. If provided, triggers immediate training.
#'   Base data samples (e.g., \code{sits_base}) are not supported.
#' @param embedding_dim Integer. Dimensionality of the latent embedding
#'   produced by the encoder (Default: 32L).
#' @param encoder_model Function or encoder factory. Defines the encoder
#'   backbone to be instantiated for SITS-BERT pretraining. SITS-BERT is a
#'   transformer method, so an attention-based backbone is expected
#'   (Default: \code{\link{sits_lighttae}()}). Must accept \code{samples} and
#'   \code{embedding_dim} and return a \code{torch::nn_module}.
#' @param decoder_width Integer. Width of the decoder MLP hidden layer.
#' @param masking_method Character. Strategy used to select the observations to
#'   corrupt. Either \code{"random"} (default) or \code{"contiguous"}.
#' @param mask_ratio Numeric in (0, 1). Fraction of timesteps to corrupt with
#'   noise. Default: 0.6.
#' @param noise_frac Numeric (>= 0). Noise standard deviation as a fraction of
#'   each band's standard deviation (in the normalized domain). Default: 0.5.
#' @param noised_bands Character vector specifying which bands are eligible for
#'   noise injection. If \code{NULL} (default), all bands are eligible.
#' @param epochs Integer. Maximum number of training epochs.
#' @param batch_size Integer. Batch size used for training and validation.
#' @param validation_split Numeric in (0, 1). Fraction of samples held out
#'   for validation loss monitoring.
#' @param optimizer Function. A \code{torch} optimizer constructor, such as
#'   \code{torch::optim_adamw}.
#' @param opt_hparams List of optimizer hyperparameters passed to
#'   \code{optimizer}. Common entries include \code{lr}, \code{eps}, and
#'   \code{weight_decay}. Only parameters supported by the chosen optimizer
#'   are accepted.
#' @param lr_decay_epochs Integer. Step size (in epochs) for learning-rate
#'   decay when using the step scheduler.
#' @param lr_decay_rate Numeric. Multiplicative decay factor applied by the
#'   learning-rate scheduler.
#' @param patience Integer. Number of epochs without improvement in
#'   validation loss before early stopping.
#' @param min_delta Numeric. Minimum decrease in validation loss required
#'   to reset the early-stopping patience counter.
#' @param verbose Logical. If \code{TRUE}, prints training progress and
#'   per-epoch losses.
#' @param seed Integer. Random seed used to initialize Torch randomness.
#'
#' @details
#' During training, the procedure:
#' \enumerate{
#' \item Normalizes inputs using quantile-based statistics derived from the
#' samples and computes per-band standard deviations in the normalized domain.
#' \item Selects observations to corrupt according to \code{masking_method} and
#' \code{mask_ratio}.
#' \item Adds zero-mean band-relative Gaussian noise (scaled by
#' \code{noise_frac}) to the selected observations on \code{noised_bands}.
#' \item Trains a transformer encoder-decoder with a masked MSE objective,
#' where loss is computed only over corrupted positions.
#' \item Applies early stopping and step learning-rate scheduling.
#' \item Discards the decoder after training, returning the pretrained encoder.
#' }
#'
#' \strong{Scope.} This is the feature-extractor variant of SITS-BERT: the
#' pretrained encoder is frozen and used to produce embeddings for a downstream
#' classifier (\code{\link{sits_encode}} followed by \code{\link{sits_train}}).
#' It does not perform end-to-end fine-tuning of the encoder together with a
#' classification head.
#'
#' @return
#' If \code{samples = NULL}, returns a training function with signature
#' \code{function(samples)} that trains a SITS-BERT encoder and returns a
#' pretrained encoder (a \code{sits_encoder} closure).
#'
#' If \code{samples} is provided, returns the result of applying the training
#' function to \code{samples}.
#'
#' @references
#' Yuan, Y., & Lin, L. (2021). Self-Supervised Pretraining of Transformers for
#' Satellite Image Time Series Classification. \emph{IEEE Journal of Selected
#' Topics in Applied Earth Observations and Remote Sensing}, 14, 474-487.
#' \doi{10.1109/JSTARS.2020.3036602}
#'
#' @author Felipe Carlos \email{efelipecarlos@@gmail.com}
#'
#' @examples
#' if (sits_run_examples()) {
#'     bert_model <- sits_pre_train(
#'         samples = samples_modis_ndvi,
#'         encoder_method = sits_bert(
#'             encoder_model = sits_lighttae(),
#'             noise_frac = 0.5
#'         )
#'     )
#'
#'     encoded <- sits_encode(samples_modis_ndvi, bert_model)
#' }
#'
#' @export
sits_bert <- function(samples = NULL,
                      embedding_dim = 32L,
                      encoder_model = sits_lighttae(),
                      decoder_width = 128L,
                      masking_method = "random",
                      mask_ratio = 0.6,
                      noise_frac = 0.5,
                      noised_bands = NULL,
                      epochs = 150L,
                      batch_size = 128L,
                      validation_split = 0.2,
                      optimizer = torch::optim_adamw,
                      opt_hparams = list(
                          lr           = 5.0e-04,
                          eps          = 1.0e-08,
                          weight_decay = 1.0e-06
                      ),
                      lr_decay_epochs = 1,
                      lr_decay_rate = 0.95,
                      patience = 20,
                      min_delta = 0.01,
                      verbose = FALSE,
                      seed = 10L) {
    # set caller for error msg
    .check_set_caller("sits_bert")
    # Verifies if 'torch' and 'luz' packages is installed
    .check_require_packages(c("torch", "luz"))
    # Documentation mode? verbose is FALSE
    verbose <- .message_verbose(verbose)
    # Band prefix for embeddings
    bands_prefix = .conf("embedding_band_prefix")
    # Check bands prefix
    .check_chr(bands_prefix, len_min = 1, lan_max = 1, allow_empty = FALSE)
    # Function that trains a torch model based on samples
    train_fun <- function(samples) {
        # does not support working with DEM or other base data
        if (inherits(samples, "sits_base")) {
            stop(.conf("messages", "sits_train_base_data"), call. = FALSE)
        }
        # Avoid add a global variable for 'self' and 'super'
        self <- NULL
        super <- NULL
        # Pre-conditions
        .check_pre_sits_bert(
            samples = samples,
            epochs = epochs,
            batch_size = batch_size,
            encoder = encoder_model,
            decoder_width = decoder_width,
            masking_method = masking_method,
            mask_ratio = mask_ratio,
            noise_frac = noise_frac,
            noised_bands = noised_bands,
            bands_prefix = bands_prefix,
            verbose = verbose
        )
        # Other pre-conditions
        .check_int_parameter(seed, allow_null = TRUE)
        # Check opt_hparams
        # Get parameters list and remove the 'param' parameter
        optim_params_function <- formals(optimizer)[-1L]
        # Check hparams
        .check_opt_hparams(opt_hparams, optim_params_function)
        # Prepare params
        optim_params_function <- utils::modifyList(
            x = optim_params_function,
            val = opt_hparams
        )
        # Samples labels
        labels <- .samples_labels(samples)
        # Samples bands
        bands <- .samples_bands(samples)
        # Samples timeline
        timeline <- .samples_timeline(samples)
        # Number of labels, bands, and number of samples (used below)
        n_labels <- length(labels)
        n_bands <- length(bands)
        n_times <- length(timeline)
        # Copy embedding_dim from parent environment to local
        embedding_dim <- embedding_dim
        # Copy bands_prefix from parent environment to local
        bands_prefix <- bands_prefix
        # If not has noised_bands set to all bands
        if (!.has(noised_bands)) {
            noised_bands <- bands
        }
        # Process samples for bert training
        ml_stats <- .samples_stats(samples)
        # Per-band standard deviation in the normalized domain. Used to scale
        # the band-relative Gaussian noise of the masked-observation pretext.
        norm_ts <- .pred_as_ts(
            data = .pred_normalize(.predictors(samples), ml_stats),
            bands = bands,
            timeline = timeline
        )
        # Generate standard deviation by band
        band_sd <- vapply(
            bands,
            function(b) stats::sd(norm_ts[[b]], na.rm = TRUE),
            numeric(1)
        )
        # Set band names
        names(band_sd) <- bands
        # Split train and validation
        idx <- sample.int(nrow(samples))
        # Define number of validation samples
        n_val <- floor(length(idx) * validation_split)
        # If there are validation samples, use it
        if (n_val > 0L) {
            val_idx <- idx[seq_len(n_val)]
            train_idx <- idx[-seq_len(n_val)]
        }
        # Otherwise, set as empty sample set
        else {
            val_idx <- integer(0)
            train_idx <- idx
        }
        # Torch dataset and dataloaders (datasets lazy)
        train_ds <- .bert_dataset_lazy(
            samples = samples,
            indices = train_idx,
            stats = ml_stats,
            bands = bands,
            timeline = timeline,
            mask_ratio = mask_ratio,
            masking_method = masking_method,
            noise_frac = noise_frac,
            noised_bands = noised_bands,
            band_sd = band_sd
        )
        # Validation dataset
        val_ds <- .bert_dataset_lazy(
            samples = samples,
            indices = val_idx,
            stats = ml_stats,
            bands = bands,
            timeline = timeline,
            mask_ratio = mask_ratio,
            masking_method = masking_method,
            noise_frac = noise_frac,
            noised_bands = noised_bands,
            band_sd = band_sd
        )
        # Create a torch seed (we define a new variable to allow users
        # to access this seed number from the model environment)
        torch_seed <- .torch_seed(seed)
        # Set torch seed
        torch::torch_manual_seed(torch_seed)
        # Set the encoder model closure
        encoder <- encoder_model(
            samples = samples,
            embedding_dim = embedding_dim
        )
        # Set the decoder model closure (reused from the MAE workflow)
        decoder <- .sits_mae_decoder_mlp(
            embedding_dim = embedding_dim,
            decoder_width = decoder_width,
            n_times = n_times,
            n_bands = n_bands
        )
        # Define full SITS-BERT pretraining model (encoder + reconstruction)
        bert_model <- torch::nn_module(
            classname = "BERT_model",
            initialize = function(encoder,
                                  decoder,
                                  n_bands = NULL,
                                  n_labels = NULL,
                                  timeline = NULL) {
                super$initialize()
                self$encoder <- encoder
                self$decoder <- decoder

                # keep metadata around for safety
                self$n_bands <- n_bands
                self$n_labels <- n_labels
                self$timeline <- timeline
            },
            forward = function(x) {
                x <- self$encoder(x)
                x <- self$decoder(x)
                torch::nnf_sigmoid(x)
            },
            predict = function(x) {
                x <- self$encoder(x)
                torch::nnf_sigmoid(x)
            }
        )
        # Loss function (masked MSE over corrupted positions only)
        bert_loss <- function(pred, target) {
            if (!is.null(target$y)) {
                y_true <- target$y
                mask <- target$mask
            } else {
                y_true <- target[[1]]
                mask <- target[[2]]
            }
            # Calc loss numerator
            num <- torch::nnf_mse_loss(
                pred * mask,
                y_true * mask,
                reduction = "sum"
            )
            # avoid divide-by-zero edge cases
            denom <- mask$sum()$clamp_min(1)
            # return loss
            num / denom
        }
        # Verify if GPU is available
        cpu_train <- .torch_cpu_train()
        # Train the model using luz
        torch_model <-
            luz::setup(
                module = bert_model,
                loss = bert_loss,
                optimizer = optimizer
            ) |>
            luz::set_hparams(
                encoder  = encoder,
                decoder  = decoder,
                n_bands  = n_bands,
                n_labels = n_labels,
                timeline = timeline
            ) |>
            luz::set_opt_hparams(
                !!!optim_params_function
            ) |>
            luz::fit(
                data = train_ds,
                epochs = epochs,
                valid_data = val_ds,
                callbacks = list(
                    luz::luz_callback_early_stopping(
                        monitor = "valid_loss",
                        mode = "min",
                        patience = patience,
                        min_delta = min_delta
                    ),
                    luz::luz_callback_lr_scheduler(
                        torch::lr_step,
                        step_size = lr_decay_epochs,
                        gamma = lr_decay_rate
                    )
                ),
                accelerator = luz::accelerator(cpu = cpu_train),
                dataloader_options = list(
                    batch_size = batch_size,
                    shuffle = TRUE
                ),
                verbose = verbose
            )

        # Serialize model
        serialized_model <- force(.torch_serialize_model(torch_model$model))

        # Function that encodes input values using the trained encoder
        predict_fun <- function(values) {
            # Verifies if torch package is installed
            .check_require_packages("torch")
            # Set torch threads to 1
            suppressWarnings(torch::torch_set_num_threads(1L))
            # Unserialize model
            torch_model$model <- .torch_unserialize_model(
                model = torch_model$model,
                raw = serialized_model
            )
            # Transform input into a 3D tensor
            # Reshape the 2D matrix into a 3D array
            n_samples <- nrow(values)
            # Performs data normalization
            values <- .pred_normalize(pred = values, stats = ml_stats)
            # Represent matrix values as array
            values <- array(
                data = as.matrix(values), dim = c(n_samples, n_times, n_bands)
            )
            # GPU or CPU classification?
            if (.torch_gpu_classification()) {
                # Get batch size
                batch_size <- sits_env[["batch_size"]]
                # Transform the input array to a dataset
                values <- .torch_as_dataset(values)
                # Transform to dataloader to use the batch size
                values <- torch::dataloader(values, batch_size = batch_size)
                # Do GPU classification
                values <- .try(
                    stats::predict(object = torch_model, values),
                    .msg_error = .conf("messages", ".check_gpu_memory_size")
                )
            } else {
                # Do CPU classification
                values <- stats::predict(object = torch_model, values)
            }
            # Convert from tensor to array
            values <- torch::as_array(values)
            # Update the columns names to labels
            colnames(values) <- paste0(bands_prefix, seq_len(ncol(values)))
            values
        }
        # Set model class
        predict_fun <- .set_class(
            predict_fun, "torch_model", "sits_encoder", class(predict_fun)
        )
    }
    # If samples is informed, train a model and return a predict function
    # Otherwise give back a train function to train model further
    .factory_function(samples, train_fun)
}
