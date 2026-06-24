#' @title Train temporal pattern classification using TWDTW
#' @name sits_twdtw
#'
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description Use the Time-Weighted Dynamic Time Warping (TWDTW) algorithm
#' to classify satellite image time series. TWDTW compares each time series
#' against one temporal pattern per land-cover class, using a dynamic time
#' warping alignment with a temporal penalty that accounts for the phenological
#' timing of the classes. Each series is assigned the class of its closest
#' pattern.
#'
#' @note
#' TWDTW operates on the raw band values of both the time series and the
#' temporal patterns. No normalization is applied. The time weight depends on
#' the day-of-year of each observation, so it is most meaningful for data with
#' a yearly phenological cycle.
#'
#' @references Maus V, Camara G, Cartaxo R, Sanchez A, Ramos F, Queiroz GR.
#' A Time-Weighted Dynamic Time Warping Method for Land-Use
#' and Land-Cover Mapping. IEEE Journal of Selected Topics in Applied
#' Earth Observations and Remote Sensing, 9(8):3729-3739,
#' August 2016. ISSN 1939-1404. \doi{10.1109/JSTARS.2016.2517118}.
#'
#' Maus, V., Câmara, G., Appel, M., & Pebesma, E. (2019).
#' dtwSat: Time-Weighted Dynamic Time Warping for Satellite Image
#' Time Series Analysis in R. Journal of Statistical Software, 88(5), 1-31.
#' \doi{10.18637/jss.v088.i05}.
#'
#' @param samples        Time series with the training samples.
#' @param patterns       Optional reference time series to compare against
#'                       (tibble of class "sits", e.g. from
#'                       \code{\link[sits]{sits_patterns}}). If \code{NULL}
#'                       (default), the references are derived from
#'                       \code{samples} according to \code{pattern_method}.
#' @param pattern_method How references are derived from \code{samples} when
#'                       \code{patterns} is not provided: "none" (default)
#'                       keeps every sample and classifies by k-nearest
#'                       neighbours, as in \code{dtwSat::twdtw_knn1} with
#'                       \code{formula = NULL}; "gam" reduces each label to a
#'                       single GAM template (faster, but discards intra-class
#'                       variability).
#' @param k              Number of nearest references per class averaged into
#'                       the class distance (default: 1, i.e. nearest
#'                       neighbour). Only relevant when more than one reference
#'                       per class exists (\code{pattern_method = "none"}).
#' @param freq           Interval in days for the GAM pattern estimates
#'                       (\code{pattern_method = "gam"}; default: 8).
#' @param formula        Formula used by the GAM pattern estimate
#'                       (\code{pattern_method = "gam"}; default:
#'                       \code{y ~ s(x)}).
#' @param weight         Time-weight function: "logistic" (default), "linear"
#'                       or "none" (plain DTW).
#' @param alpha          Steepness of the logistic time weight (default: 0.1).
#'                       Larger time differences must yield larger penalties,
#'                       which requires a positive steepness with the logistic
#'                       form \eqn{1/(1+e^{-\alpha (g-\beta)})}.
#' @param beta           Midpoint, in days, of the logistic time weight
#'                       (default: 100).
#' @param a              Slope of the linear time weight (default: 0).
#' @param b              Intercept of the linear time weight (default: 0).
#' @param temperature    Softmax temperature used to convert TWDTW distances
#'                       into class probabilities (default: 1).
#' @param ...            Other parameters passed to
#'                       \code{\link[sits]{sits_patterns}}.
#' @return               Model fitted to input data
#'                      (to be passed to \code{\link[sits]{sits_classify}}).
#'
#' @examples
#' if (sits_run_examples()) {
#'     # Train a TWDTW model
#'     twdtw_model <- sits_train(samples_modis_ndvi, ml_method = sits_twdtw())
#'     # Classify a point
#'     point_ndvi <- sits_select(point_mt_6bands, bands = "NDVI")
#'     point_class <- sits_classify(
#'         data = point_ndvi, ml_model = twdtw_model
#'     )
#'     # Plot
#'     plot(point_class)
#' }
#' @export
sits_twdtw <- function(samples = NULL,
                       patterns = NULL,
                       pattern_method = "none",
                       k = 1L,
                       freq = 8L,
                       formula = y ~ s(x),
                       weight = "logistic",
                       alpha = 0.1,
                       beta = 100,
                       a = 0,
                       b = 0,
                       temperature = 1.0, ...) {
    .check_set_caller("sits_twdtw")
    # Function that trains a TWDTW model
    train_fun <- function(samples) {
        # Pre-conditions
        .check_chr_within(weight, c("logistic", "linear", "none"))
        .check_chr_within(pattern_method, c("none", "gam"))
        .check_int_parameter(k, min = 1L)
        .check_num_parameter(temperature, exclusive_min = 0.0)
        # Pre-condition: are predictors valid?
        .check_predictors(.predictors(samples), samples)
        # Get labels (used to fix column order in the result matrix)
        labels <- .samples_labels(samples)
        # Bands ordered as in the predictors layout (time-series bands only)
        bands <- .samples_bands(samples, include_base = FALSE)
        # Number of observations
        n_times <- .samples_ntimes(samples)
        # Day-of-year of the samples timeline (reference for the query series)
        query_doy <- .twdtw_doy(.samples_timeline(samples))
        # Define numeric code of the chosen time-weight function
        weight_type <- switch(
            weight,
            none = 0,
            logistic = 1,
            linear = 2
        )
        # Build (or validate) the reference series and convert to matrices.
        # By default every sample is a reference (k-NN). References carry a
        # per-reference label used to aggregate distances per class.
        refs <- .twdtw_references(
            samples = samples,
            method = pattern_method,
            freq = freq,
            formula = formula,
            patterns = patterns, ...
        )
        # Create reference matrix
        refs_matrix <- .twdtw_reference_matrices(
            references = refs,
            bands = bands
        )
        # Define distance power ("2" means euclidean)
        dist_power <- 2
        # Function that predicts results
        predict_fun <- function(values) {
            # Get number of input pixels
            input_pixels <- nrow(values)
            # Prepare values as features matrix
            values <- as.matrix(.pred_features(values))
            # Compute TWDTW distances
            distances <- C_twdtw_distances(
                values = values,
                patterns = refs_matrix[["matrices"]],
                query_doy = query_doy,
                pattern_doy = refs_matrix[["doys"]],
                n_bands = length(bands),
                dist_power = dist_power,
                weight_type = weight_type,
                alpha = alpha,
                beta = beta,
                a = a,
                b = b
            )
            # Aggregate per-reference distances into per-class distances (k-NN)
            distances <- .twdtw_class_distances(
                distances, refs_matrix[["labels"]], labels, k
            )
            # Convert distances into class probabilities
            values <- .twdtw_probs(distances, labels, temperature)
            # Are the results consistent with the data input?
            .check_processed_values(values, input_pixels)
            # Reorder matrix columns if needed
            if (any(labels != colnames(values))) {
                values <- values[, labels]
            }
            # Return!
            values
        }
        # Set model class
        predict_fun <- .set_class(
            predict_fun, "twdtw_model", "sits_model", class(predict_fun)
        )
        predict_fun
    }
    # If samples is informed, train a model and return a predict function
    # Otherwise give back a train function to train model further
    .factory_function(samples, train_fun)
}
