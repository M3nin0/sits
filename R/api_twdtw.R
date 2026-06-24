#' @title Build (or validate) the reference time series for TWDTW
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @keywords internal
#' @noRd
#' @description
#' Returns the set of labelled reference time series the query series are
#' compared against. By default (\code{method = "none"}) every training
#' sample is kept as a reference, so prediction is a k-nearest-neighbour
#' search over all samples (as in \code{dtwSat::twdtw_knn1} with
#' \code{formula = NULL}). With \code{method = "gam"} the samples of each
#' label are reduced to a single GAM template (one reference per class),
#' which is faster but discards intra-class variability.
#'
#' @param samples  Training samples.
#' @param method   "none" (all samples) or "gam" (one template per label).
#' @param freq     Interval in days for the GAM estimates (method "gam").
#' @param formula  Formula used by the GAM estimate (method "gam").
#' @param patterns Optional user-supplied references.
#' @param ...      Additional parameters passed to \code{sits_patterns}.
#' @return         References tibble.
.twdtw_references <- function(samples,
                              method,
                              freq,
                              formula,
                              patterns = NULL,
                              ...) {
    # Samples labels
    labels <- .samples_labels(samples)
    # User-supplied references take precedence, so if user
    # defined patterns, use them
    if (.has(patterns)) {
        # Pre-condition: Check if the pattern is a valid sits tibble
        .check_that(inherits(patterns, "sits"))
        # Pre-condition: All training labels must have a corresponding reference
        .check_that(all(labels %in% .samples_labels(patterns)))
        # Pre-condition: References must contain the training bands
        .check_that(all(.samples_bands(samples) %in% .samples_bands(patterns)))
        # Return!
        return(patterns)
    }
    # k-NN over all samples (no reduction)
    if (method == "none") {
        return(samples)
    }
    # One GAM template per label
    .check_that(method == "gam")
    # Generate and return!
    sits_patterns(
        data = samples,
        freq = freq,
        formula = formula, ...
    )
}

#' @title Convert reference time series into matrices (one per reference)
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @keywords internal
#' @noRd
#' @param references References tibble, one row per reference.
#' @param bands      Bands giving the desired column order.
#' @return           List with \code{matrices} (one numeric matrix per
#'                   reference, time x bands), \code{doys} (day-of-year vector
#'                   per reference) and \code{labels} (label of each reference).
.twdtw_reference_matrices <- function(references, bands) {
    # Get time-series
    ts_list <- references[["time_series"]]
    # Order values by bands
    matrices <- purrr::map(ts_list, function(ts) {
        as.matrix(ts[, bands, drop = FALSE])
    })
    # Generate day-of-year
    doys <- purrr::map(ts_list, function(ts) {
        .twdtw_doy(ts[["Index"]])
    })
    # Prepare result object
    list(
        matrices = matrices,
        doys = doys,
        labels = references[["label"]]
    )
}

#' @title Day-of-year of a vector of dates
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @keywords internal
#' @noRd
#' @param dates A vector of dates.
#' @return      Numeric vector with the day-of-year of each date.
.twdtw_doy <- function(dates) {
    as.numeric(lubridate::yday(lubridate::as_date(dates)))
}

#' @title Aggregate per-reference distances into per-class distances (k-NN)
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @keywords internal
#' @noRd
#' @description
#' For each query series and each class, summarizes the TWDTW distances to the
#' references of that class by the mean of the \code{k} smallest distances
#' (\code{k = 1} gives the nearest-neighbour distance). With one reference per
#' class this reduces to that reference's distance. The class with the smallest
#' aggregated distance is the k-NN prediction.
#'
#' @param distances  Matrix (samples x references) of TWDTW distances.
#' @param ref_labels Label of each reference (length = number of columns).
#' @param labels     Class labels giving the column order of the result.
#' @param k          Number of nearest references per class to average.
#'
#' @return           Matrix (samples x classes) of aggregated distances.
.twdtw_class_distances <- function(distances, ref_labels, labels, k) {
    # Define class matrix
    class_dist <- matrix(
        NA_real_,
        nrow = nrow(distances),
        ncol = length(labels),
        dimnames = list(NULL, labels)
    )
    # Process by label
    for (j in seq_along(labels)) {
        # Get columns
        cols <- which(ref_labels == labels[[j]])
        # Get distances of given label
        cols_dist <- distances[, cols, drop = FALSE]
        # Define distances seq
        cols_dist_seq <- seq_len(nrow(cols_dist))
        # Define number of neighbors
        neighbors <- min(k, ncol(cols_dist))
        # Define class value
        if (neighbors == 1L) {
            # Nearest-neighbour distance to this class
            class_val <- purrr::map_dbl(cols_dist_seq, function(i) {
                min(cols_dist[i, ], na.rm = TRUE)
            })
        } else {
            # Mean of the k smallest distances to this class
            class_val <- purrr::map_dbl(cols_dist_seq, function(i) {
                mean(sort(cols_dist[i, ])[seq_len(neighbors)])
            })
        }
        # Save in the results matrix
        class_dist[, j] <- class_val
    }
    # Return!
    class_dist
}

#' @title Convert TWDTW distances into class probabilities
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @keywords internal
#' @noRd
#'
#' @param distances   Matrix of TWDTW distances.
#' @param labels      Class labels (column names of the result).
#' @param temperature Softmax temperature (smaller = harder assignment).
#' @return            Matrix of probabilities that sum to one per row, with
#'                    columns named after \code{labels}.
#'
#' @note The class chosen by \code{which.max} of this softmax is exactly the
#'   minimum-distance (k-NN) class for any \code{temperature > 0}. Temperature
#'   only controls how peaked the probabilities are, not the predicted label.
.twdtw_probs <- function(distances, labels, temperature) {
    # Define logits
    # Subtract the row-wise max of (-distance / temperature) for numerical
    # stability (shift-invariance of softmax).
    logits <- -distances / temperature
    logits <- logits - apply(logits, 1L, max)
    # Log for probabilities
    probs <- exp(logits)
    probs <- probs / rowSums(probs)
    # Name probs
    colnames(probs) <- labels
    # Softmax over negative distances:
    #  > Closer reference means higher probability.
    probs
}
