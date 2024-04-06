#' @title Detect changes in time-series using various methods.
#' @name .detect_change_ts
#' @keywords internal
#' @noRd
.detect_change_ts <- function(samples,
                              cd_method,
                              filter_fn,
                              multicores,
                              progress) {
    # Start parallel workers
    .parallel_start(workers = multicores)
    on.exit(.parallel_stop(), add = TRUE)
    # Get bands from model
    bands <- .ml_bands(cd_method)  # TODO: Should we use `ml_bands` ?
    # Update samples bands order
    if (any(bands != .samples_bands(samples))) {
        samples <- .samples_select_bands(samples = samples,
                                         bands = bands)
    }
    # Apply time series filter
    if (.has(filter_fn)) {
        samples <- .apply_across(data = samples,
                                 fn = filter_fn)
    }
    # Divide samples in chunks to parallel processing
    # TODO: Check generalization of the `.pred_create_partition` function
    parts <- .pred_create_partition(pred = samples, partitions = multicores)
    # Detect changes!
    detections <- .jobs_map_parallel_dfr(parts, function(part) {
        # TODO: Check generalization of the `.pred_part` function
        # Get samples
        values <- .pred_part(part)
        # Detect changes!
        detections <- cd_method(values[["time_series"]])
        detections <- tibble::tibble(detections)
        # Prepare and return results
        result <- tibble::tibble(data.frame(values, detections = detections))
        # TODO: Is this required ?
        class(result) <- class(values)
        result

    }, progress = progress)

    return(detections)
}
