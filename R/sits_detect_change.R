
#' @title Detect changes in time series
#' @name sits_detect_change
#' @description Given a set of time series or an image, this function compares each time series
#'                         with a set of change/no-change patterns, and indicates places and dates where
#'                          change has been detected
#' @param data         Set of time series
#' @param cd_method    Change detection method (with parameters)
#' @return             Set of time series where significant change has been detected
#' @export
sits_detect_change <- function(data, cd_method) {
    labels <- environment(cd_method)[["labels"]]

    slider::slide_dfr(data, function(data_row) {
        # Extract row time-series
        row_ts <- data_row$time_series[[1]]

        # Detect changes
        row_detections <-
            row_ts |> dplyr::select(-Index) |> as.matrix() |> cd_method()

        # Create results
        row_detections <-
            purrr::map_dfr(1:nrow(row_ts), function(ts_index) {
                row_index <- row_detections[[1]][[ts_index]]$index

                # Window defined by users can produce empty results. Validates
                # if this is the case for the current row.
                if (is.null(row_index)) {
                    return(data.frame())
                }

                # Extract indices of the `from` and `to` dates.
                row_detection_label <- NA
                row_detection_from <- row_index[[1]]
                row_detection_to <-
                    row_index[[length(row_index)]]

                # Define dates for the `from` and `to` dates.
                row_detection_from <-
                    row_ts[row_detection_from, ]$Index
                row_detection_to <- row_ts[row_detection_to, ]$Index

                # Extract distances calculated for the current time-series.
                distances <-
                    purrr::map(1:length(row_detections), function(row_detection_index) {
                        detection_distance <-
                            row_detections[[row_detection_index]][[ts_index]]$distance

                        ifelse(is.null(detection_distance),
                               NA,
                               detection_distance)
                    })

                # Get the min distance and use its associated pattern as the
                # change detected.
                distances_min <- which.min(distances)

                if (any(distances_min)) {
                    row_detection <-
                        row_detections[[distances_min]][[ts_index]]

                    if (!is.null(row_detection) &&
                        !is.na(distances[[distances_min]])) {
                        row_detection_label <- labels[distances_min]
                    }
                }

                # Create the data.frame with the structure:
                # `from`, `to`, `change`
                # `labels` are populated dynamically.
                stats::setNames(
                    data.frame(
                        from = row_detection_from,
                        to = row_detection_to,
                        change = row_detection_label,
                        distances
                    ),
                    c("from", "to", "change", labels)
                )
            })

        data_row$detections <- list(row_detections)
        data_row
    })
}
