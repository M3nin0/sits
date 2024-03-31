

#' @title Train dtw
#' @name sits_dtw
#' @export
#'
sits_dtw <-
    function(samples = NULL,
             distance_threshold,
             window_size_left = 2,
             window_size_right = 2,
             window_step = 4,
             window_complete = T) {
        train_fun <-
            function(samples) {
                # TODO: Temporary solution
                train_samples <- .predictors(samples)

                labels <- .samples_labels(samples)
                # Get bands available in the samples
                train_samples_bands <- paste0("(",
                                              paste(sits_bands(samples), collapse = "|"),
                                              ")")
                # Get predictors names
                train_samples_bands_preds <- names(train_samples)
                train_samples_bands_preds <-
                    train_samples_bands_preds[3:length(train_samples_bands_preds)]
                # Get samples patterns (temporal median)
                train_samples_patterns <-
                    dplyr::group_by(samples, label) |>
                    dplyr::group_map(function(data, name) {
                        dplyr::bind_rows(data$time_series) |>
                            dplyr::group_by(Index) |>
                            dplyr::summarize(dplyr::across(dplyr::everything(),
                                                           stats::median, na.rm = TRUE)) |>
                            dplyr::select(-Index) |>
                            as.matrix()
                    })

                # TODO values = 1 time-series as matrix (with multiple columns)
                # TODO review the definition of values
                detect_change_fun <- function(values) {
                    # Search for each pattern in the `values` time-series
                    purrr::map(1:length(train_samples_patterns), function(pattern_index) {
                        pattern_name <- labels[pattern_index]

                        pattern_ts <-
                            train_samples_patterns[pattern_index]
                        pattern_ts <- pattern_ts[[1]]

                        slider::slide(
                            1:length(values),
                            .f = function(index_in_window) {
                                values_in_window <- values[index_in_window, ]

                                ts_window <-
                                    as.matrix(values_in_window)

                                distance_from_pattern <-
                                    distance_dtw(ts_window, pattern_ts)

                                if (distance_from_pattern <= distance_threshold) {
                                    return(
                                        list(
                                            index = index_in_window,
                                            change = pattern_name,
                                            distance = distance_from_pattern
                                        )
                                    )
                                }

                                list(
                                    index = index_in_window,
                                    change = NULL,
                                    distance = NULL
                                )
                            },
                            .before   = window_size_left,
                            .after    = window_size_right,
                            .step     = window_step,
                            .complete = window_complete
                        )
                    })
                }

                # Set model class
                detect_change_fun <- .set_class(detect_change_fun,
                                                "dtw_model",
                                                "sits_model",
                                                class(detect_change_fun))
                return(detect_change_fun)
            }
        # If samples is informed, train a model and return a predict function
        # Otherwise give back a train function to train model further
        result <- .factory_function(samples, train_fun)
        return(result)
    }

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

                if (is.null(row_index)) {
                    return(data.frame())
                }

                row_detection_label <- NA
                row_detection_from <- row_index[[1]]
                row_detection_to <-
                    row_index[[length(row_index)]]

                row_detection_from <-
                    row_ts[row_detection_from, ]$Index
                row_detection_to <- row_ts[row_detection_to, ]$Index

                distances <-
                    purrr::map(1:length(row_detections), function(row_detection_index) {
                        detection_distance <-
                            row_detections[[row_detection_index]][[ts_index]]$distance

                        ifelse(is.null(detection_distance),
                               NA,
                               detection_distance)
                    })

                distances_min <- which.min(distances)

                if (any(distances_min)) {
                    row_detection <-
                        row_detections[[distances_min]][[ts_index]]

                    if (!is.null(row_detection) &&
                        !is.na(distances[[distances_min]])) {
                        row_detection_label <- labels[distances_min]
                    }
                }

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
