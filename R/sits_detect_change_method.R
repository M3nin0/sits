#' @title Create detect change method.
#' @name sits_detect_change_method
#'
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#'
#' @description Prepare detection change method. Currently, sits supports the
#' following methods: 'dtw' (see \code{\link[sits]{sits_dtw}})
#'
#' @param  samples          Time series with the training samples.
#' @param  cd_method        Change detection method.
#' @return                  Change detection method prepared
#'                          to be passed to
#'                          \code{\link[sits]{sits_detect_change}}
#' @export
#'
sits_detect_change_method <-
    function(samples, cd_method = sits_dtw()) {
        # set caller to show in errors
        .check_set_caller("sits_detect_change_method")
        # check if samples are valid
        .check_samples_train(samples)

        # is the train method a function?
        .check_that(x = inherits(cd_method, "function"),
                    msg = "cd_method is not a valid function")
        # are the timelines OK?
        .check_that(
            x = .timeline_check(samples) == TRUE,
            msg = paste0(
                "Samples have different timeline lengths",
                "\n",
                "Use .tibble_prune or sits_fix_timeline"
            )
        )
        # compute the training method by the given data
        result <- cd_method(samples)
        # return a valid machine learning method
        return(result)
    }


#' @title Dynamic Time Warping for Detect changes.
#' @name sits_dtw
#'
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#'
#' @description Create a Dynamic Time Warping (DTW) method for the
#' \code{\link[sits]{sits_detect_change_method}}.
#'
#' @param  samples            Time series with the training samples.
#' @param  distance_threshold Threshold used to define if an event was detected.
#' @param  window_size_left   Size of window in the left side of the current
#'                            position.
#' @param  window_size_right  Size of window in the right side of the current
#'                            position.
#' @param  window_step        Window step size.
#' @param  window_complete    Flag indicating if only complete window, including
#'                            all values defined in `window_size_left` and
#'                            `window_size_right`, should be used.
#' @return                  Change detection method prepared
#'                          to be passed to
#'                          \code{\link[sits]{sits_detect_change_method}}
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
                # Sample labels
                labels <- .samples_labels(samples)
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

                # TODO review the definition of values
                # The current definition of it is:
                #   1 time-series as matrix (with multiple columns)
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
                                values_in_window <- values[index_in_window,]

                                ts_window <-
                                    as.matrix(values_in_window)

                                # DTW distance between the pattern and the current
                                # window time-series.
                                distance_from_pattern <-
                                    distance_dtw(ts_window, pattern_ts)

                                # Validate if the pattern is under the user-defined
                                # threshold. If yes, it is a valid detection.
                                # Otherwise, just ignore it.
                                if (distance_from_pattern <= distance_threshold) {
                                    return(
                                        list(
                                            index = index_in_window,
                                            change = pattern_name,
                                            distance = distance_from_pattern
                                        )
                                    )
                                }

                                # Empty list is returned to keep the consistency
                                # of the results.
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
