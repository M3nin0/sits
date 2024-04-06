#' @title Dynamic Time Warping for Detect changes.
#' @name sits_dtw
#'
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#' @author Rolf Simoes, \email{rolf.simoes@@inpe.br}
#'
#' @description Create a Dynamic Time Warping (DTW) method for the
#' \code{\link[sits]{sits_detect_change_method}}.
#'
#' @param  samples            Time series with the training samples.
#' @param  threshold          Threshold used to define if an event was detected.
#' @param  start_date         Initial date of the interval used to calculate DTW.
#' @param  end_date           Final date of the interval used to calculate DTW.
#' @param  step               ISO8601-compliant time period used as step for the
#'                            DTW moving window, with number and unit,
#'                            where "D", "M" and "Y" stands for days, month and
#'                            year; e.g., "P16D" for 16 days.
#' @param  overlap            ISO8601-compliant time period of the overlap
#'                            between two moving window steps.
#' @return                  Change detection method prepared
#'                          to be passed to
#'                          \code{\link[sits]{sits_detect_change_method}}
#' @export
#'
sits_dtw <-
    function(samples    = NULL,
             threshold  = 2,
             start_date = NULL,
             end_date   = NULL,
             window     = NULL,
             step       = NULL) {
        train_fun <-
            function(samples) {
                # TODO: Validate `start_date` and `end_date` (both can be null)
                # TODO: Validate `window` (can't be null)
                # TODO: Validate `step` (can't be null)
                .period_check(period = window)
                .period_check(period = step)
                # TODO: Validate if `window (period unit)` = `step (period unit)` ?

                # Sample labels
                labels <- .samples_labels(samples)
                # Get samples patterns (temporal median)
                train_samples_patterns <-
                    samples |> sits_select(start_date = start_date, end_date = end_date) |>
                    dplyr::group_by(label) |>
                    dplyr::group_map(function(data, name) {
                        dplyr::bind_rows(data[["time_series"]]) |>
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
                    # TODO: Validate values time-line (should be inside period + step)
                    # Extract values
                    min_date <- min(values[[1]][["Index"]])
                    max_date <- max(values[[1]][["Index"]])
                    # Create comparison windows
                    comparison_windows <- .period_windows(
                        period = window,
                        step = step,
                        start_date = min_date,
                        end_date = max_date
                    )
                    # Do the change detection for each time-series
                    purrr::map(values, function(value_row) {
                        # Search for the patterns
                        patterns_distances <- purrr::map_dfc(1:length(train_samples_patterns), function(pattern_index) {
                            # Get pattern metadata
                            pattern_ts <-
                                train_samples_patterns[pattern_index]
                            pattern_ts <- pattern_ts[[1]]
                            pattern_label <- labels[pattern_index]
                            # Calculate distances in each comparison window
                            distances <- purrr::map_df(comparison_windows, function(comparison_window) {
                                # Get time-series in the window
                                value_row_in_window <-
                                    dplyr::filter(value_row,
                                                  .data[["Index"]] >= comparison_window[["start"]],
                                                  .data[["Index"]] <= comparison_window[["end"]])
                                # Remove the time reference column
                                value_row_in_window <-
                                    dplyr::select(value_row_in_window, -.data[["Index"]])
                                # Transform values in matrix for the cpp code
                                # TODO: Validate order of bands
                                value_row_in_window <-
                                    as.matrix(value_row_in_window)
                                # DTW distance between the pattern and the
                                # current window time-series.
                                distance_from_pattern <-
                                    distance_dtw(value_row_in_window, as.matrix(pattern_ts))
                                # Create result object
                                data.frame(distance = distance_from_pattern)
                            })
                            # Associate the pattern name with the distance and
                            # avoid extra manipulations in the next steps
                            stats::setNames(distances, pattern_label)
                        })
                        # Remove distances out the user-defined threshold
                        patterns_distances[patterns_distances > threshold] <- NA
                        # Select the min value, not null, as the change detected
                        patterns_distances[["change"]] <- purrr::pmap(patterns_distances, function(...) {
                            row <- data.frame(...)

                            min_idx <- which.min(row)
                            min_val <- row[min_idx]

                            ifelse(is.null(min_val), NA, names(min_val))
                        })

                        # Result associating the window dates and the changes
                        # detected
                        data.frame(
                            # Date columns
                            stats::setNames(
                                purrr::map_df(comparison_windows, unlist),
                                            c("from", "to")
                            ),
                            # changes
                            change = unlist(patterns_distances[["change"]]),
                            # distances
                            dplyr::select(patterns_distances, -.data[["change"]])
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
