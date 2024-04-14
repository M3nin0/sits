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
#' @param  start_date         Initial date of the interval used to calculate DTW.
#' @param  end_date           Final date of the interval used to calculate DTW.
#' @param  step               ISO8601-compliant time period used as step for the
#'                            DTW moving window, with number and unit,
#'                            where "D", "M" and "Y" stands for days, month and
#'                            year; e.g., "P16D" for 16 days.
#' @param  threshold          Distance threshold used to define if an event was
#'                            detected.
#' @param  ...                Parameters to configure the
#'                            \code{\link[sits_som]{sits_detect_change_method}}
#'                            used to create the temporal patterns of samples.
#' @return                  Change detection method prepared
#'                          to be passed to
#'                          \code{\link[sits]{sits_detect_change_method}}
#' @export
#'
sits_dtw <-
    function(samples       = NULL,
             start_date    = NULL,
             end_date      = NULL,
             window        = NULL,
             step          = NULL,
             threshold     = 2,
             ...) {
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
                timeline <- sits_timeline(samples)
                # Create samples patterns (SOM)
                samples_som_grid <- sits_som_map(
                    data      = samples,
                    ...
                )
                # Extract properties from SOM to create the samples patterns
                som_properties <- samples_som_grid[["som_properties"]]
                # Codebook labels
                som_labels <- som_properties[["neuron_label"]]
                som_labels_names <- unique(som_labels)
                # Codebook attributes (e.g., bands, spectral indices)
                som_codes <- som_properties[["codes"]]
                som_codes_attributes <- names(som_codes)
                som_codes_attributes <- unique(som_codes_attributes)
                # Create patterns
                train_samples_patterns <- purrr::map_df(som_labels_names,
                                                         function(label) {
                     # Define which rows are associated with the current label
                     label_rows <- which(som_labels == label, TRUE)
                     # Create a table with the label rows and their respective
                     # data attributes
                     purrr::map_df(label_rows, function(label_row) {
                         # Define the label of the current row
                         row_label <- som_labels[label_row]
                         # Extract attributes data associated with the current row
                         row_data <- purrr::map_dfc(som_codes_attributes,
                                                    function(attribute) {
                             som_code_data <- som_codes[[attribute]]
                             stats::setNames(
                                 data.frame(som_code_data[label_row,]),
                                 attribute
                             )
                         })
                         row_data <- tibble::tibble(
                             data.frame(Index = timeline, row_data)
                         )
                         row_data <- list(row_data)
                         # Prepare the results following the organization
                         row_result <- data.frame(label = row_label)
                         row_result[["patterns"]] <- row_data
                         row_result
                     })
                })
                # Define patterns name to facilitate the detection phase (Below)
                patterns_names <- unique(train_samples_patterns[["label"]])
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
                        patterns_distances <- purrr::map_dfc(patterns_names, function(pattern_name) {
                            # Get pattern object
                            pattern_obj <- dplyr::filter(train_samples_patterns,
                                                         .data[["label"]] == pattern_name)
                            # Calculate the distance between each pattern
                            # behavior and the input data
                            pattern_distances <- purrr::map_dfc(seq_len(nrow(pattern_obj)), function(row_index) {
                                pattern_row <- pattern_obj[row_index,]
                                # Get time-series
                                pattern_row_ts <- pattern_row[["patterns"]][[1]]
                                # Filter in the user-defined interval
                                pattern_row_ts <- dplyr::filter(
                                    pattern_row_ts,
                                    .data[["Index"]] >= start_date,
                                    .data[["Index"]] <= end_date
                                )
                                # Get data as matrix to use the DTW distance
                                pattern_row_ts <- as.matrix(
                                    dplyr::select(
                                        pattern_row_ts,
                                        -.data[["Index"]]
                                    )
                                )
                                pattern_distance_windows <-
                                    purrr::map(comparison_windows, function(comparison_window) {
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
                                    distance_dtw(value_row_in_window, pattern_row_ts)
                                })
                                pattern_distance_windows <- unlist(
                                    pattern_distance_windows
                                )

                                stats::setNames(
                                    data.frame(distances = unlist(pattern_distance_windows)),
                                    paste0(pattern_name, row_index)
                                )
                            })
                            # Calculate the median distance in each window,
                            # for all pattern behaviors available
                            pattern_distance <- purrr::pmap_dbl(
                                pattern_distances,
                                ~stats::median(c(...))
                            )
                            # Save pattern distance as data.frame
                            pattern_distance <- data.frame(
                                distances = pattern_distance
                            )
                            # Associate the pattern name with the distance and
                            # avoid extra manipulations in the next steps
                            stats::setNames(pattern_distance, pattern_name)
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

#' @export
#'
plot.dtw_model <- function(x) {
    som_obj <- environment(x)[["samples_som_grid"]]
    plot(som_obj)
}
