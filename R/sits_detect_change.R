#' @title Detect changes in time series
#' @name sits_detect_change
#' @description Given a set of time series or an image, this function compares each time series
#'                         with a set of change/no-change patterns, and indicates places and dates where
#'                          change has been detected
#' @param data         Set of time series
#' @param cd_method    Change detection method (with parameters)
#' @return             Set of time series where significant change has been detected
#' @export
sits_detect_change <- function(data,
                               cd_method,
                               ...,
                               filter_fn = NULL,
                               multicores = 2L,
                               progress = TRUE) {
    UseMethod("sits_detect_change", data)
}

#' @rdname sits_detect_change
#' @export
sits_detect_change.sits <- function(data,
                                    cd_method,
                                    ...,
                                    filter_fn = NULL,
                                    multicores = 2L,
                                    progress = TRUE) {
    # Pre-conditions
    data <- .check_samples_ts(data)
    # TODO:
    #  |> Create function to check the `cd_method` or generalize
    #  |> `.check_is_sits_model`.
    # .check_is_sits_model(ml_model)
    .check_multicores(multicores, min = 1, max = 2048)
    .check_progress(progress)
    # To change detection
    detections <- .detect_change_ts(
        samples = data,
        cd_method = cd_method,
        filter_fn = filter_fn,
        multicores = multicores,
        progress = progress
    )
    return(detections)
}

#' @rdname sits_detect_change
#' @export
sits_detect_change.default <- function(data, cd_method, ...) {
    stop("Input should be a sits tibble or a data cube")
}
