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
        # return a valid detect change method
        return(result)
    }
