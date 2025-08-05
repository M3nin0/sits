#' @title Starts transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name for the transition
#' @return expression that creates a STARTS transition call
#' @export
Starts <- function(target_class) {
    substitute(STARTS(values, target_class),
               list(target_class = target_class))
}
#' @title Ends transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name for the transition
#' @return expression that creates a ENDS transition call
#' @export
Ends <- function(target_class) {
    substitute(ENDS(values, target_class),
               list(target_class = target_class))
}
#' @title Edges transition function (start and end are equal in the series)
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name for the transition
#' @return expression that creates a EDGES transition call
#' @export
Edges <- function(target_class) {
    substitute(EDGES(values, target_class),
               list(target_class = target_class))
}
#' @title Dominates transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name for the transition
#' @param size             Neighborhood size
#' @return expression that creates a DOMINATES transition call
#' @export
Dominates <- function(target_class, size) {
    substitute(DOMINATES(values, target_class, size),
               list(target_class = target_class, size = size))
}
#' @title Persist transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name for the transition
#' @param size             Temporal interval size
#' @return expression that creates a PERSIST transition call
#' @export
Persist <- function(target_class, size) {
    substitute(PERSIST(values, target_class, size),
               list(target_class = target_class, size = size))
}
#' @title Peaks transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name for the transition
#' @return expression that creates a PEAKS transition call
#' @export
Peaks <- function(target_class) {
    substitute(PEAKS(values, target_class),
               list(target_class = target_class))
}
#' @title Recur transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name for the transition
#' @return expression that creates a RECUR transition call
#' @export
Recur <- function(target_class) {
    substitute(RECUR(values, target_class),
               list(target_class = target_class))
}
#' @title Convert transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param source_class     Source class name for the transition
#' @param target_class     Target class name for the transition
#' @return expression that creates a CONVERT transition call
#' @export
Convert <- function(source_class, target_class) {
    substitute(
        CONVERT(values, source_class, target_class),
        list(source_class = source_class, target_class = target_class)
    )
}
#' @title Evolve transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param source_class     Source class name for the transition
#' @param target_class     Target class name for the transition
#' @return expression that creates an EVOLVE transition call
#' @export
Evolve <- function(source_class, target_class) {
    substitute(
        EVOLVE(values, source_class, target_class),
        list(source_class = source_class, target_class = target_class)
    )
}
#' @title Keeps transition function
#' @keywords internal
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param target_class     Target class name to keep
#' @return expression that creates a KEEPS transition call
#' @export
Keeps <- function(target_class) {
    substitute(
        KEEPS(values, target_class),
        list(target_class = target_class)
    )
}
#' @title Reclassify a classified cube
#' @name sits_transitions
#'
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @description
#' Apply a set of named expressions to reclassify a classified image.
#' The expressions should use character values to refer to labels in
#' logical expressions.
#'
#' @param cube        Image cube to be reclassified (class = "class_cube")
#' @param  ...        Other parameters for specific functions.
#' @param rules       Expressions to be evaluated (named list).
#' @param memsize     Memory available for classification in GB
#'                    (integer, min = 1, max = 16384).
#' @param multicores  Number of cores to be used for classification
#'                    (integer, min = 1, max = 2048).
#' @param output_dir  Directory where files will be saved
#'                    (character vector of length 1 with valid location).
#' @param version    Version of resulting image (character).
#' @param progress    Set progress bar??
#'
#' @note
#'
#' Reclassification of a remote sensing map refers
#' to changing the classes assigned to different pixels in the image.
#' Reclassification involves assigning new classes to pixels based
#' on additional information from a reference map.
#' Users define rules according to the desired outcome.
#' These rules are then applied to the classified map to produce
#' a new map with updated classes.
#'
#' \code{sits_reclassify()} allow any valid R expression to compute
#' reclassification. User should refer to \code{cube} and \code{mask}
#' to construct logical expressions.
#' Users can use can use any R expression that evaluates to logical.
#' \code{TRUE} values will be relabeled to expression name.
#' Updates are done in asynchronous manner, that is, all expressions
#' are evaluated using original classified values. Expressions are
#' evaluated sequentially and resulting values are assigned to
#' output cube. Last expressions has precedence over first ones.
#'
#' @return An object of class "class_cube" (reclassified cube).
#'
#' @export
sits_transitions <- function(cube, ...) {
    .check_set_caller("sits_transitions")
    UseMethod("sits_transitions", cube)
}

#' @rdname sits_transitions
#' @export
sits_transitions.class_cube <- function(cube, ...,
                                       rules = NULL,
                                       memsize = 4L,
                                       multicores = 2L,
                                       output_dir,
                                       version = "v1",
                                       progress = TRUE) {
    # Preconditions
    .check_is_class_cube(cube)
    # check other params
    .check_int_parameter(memsize, min = 1L, max = 16384L)
    .check_int_parameter(multicores, min = 1L, max = 2048L)
    .check_output_dir(output_dir)
    # Check version and progress
    version <- .message_version(version)
    progress <- .message_progress(progress)
    # Define optimal parameters for parallel processing
    # Get block size
    block <- .raster_file_blocksize(.raster_open_rast(.tile_path(cube)))
    # Check minimum memory needed to process one block
    job_block_memsize <- .jobs_block_memsize(
        block_size = .block_size(block = block, overlap = 0L),
        npaths = 2L,
        nbytes = 8L, proc_bloat = .conf("processing_bloat_cpu")
    )
    # Update multicores parameter
    multicores <- .jobs_max_multicores(
        job_block_memsize = job_block_memsize,
        memsize = memsize,
        multicores = multicores
    )
    # Update block parameter based on the size of memory and number of cores
    block <- .jobs_optimal_block(
        job_block_memsize = job_block_memsize,
        block = block,
        image_size = .tile_size(.tile(cube)),
        memsize = memsize,
        multicores = multicores
    )
    # Prepare parallelization
    .parallel_start(workers = multicores)
    on.exit(.parallel_stop(), add = TRUE)
    # Capture expression
    rules <- as.list(substitute(rules, environment()))[-1L]
    # Prepare transitions function
    transitions_fn <- .transitions_fn_expr(
        rules = rules,
        labels_cube = .cube_labels(cube)
    )
    # Process each tile sequentially
    class_cube <- .cube_foreach_tile(cube, function(tile) {
        # Get new labels from cube and pre-defined rules from user
        cube_labels <- .transitions_new_labels(cube, rules)
        # Transitions reasoning for reclassification
        .transitions_tile(
            tile = tile,
            band = "class",
            labels = cube_labels,
            transitions_fn = transitions_fn,
            block = block,
            output_dir = output_dir,
            version = version,
            progress = progress
        )
    })
    class(class_cube) <- c("class_cube", class(class_cube))
    return(class_cube)
}
#' @rdname sits_transitions
#' @export
sits_transitions.default <- function(cube, ...) {
    stop(.conf("messages", "sits_transitions"))
}
