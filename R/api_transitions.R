#' @title Replace expression class with numeric index
#' @keywords internal
#' @noRd
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param expr             Expression to be processed
#' @param labels           Named vector of labels where names are indices
#' @return numeric index if expr matches a label, otherwise the original expr
.transitions_replace_expr_class <- function(expr, labels) {
    if (is.character(expr) && length(expr) == 1) {
        idx <- names(labels)[labels == expr]
        if (length(idx) == 1) return(as.integer(idx))
    }
    return(expr)
}
#' @title Expand and transform expressions for transition rules
#' @keywords internal
#' @noRd
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param expr    Expression to be expanded (can be a call, symbol, or literal)
#' @param labels  Named vector of labels where names are indices
#' @return expanded expression with string labels replaced by numeric indices
#'
#' @description
#' This function recursively processes expressions used in transition rules.
#' It handles different types of expressions:
#' - Function calls (including special macros like Convert, Recur, Keeps)
#' - Binary operators (|, &)
#' - Unary operators (!)
#' - Parentheses
#' - String literals (replaced with numeric indices)
.transitions_expand_expr <- function(expr, labels) {
    # Define name of predicate names available in C++
    available_predictors <- c(
        "Convert", "Recur", "Keeps", "Evolve",
        "Starts", "Ends", "Persist", "Peaks",
        "Edges"
    )
    # Check if expression is a function call
    if (is.call(expr)) {
        # Extract function name from the call
        fname <- as.character(expr[[1]])
        # Handle special macro functions: Convert, Recur, Keeps
        # These are transformed to uppercase and get 'values' as first argument
        if (fname %in% available_predictors) {
            # Convert to uppercase (e.g., Convert -> CONVERT)
            new_fname <- toupper(fname)
            # Extract function arguments (remove function name)
            args <- as.list(expr)[-1]
            # Recursively expand each argument and replace string labels
            args <- lapply(args, function(arg)
                .transitions_replace_expr_class(
                    .transitions_expand_expr(arg, labels), labels)
                )
            # Reconstruct call with uppercase function name and 'values'
            return(as.call(c(
                as.name(new_fname), quote(values), args
            )))
            # Handle binary logical operators: | (OR), & (AND)
        } else if (fname %in% c("|", "&")) {
            # Recursively expand left and right operands
            return(as.call(
                list(
                    as.name(fname),
                    .transitions_expand_expr(expr[[2]], labels),
                    # Left operand
                    .transitions_expand_expr(expr[[3]], labels)
                )
            )) # Right operand

        # Handle unary logical operator: ! (NOT)
        } else if (fname == "!") {
            # Recursively expand the operand
            return(as.call(list(
                as.name("!"),
                .transitions_expand_expr(expr[[2]], labels)
            )))

        # Handle parentheses for grouping expressions
        } else if (fname == "(") {
            # Recursively expand the expression inside parentheses
            return(as.call(list(
                as.name("("),
                .transitions_expand_expr(expr[[2]], labels)
            )))
        # Handle all other function calls generically
        } else {
            # Recursively expand all arguments of the function call
            return(as.call(
                lapply(expr, function(e)
                    .transitions_expand_expr(e, labels))
            ))
        }
    } else {
        # For non-call expressions (symbols, literals, etc.)
        # Replace string literals with numeric indices if they match labels
        return(.transitions_replace_expr_class(expr, labels))
    }
}

#' @title Reasoning transitions in a tile
#' @keywords internal
#' @noRd
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param  tile.           Subset of a data cube
#' @param  mask            Reclassification mask
#' @param  band            Output band
#' @param  labels          Output labels
#' @param  transitions_fn  Function to be applied for transition reasoning
#' @param  block           Image block to be processed
#' @param  output_dir      Directory where image will be save
#' @param  version         Version of result.
#' @param  progress        Show progress bar?
#' @return transitions reasoning results tile
.transitions_tile <- function(tile, band, labels, transitions_fn, block,
                             output_dir, version, progress) {
    # Output files
    out_file <- .file_derived_name(
        tile = tile, band = band, version = version, output_dir = output_dir
    )
    # Resume feature
    if (file.exists(out_file)) {
        .check_recovery()
        class_tile <- .tile_derived_from_file(
            file = out_file,
            band = band,
            base_tile = tile,
            derived_class = "class_cube",
            labels = labels,
            update_bbox = FALSE
        )
        # Update tile labels
        class_tile <- .tile_update_label(class_tile, labels)
        return(class_tile)
    }
    # Create chunks as jobs
    chunks <- .tile_chunks_create(tile = tile, block = block, overlap = 0L)
    # start parallel process
    block_files <- .jobs_map_parallel_chr(chunks, function(chunk) {
        # Get job block
        block <- .block(chunk)
        # Output file name
        block_file <- .file_block_name(
            pattern = .file_pattern(out_file),
            block = block,
            output_dir = output_dir
        )
        # Resume processing in case of failure
        if (.raster_is_valid(block_file)) {
            return(block_file)
        }
        # Project mask block to template block
        # Get band conf missing value
        band_conf <- .conf_derived_band(
            derived_class = "class_cube", band = band
        )
        # Read and preprocess values
        values <- .tile_read_block(
            tile = tile, band = .tile_bands(tile), block = block
        )
        # Evaluate expressions
        values <- transitions_fn(values = values)
        # Does values is valid? In case of a matrix with integer(0) values
        if (.has_not(values)) {
            values <- rep(NA, .block_size(block))
        }
        offset <- .offset(band_conf)
        if (.has(offset) && offset != 0.0) {
            values <- values - offset
        }
        scale <- .scale(band_conf)
        if (.has(scale) && scale != 1.0) {
            values <- values / scale
        }
        # Prepare and save results as raster
        .raster_write_block(
            files = block_file, block = block, bbox = .bbox(chunk),
            values = values, data_type = .data_type(band_conf),
            missing_value = .miss_value(band_conf),
            crop_block = NULL
        )
        # Free memory
        gc()
        # Returned value
        block_file
    }, progress = progress)
    # Merge blocks into a new class_cube tile
    class_tile <- .tile_derived_merge_blocks(
        file = out_file,
        band = band,
        labels = labels,
        base_tile = tile,
        block_files = block_files,
        derived_class = "class_cube",
        multicores = .jobs_multicores(),
        update_bbox = FALSE
    )
    # Update tile labels
    class_tile <- .tile_update_label(class_tile, labels)
    # Return class tile
    class_tile
}

#' @title Reclassify function
#' @keywords internal
#' @noRd
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param  rules           Rules to be applied
#' @param  labels_cube     Labels of input cube
#' @return function to be applied for transition reasoning
.transitions_fn_expr <- function(rules, labels_cube) {
    .check_set_caller(".transitions_fn_expr")
    # Check if rules are named
    .check_that(all(.has_name(rules)))
    # Get output labels
    labels_rule <- setdiff(names(rules), labels_cube)
    names(labels_rule) <- max(.as_int(names(labels_cube))) +
        seq_along(labels_rule)
    labels <- c(labels_cube, labels_rule)
    labels_code <- .as_int(names(labels))
    # Define reclassify function
    transitions_fn <- function(values) {
        # Used to check values (below)
        input_pixels <- nrow(values)
        # New evaluation environment
        env <- list2env(list(
            values = values
        ))
        # Get values as character
        values <- rep(0, input_pixels)
        # Evaluate each expression
        for (label in names(rules)) {
            # Get expression
            expr <- rules[[label]]
            # Expand expression with class names and C++ function names
            expr <- .transitions_expand_expr(expr, labels_cube)
            # Evaluate
            result <- eval(expr, envir = env)
            # Update values
            if (!is.logical(result)) {
                stop(.conf("messages", ".transitions_fn_results"))
            }
            values[result] <- label
        }
        # Get values as numeric
        values <- matrix(
            data = labels_code[match(values, labels)],
            nrow = input_pixels
        )
        # Are the results consistent with the data input?
        .check_processed_values(values, input_pixels)
        # Return values
        values
    }
    # Return closure
    transitions_fn
}
#' @title Obtain new labels on reclassification operation
#' @keywords internal
#' @noRd
#' @author Felipe Carlos, \email{efelipecarlos@@gmail.com}
#' @author Felipe Carvalho, \email{felipe.carvalho@@inpe.br}
#' @author Gilberto Camara, \email{gilberto.camara@@inpe.br}
#'
#' @param  cube            Labelled data cube
#' @param  rules           Rules to be applied
#' @return new labels to be applied to the cube
.transitions_new_labels <- function(cube, rules) {
    # Get cube labels
    cube_labels <- .cube_labels(cube)
    # Get rules new labels
    new_labels <- setdiff(names(rules), cube_labels)
    # Does rules has new labels in the composition?
    if (.has(new_labels)) {
        # Get the next index
        next_idx <- max(as.numeric(names(cube_labels))) + 1L
        idx_values <- seq.int(
            from = next_idx, to = next_idx + length(new_labels) - 1L
        )
        names(new_labels) <- as.character(idx_values)
    }
    c(cube_labels, new_labels)
}
