
#' @title Extract datetime from a STAC Query.
#' @keywords internal
#' @noRd
#'
#' @param stac_query Query that follows the STAC protocol.
#' @return List with `start_date` and `end_date` properties.
#' @export
.stac_cdse_datetime_as_dates <- function(stac_query) {
    query_datetime <- stringr::str_split(
        stac_query[["params"]][["datetime"]], "/"
    )
    list(
        start_date = query_datetime[[1]][1],
        end_date = query_datetime[[1]][2]
    )
}

#' @title Extract bounding box from a STAC Query.
#' @keywords internal
#' @noRd
#'
#' @param stac_query Query that follows the STAC protocol.
#' @return List with `bbox` property.
#' @export
.stac_cdse_intersects_as_bbox <- function(stac_query) {
    intersects <- stac_query[["params"]][["intersects"]]
    coordinates <- intersects[["coordinates"]]
    # Extract x-coordinates and y-coordinates
    coordinates_x <- coordinates[,,1]
    coordinates_y <- coordinates[,,2]
    # Calculate bounding box
    min_x <- min(coordinates_x)
    max_x <- max(coordinates_x)
    min_y <- min(coordinates_y)
    max_y <- max(coordinates_y)
    # Create bbox object
    list(bbox = c(min_y, min_x, max_y, max_x))
}

#' @title Fix STAC Items from CDSE with assets metadata.
#' @keywords internal
#' @noRd
#'
#' @description
#'  This auxiliary function creates `STAC Assets` for each of the bands selected
#'  by user. This is necessary since CDSE STAC does not return band files as
#'  `STAC Assets` The information is extracted from files available in the
#'  CDSE S3 API.
#' @param source     Data source.
#' @param items      STAC items
#' @param bands      Band names
#' @param collection Image collection.
#' @return           STAC Items updated with `assets` property.
#' @export
.stac_cdse_fix_items <- function(source, items, bands, collection) {
    # Define path used to extract the product prefix in CDSE S3
    s3_path <- c("assets", "PRODUCT", "alternate", "s3", "href")
    # Define name of the CDSE products bucket
    s3_bucket <- "/eodata/"
    # Define protocol to access S3 data
    s3_protocol <- "/vsis3"
    # CDSE does not provide files directly in the `assets` property`. As a
    # workaround, we build the `assets` property using files in the CDSE S3 API
    features_fixed <- purrr::map(items[["features"]], function(feature) {
        # Prepare item path in CDSE S3
        item_s3_path <- rstac::items_reap(feature, s3_path) |>
            stringr::str_replace(s3_bucket, "")
        # Read the item's content from CDSE S3
        item_s3_content <- aws.s3::get_bucket_df(
            bucket = "eodata",
            use_https = TRUE,
            region = "",
            prefix = item_s3_path
        )
        # Extract the address of the files associated with bands selected by
        # users.
        item_bands <- purrr::map_df(bands, function(band) {
            band_conf <- .conf_eo_band(source, collection, band)
            # Create pattern to select file associated with the band.
            band_pattern <- band_conf[["pattern"]]
            # Filter the S3 content to get files from the band
            band_item <-
                dplyr::filter(item_s3_content,
                              stringr::str_detect(.data[["Key"]], band_pattern))
            # Check if the correct file was selected.
            .check_that(nrow(band_item) == 1,
                        msg = paste0("Invalid band: ", band))
            # Prepare the file address
            band_path_s3 <- paste0(s3_protocol, s3_bucket, band_item[["Key"]])
            # Prepare result and return it
            # As this auxiliary function only needs to provide the right content
            # to other parts of `sits`, only the `href` of the image is returned.
            # The other necessary actions are managed by `sits.`
            stats::setNames(list(band = list(href = band_path_s3)), band)
        })
        feature[["assets"]] <- item_bands
        feature
    })
    items[["features"]] <- features_fixed
    items
}

#' @title Query scenes available in the CDSE Open Search.
#' @keywords internal
#' @noRd
#'
#' @description
#'  This auxiliary function is used to get scenes available in CDSE using
#'  Open Search. This is required as the current version of the CDSE STAC does
#'  not support fields / advanced search (e.g., search by product type, cloud
#'  coverage).
#' @param source     Data source
#' @param collection Open Search collection endpoint.
#' @param bands      Band names
#' @param bbox       Bounding box of the area from data must be from
#' @param start_date Start date.
#' @param end_date   End date.
#' @param product_type Type of the CDSE Product (e.g., S2MSI2A)
#' @param filter_fn  Function to add an extra filter to the scenes found.
#' @return           List of Scene IDS
#' @export
.opensearch_cdse_query_scenes <- function(source,
                                          collection,
                                          start_date,
                                          end_date,
                                          bbox,
                                          product_type,
                                          limit = 1000,
                                          filter_fn = NULL) {
    # CDSE Open Search configurations
    cdse_opensearch_base_url <- .conf(
        "sources",
        source,
        "service_auxiliary",
        "opensearch"
    )
    cdse_opensearch_max_items <- limit
    cdse_opensearch_endpoint <- "search.json"
    # Create the Open Search endpoint for the collection
    # Selected by user
    collection_url <- paste(
        cdse_opensearch_base_url,
        collection,
        cdse_opensearch_endpoint,
        sep = "/"
    )
    # Define features to save content from Open Search
    features <- c()
    # Define variables to support the pagination in the Open Search
    current_page <- 1
    is_to_fetch_more <- TRUE
    # Prepare bounding box in the format required by Open Search
    bbox <- paste(bbox, collapse = ",")
    # Get items from Open Search (with pagination)
    while(is_to_fetch_more) {
        # Get raw content from Open Search API
        response <- httr::GET(url = collection_url, query = list(
            startDate      = start_date,
            completionDate = end_date,
            maxRecords     = cdse_opensearch_max_items,
            page           = current_page,
            box            = bbox,
            productType    = product_type,
            status         = "ONLINE"
        ))
        .check_int_parameter(httr::status_code(response),
                             min = 200,
                             max = 200,
                             msg = "Failed to fetch data")
        # Extract data from the response
        page_data <- httr::content(response, "parsed")
        # Save results
        features <- c(features, page_data$features)
        # Check if is required to fetch more
        # if less than `1000` means the pagination ends
        is_to_fetch_more <- length(page_data$features) >=
            cdse_opensearch_max_items
        # Prepare next page fetch
        current_page <- current_page + as.numeric(is_to_fetch_more)
    }
    # Apply user-defined filter if available.
    if (!is.null(filter_fn)) {
        features_valid <- unlist(purrr::map(features, filter_fn))
        features <- features[features_valid]
    }
    # Extract scene IDs
    scene_ids <- unlist(purrr::map(features, function(feature) {
        feature[["properties"]][["title"]]
    }))
    unique(scene_ids)
}

#' @title Test access to STAC collection
#' @keywords internal
#' @noRd
#'
#' @description
#' These functions provide an API to handle/retrieve data from source's
#' collections.
#'
#' @param source     Data source.
#' @param collection Image collection.
#' @param bands      Band names
#' @param ...        Other parameters to be passed for specific types.
#' @param start_date Start date.
#' @param end_date   End date.
#' @param dry_run    TRUE/FALSE
#' @return           Called for side effects
#' @export
.source_collection_access_test.cdse_cube <- function(source, collection,
                                                     bands, ...,
                                                     start_date = NULL,
                                                     end_date = NULL,
                                                     dry_run = FALSE) {
    # Currently CDSE STAC filters are limited. As there is no possibility of
    # using specific selections, sometimes using the first item returned from
    # STAC can be a problem (e.g., Auxiliary products). To avoid this problem,
    # we use a pre-defined image.
    item_for_test <- .conf(
        "sources",
        source,
        "collections",
        collection,
        "item_test"
    )
    # require package
    .check_require_packages("rstac")
    # create a query
    items_query <- .stac_create_items_query(
        source = source,
        collection = collection,
        limit = 1
    )
    items_query <- rstac::stac_search(items_query, ids = c(item_for_test))
    # assert that service is online
    tryCatch(
        {
            items <- rstac::post_request(items_query, ...)
        },
        error = function(e) {
            stop(
                paste(
                    ".source_collection_access_test.stac_cube: service is",
                    "unreachable\n", e$message
                ),
                call. = FALSE
            )
        }
    )
    .check_stac_items(items)
    items <- .source_items_bands_select(
        source = source,
        items  = items,
        bands  = bands[[1]],
        collection = collection, ...
    )
    href <- .source_item_get_hrefs(
        source = source,
        item = items$feature[[1]],
        collection = collection, ...
    )
    # assert that token and/or href is valid
    if (dry_run) {
        tryCatch(
            {
                .raster_open_rast(href)
            },
            error = function(e) {
                stop(paste(
                    ".source_collection_access_test.stac_cube: cannot",
                    "open url\n", href, "\n", e$message
                ), call. = FALSE)
            }
        )
    }
    return(invisible(source))
}

#' @title Select bands from a STAC item
#' @keywords internal
#' @noRd
#'
#' @param source     Data source
#' @param items      STAC items
#' @param bands      Bands to be selected in the collection.
#' @param collection Image collection
#' @param ...        Additional parameters.
#' @return List of STAC items
#' @export
.source_items_bands_select.cdse_cube <- function(source,
                                                 items,
                                                 bands,
                                                 collection, ...) {
    # CDSE does not provide files in the `assets` property. So, it is
    # required to fix this using content from CDSE S3 API.
    items <- .stac_cdse_fix_items(source, items, bands, collection)
    # With the correct metadata, it is possible to use the base workflow.
    items <- .source_items_bands_select.stac_cube(
        source,
        items,
        bands,
        collection,
        ...
    )
    return(items)
}

#' @title Transform an items object in a CDSE cube
#' @keywords internal
#' @noRd
#' @description \code{.source_items_new()} this function is called to create
#' an items object. In case of Web services, this function is responsible for
#' making the Web requests to the server.
#' @param source     Name of the STAC provider.
#' @param ...        Other parameters to be passed for specific types.
#' @param collection Collection to be searched in the data source.
#' @param stac_query Query that follows the STAC protocol
#' @param tiles      Selected tiles (not supported)
#' @param platform   Satellite platform (optional).
#' @return An object referring the images of a sits cube.
#'
#' @export
.source_items_new.cdse_cube <- function(source, ...,
                                        collection,
                                        stac_query,
                                        tiles = NULL, # TODO: Not supported
                                        platform = NULL) {
    # set caller to show in errors
    .check_set_caller(".source_items_new.cdse_cube")
    # define the maximum number of records per request
    cdse_query_limit <- 1000
    # as CDSE STAC returns many types of items in the same collection,
    # it is required to filter the content by a specific type.
    item_type <- .conf(
        "sources",
        source,
        "collections",
        collection,
        "item_type"
    )
    # extract collection endpoint
    collection_endpoint <- .conf(
        "sources",
        source,
        "collections",
        collection,
        "collection_name"
    )
    # extract query parameters from STAC Query
    query_bbox <- .stac_cdse_intersects_as_bbox(stac_query)
    query_date <- .stac_cdse_datetime_as_dates(stac_query)
    # Before using STAC, we go to Open Search to get all Scene IDs matching
    # the user-defined query. This is required as CDSE STAC does not provide
    # fields / advanced search.
    valid_scene_ids <- .opensearch_cdse_query_scenes(
        source       = source,
        collection   = collection_endpoint,
        start_date   = query_date$start_date,
        end_date     = query_date$end_date,
        bbox         = query_bbox$bbox,
        product_type = item_type,
        limit        = cdse_query_limit,
        filter_fn    = function(feature) {
            # Function used to remove `DESCENDING` data from results.
            # This option is available in the
            is_valid <- TRUE
            feature_properties <- feature[["properties"]]
            if ("orbitDirection" %in% feature_properties) {
                is_valid <-
                    feature_properties[["orbitDirection"]] != "DESCENDING"
            }
            is_valid
        }
    )
    # CDSE presented some errors when many IDs were checked at the same time.
    # To avoid errors, we are requesting the scenes in chunks
    chunk_size <- 250  # value identified as reasonable
    # Check using chunks only if many scenes were requested
    if (length(valid_scene_ids) > chunk_size) {
        chunks <- (seq_along(valid_scene_ids) - 1) %/% chunk_size
        chunks <- split(valid_scene_ids, chunks)

        items_info <- purrr::map(chunks, function(chunk) {
            # create a new STAC search for each chunk
            stac_query <- rstac::stac(stac_query[["base_url"]])
            # prepare a new STAC query using the IDS extracted from Open Search.
            stac_query <- rstac::stac_search(stac_query,
                                             ids = chunk,
                                             limit = cdse_query_limit)
            rstac::post_request(q = stac_query, ...)
        })
        # use first result as items
        items <- items_info[[1]]
        # extract features
        items[["features"]] <- purrr::reduce(
            purrr::map(items_info, purrr::pluck, "features"), c
        )
        # sanity check
        .check_that(length(items[["features"]]) == length(valid_scene_ids),
                    msg = "It was not possible to process the results of your
                          query.")

        items_info <- items
    } else {
        # prepare a new STAC query using the IDS extracted from Open Search.
        stac_query <- rstac::stac(stac_query[["base_url"]])
        stac_query <- rstac::stac_search(stac_query,
                                         ids = valid_scene_ids,
                                         limit = cdse_query_limit)
        # get metadata
        items_info <- rstac::post_request(q = stac_query, ...)
    }
    .check_stac_items(items_info)
    return(items_info)
}

#' @title Organizes items by tiles for CDSE collections
#' @param source     Name of the STAC provider.
#' @param ...        Other parameters to be passed for specific types.
#' @param items      \code{STACItemcollection} object from rstac package.
#' @param collection Collection to be searched in the data source.
#' @return A list of STAC items.
#' @keywords internal
#' @noRd
#' @export
.source_items_tile.cdse_cube <- function(source, ...,
                                         items,
                                         collection = NULL) {
    rstac::items_reap(items, field = c("properties", "tileId"))
}
