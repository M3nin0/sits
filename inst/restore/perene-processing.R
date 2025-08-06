reclassify_timeseries_chunk <- function(files,
                                        reference_class_number,
                                        neighbor_class_number,
                                        version,
                                        multicores,
                                        memsize,
                                        output_dir) {
    output_dir <- fs::path(output_dir)
    fs::dir_create(output_dir)
    stopifnot(is.character(output_dir))
    # block <- c("nrows" = 5000, "ncols" = 3000)

    out_filename <- paste0("transitions-analysis", "-", version, ".tif")
    out_file <- fs::path(output_dir) / out_filename

    if (file.exists(out_file)) {
        return(out_file)
    }

    rast_template <- sits:::.raster_open_rast(files)
    image_size <- list(
        nrows = sits:::.raster_nrows(rast_template),
        ncols = sits:::.raster_ncols(rast_template)
    )

    block <- sits:::.raster_file_blocksize(sits:::.raster_open_rast(files))
    # Check minimum memory needed to process one block
    job_block_memsize <- sits:::.jobs_block_memsize(
        block_size = sits:::.block_size(block = block, overlap = 0),
        npaths = (length(files) * terra::nlyr(rast_template)),
        nbytes = 8,
        proc_bloat = sits:::.conf("processing_bloat")
    )
    # Update multicores parameter based on size of a single block
    multicores <- sits:::.jobs_max_multicores(
        job_block_memsize = job_block_memsize,
        memsize = memsize,
        multicores = multicores
    )
    # Update block parameter based on the size of memory and number of cores
    block <- sits:::.jobs_optimal_block(
        job_block_memsize = job_block_memsize,
        block = block,
        image_size = image_size,
        memsize = memsize,
        multicores = multicores
    )
    # Create chunks
    chunks <- sits:::.chunks_create(
        block = block,
        overlap = 0,
        image_size = image_size,
        image_bbox = sits:::.bbox(
            sits:::.raster_bbox(rast_template),
            default_crs = terra::crs(rast_template)
        )
    )
    # Start workers
    sits:::.parallel_start(workers = multicores)
    on.exit(sits:::.parallel_stop(), add = TRUE)
    # Process data!
    block_files <- sits:::.jobs_map_parallel_chr(chunks, function(chunk) {
        block <- sits:::.block(chunk)
        block_file <- sits:::.file_block_name(
            pattern = tools::file_path_sans_ext(out_filename),
            block = block,
            output_dir = output_dir
        )

        if (file.exists(block_file)) {
            return(block_file)
        }

        values <- sits:::.raster_read_rast(files = files, block = block)

        values <- sits:::transition_neighbor_analysis(
            data = values,
            reference_class = reference_class_number,
            neighbor_class = neighbor_class_number
        )

        sits:::.raster_write_block(
            files = block_file,
            block = block,
            bbox = sits:::.bbox(chunk),
            values = values,
            data_type = "INT1U",
            missing_value = 255,
            crop_block = NULL
        )
        block_file
    }, progress = TRUE)

    sits:::.raster_merge_blocks(
        out_files = out_file,
        base_file = files,
        block_files = block_files,
        data_type = "INT1U",
        missing_value = 255,
        multicores = multicores
    )

    # remove block files
    unlink(block_files)

    return(out_file)
}

files <- c(
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2015/LANDSAT_OLI_MOSAIC_2015-01-01_2015-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2016/LANDSAT_OLI_MOSAIC_2016-01-01_2016-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2017/LANDSAT_OLI_MOSAIC_2017-01-01_2017-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2018/LANDSAT_OLI_MOSAIC_2018-01-01_2018-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2019/LANDSAT_OLI_MOSAIC_2019-01-01_2019-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2020/LANDSAT_OLI_MOSAIC_2020-01-01_2020-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2021/LANDSAT_OLI_MOSAIC_2021-01-01_2021-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2022/LANDSAT_OLI_MOSAIC_2022-01-01_2022-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2023/LANDSAT_OLI_MOSAIC_2023-01-01_2023-12-01_class_mask-clean-step8.tif",
    "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/2024/LANDSAT_OLI_MOSAIC_2024-01-01_2024-12-01_class_mask-clean-step8.tif"
)

# purrr::map(files, function(file) {
#     rst <- terra::rast(file)
#
#     paste0("nrow = ", terra::nrow(rst), " | ", "ncol = ", terra::ncol(rst))
# })

output_dir <- "/data/experiments/water-mask-variations/data/derived/masks/mask-mcti-v3/transitions/test-reclassify-temporal"
fs::dir_create(output_dir)

ismain <- TRUE

if (ismain) {
    # if "Ag_perene" -> "vegetacao_secundaria" -> "Ag_perene", then transform "vegetacao_secundaria" to "Ag_perene"
    reclassified_file <- reclassify_timeseries_chunk(
        files = files,
        reference_class_number = 12, # 12 = "vegetacao_secundaria"
        neighbor_class_number = 2, # 2 = "Ag_perene"
        multicores = 16,
        memsize = 100,
        version = "v3",
        output_dir = output_dir
    )

    reclassified_raster <- terra::rast(reclassified_file)

    for (idx in seq_len(length(files))) {
        file_path <- files[[idx]]
        file_out_path <- stringr::str_replace(file_path, ".tif", "-perene-reclass.tif")

        message("Processing: ",
                basename(reclassified_file),
                " → ",
                basename(file_out_path))

        sf::gdal_utils(
            util = "translate",
            source = as.character(fs::path_expand(reclassified_file)),
            destination = file_out_path,
            options = sits:::.gdal_params(
                list(
                    "-b"     = as.character(idx),
                    "-of"    = "GTiff",
                    "-co"    = "TILED=YES",
                    "-co"    = "COMPRESS=LZW",
                    "-co"    = "INTERLEAVE=BAND",
                    "-co"    =  "PREDICTOR=2"
                )
            ),
            quiet = FALSE
        )

        sf::gdal_addo(file_out_path)
    }

    # cmd <- "gdal_translate"
    # args <- c(
    #     "-b", as.character(idx),
    #     "-of", "GTiff",
    #     "-co", "TILED=YES",
    #     "-co", "COMPRESS=LZW",
    #     "-co", "COPY_SRC_OVERVIEWS=YES",
    #     shQuote(as.character(reclassified_file)),
    #     shQuote(file_out_path)
    # )
    # system2(cmd, args)
}
