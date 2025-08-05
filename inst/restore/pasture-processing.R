reclassify_timeseries_chunk <- function(files, reference_class_number, neighbor_class_number, multicores, output_dir) {
    output_dir <- fs::path(output_dir)
    fs::dir_create(output_dir)

    stopifnot(is.character(output_dir))
    block <- c("nrows" = 1024, "ncols" = 1024)

    rast_template <- sits:::.raster_open_rast(files)

    image_size <- list(nrows = sits:::.raster_nrows(rast_template),
                       ncols = sits:::.raster_ncols(rast_template))

    chunks <- sits:::.chunks_create(
        block = block,
        overlap = 0,
        image_size = image_size,
        image_bbox = sits:::.bbox(sits:::.raster_bbox(rast_template))
    )

    sits:::.parallel_start(workers = multicores)
    on.exit(sits:::.parallel_stop(), add = TRUE)

    block_files <- sits:::.jobs_map_parallel_chr(chunks, function(chunk) {
        block <- sits:::.block(chunk)
        values <- sits:::.raster_read_rast(
            files = files, block = block
        )

        values <- sits:::transition_neighbor_analysis(
            data = values,
            reference_class = reference_class_number,
            neighbor_class = neighbor_class_number
        )

        block_file <- sits:::.file_block_name(
            pattern = tools::file_path_sans_ext("transitions-analysis.tif"),
            block = block,
            output_dir = output_dir
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

    out_file <- fs::path(output_dir) / "transitions-analysis.tif"

    sits:::.raster_merge_blocks(
        out_files = out_file,
        base_file = file,
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
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2015/map-eco3-2015.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2016/map-eco3-2016.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2017/map-eco3-2017.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2018/map-eco3-2018.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2019/map-eco3-2019.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2020/map-eco3-2020.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2021/map-eco3-2021.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2022/map-eco3-2022.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2023/map-eco3-2023.tif",
    "~/Dropbox/projects/01_restore-plus/11_classification-eco3-mcti-v2/2024/map-eco3-2024.tif"
)

reclassified_raster <- reclassify_timeseries_chunk(
    files = files,
    reference_class_number = 2,
    neighbor_class_number = 12,
    multicores = 10,
    output_dir = "~/Downloads/reclassify-pasture"
)

