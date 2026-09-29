library(dplyr)
library(readxl)
library(terra)

rm(list = ls(all = TRUE))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
package_dir <- dirname(dirname(script_file))
data_dir <- file.path(package_dir, "data")
input_dir <- file.path(package_dir, "input_data")

allocation_dir <- file.path(input_dir, "indicator_allocations")
statistic_file <- file.path(input_dir, "statistics", "chn_water_save_irr_statistics_province_level_mean.xlsx")
province_shp <- file.path(input_dir, "boundaries", "province", "province.shp")
irrigation_mask_dir <- file.path(input_dir, "irrigation_mask")
output_dir <- Sys.getenv("FIGSHARE_DATA_DIR", data_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

years <- c(2000, 2005, 2010, 2015, 2020)
indicator_order <- c("det", "iwu", "dkndvi_by_det_bycrop", "det_by_iwu_bycrop_q90")
priority_weight <- 2 ^ rev(seq_along(indicator_order) - 1)
names(priority_weight) <- indicator_order
priority_max <- length(indicator_order) * 100 + sum(priority_weight)

standardize_name <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("维吾尔自治区$", "", x)
    x <- sub("壮族自治区$", "", x)
    x <- sub("回族自治区$", "", x)
    x <- sub("自治区$", "", x)
    x <- sub("特别行政区$", "", x)
    x <- sub("省$", "", x)
    x <- sub("市$", "", x)
    x
}

share_to_ratio <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[is.na(x)] <- 0
    pmax(pmin(x / 100, 1), 0)
}

target_count <- function(n, share) {
    if (!is.finite(n) || n <= 0 || !is.finite(share) || share <= 0) return(0L)
    as.integer(min(n, ceiling(n * share)))
}

indicator_file <- function(indicator, year) {
    file.path(allocation_dir, indicator, sprintf("irr_type_%s_%d.tif", indicator, year))
}

prepare_statistic <- function() {
    statistic_df <- as.data.frame(read_excel(statistic_file))
    required_cols <- c("province_name", "mapped_year", "drip_frac", "sprinkler_frac")
    missing_cols <- setdiff(required_cols, names(statistic_df))
    if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

    statistic_df %>%
        mutate(
            province_name = standardize_name(province_name),
            year = as.integer(mapped_year),
            drip_share = share_to_ratio(drip_frac),
            sprinkler_share = share_to_ratio(sprinkler_frac)
        ) %>%
        select(province_name, year, drip_share, sprinkler_share)
}

prepare_province <- function() {
    province <- terra::vect(province_shp)
    province$province_code <- as.integer(province[["省代码"]][, 1])
    province$province_name <- standardize_name(province[["省"]][, 1])
    lookup <- unique(as.data.frame(province)[, c("province_code", "province_name")])
    lookup <- lookup[!is.na(lookup$province_code), ]
    lookup <- lookup[!lookup$province_code %in% c(710000L, 810000L, 820000L), ]
    lookup <- lookup[order(lookup$province_code), ]
    list(province = province, lookup = lookup)
}

calc_priority <- function(hit_mat) {
    support_count <- rowSums(hit_mat, na.rm = TRUE)
    tie_break_score <- as.numeric(hit_mat %*% priority_weight)
    support_count * 100 + tie_break_score
}

rank_cells_by_priority <- function(cells, priority, n_select, require_positive = TRUE) {
    if (n_select <= 0 || length(cells) == 0) return(integer(0))
    keep <- is.finite(priority[cells])
    if (require_positive) keep <- keep & priority[cells] > 0
    cells <- cells[keep]
    if (length(cells) == 0) return(integer(0))
    order_idx <- order(priority[cells], cells, decreasing = c(TRUE, FALSE), na.last = NA)
    cells[order_idx][seq_len(min(length(order_idx), n_select))]
}

select_cells_with_supplement <- function(preferred_cells, supplement_cells, priority, n_select) {
    selected_cells <- rank_cells_by_priority(preferred_cells, priority, n_select, TRUE)
    if (length(selected_cells) < n_select) {
        extra_cells <- rank_cells_by_priority(
            setdiff(supplement_cells, selected_cells), priority,
            n_select - length(selected_cells), FALSE
        )
        selected_cells <- c(selected_cells, extra_cells)
    }
    selected_cells
}

allocate_one_year <- function(year, statistic_df, province_data) {
    message("Processing ", year)
    input_files <- vapply(indicator_order, indicator_file, character(1), year = year)
    missing_files <- input_files[!file.exists(input_files)]
    if (length(missing_files) > 0) stop("Missing indicator allocation files: ", paste(missing_files, collapse = ", "))

    allocation_stack <- terra::rast(input_files)
    names(allocation_stack) <- indicator_order
    template <- allocation_stack[[1]]
    province <- province_data$province
    if (!identical(terra::crs(province), terra::crs(template))) {
        province <- terra::project(province, terra::crs(template))
    }
    province_r <- terra::rasterize(province, template, field = "province_code", touches = TRUE)

    allocation_mat <- terra::values(allocation_stack, mat = TRUE)
    province_vals <- as.integer(terra::values(province_r, mat = FALSE))
    valid_cell <- !is.na(province_vals) & rowSums(!is.na(allocation_mat)) == length(indicator_order)
    drip_hit <- allocation_mat == 3
    sprinkler_hit <- allocation_mat == 2
    drip_hit[is.na(drip_hit)] <- FALSE
    sprinkler_hit[is.na(sprinkler_hit)] <- FALSE
    drip_priority <- calc_priority(drip_hit)
    sprinkler_priority <- calc_priority(sprinkler_hit)

    out_vals <- rep(NA_integer_, nrow(allocation_mat))
    out_vals[valid_cell] <- 1L
    year_statistic <- statistic_df %>% filter(year == !!year)
    summary_rows <- vector("list", nrow(province_data$lookup))

    for (i in seq_len(nrow(province_data$lookup))) {
        province_code <- province_data$lookup$province_code[i]
        province_name <- province_data$lookup$province_name[i]
        province_cells <- which(valid_cell & province_vals == province_code)
        n_total <- length(province_cells)
        statistic_row <- year_statistic %>% filter(province_name == !!province_name)

        if (nrow(statistic_row) == 0) {
            drip_share <- NA_real_
            sprinkler_share <- NA_real_
            drip_target <- 0L
            sprinkler_target <- 0L
        } else {
            drip_share <- statistic_row$drip_share[1]
            sprinkler_share <- statistic_row$sprinkler_share[1]
            drip_target <- target_count(n_total, drip_share)
            sprinkler_target <- target_count(n_total, sprinkler_share)
        }

        drip_preferred <- province_cells[drip_priority[province_cells] > sprinkler_priority[province_cells]]
        drip_cells <- select_cells_with_supplement(drip_preferred, province_cells, drip_priority, drip_target)
        out_vals[drip_cells] <- 3L
        remaining_cells <- setdiff(province_cells, drip_cells)
        sprinkler_preferred <- remaining_cells[sprinkler_priority[remaining_cells] > drip_priority[remaining_cells]]
        sprinkler_cells <- select_cells_with_supplement(
            sprinkler_preferred, remaining_cells, sprinkler_priority, sprinkler_target
        )
        out_vals[sprinkler_cells] <- 2L

        summary_rows[[i]] <- data.frame(
            year = year, province_code = province_code, province_name = province_name,
            total_irrigated_pixels = n_total,
            surface_pixels = sum(out_vals[province_cells] == 1L, na.rm = TRUE),
            sprinkler_pixels = length(sprinkler_cells), micro_pixels = length(drip_cells)
        )
    }

    classified <- template
    terra::values(classified) <- out_vals
    irrigation_mask <- terra::rast(file.path(
        irrigation_mask_dir, sprintf("map_irr_rf_class_1km_%d_AEI65.tif", year)
    ))[[1]]
    if (!terra::compareGeom(template, irrigation_mask, stopOnError = FALSE)) {
        stop("Irrigation mask does not share the classification grid in ", year)
    }

    mask_vals <- terra::values(irrigation_mask, mat = FALSE)
    method_vals <- rep(NA_integer_, length(out_vals))
    method_vals[which(mask_vals == 0)] <- 0L
    irrigated_idx <- which(mask_vals == 1)
    method_vals[irrigated_idx] <- out_vals[irrigated_idx]
    method_vals[irrigated_idx[is.na(method_vals[irrigated_idx])]] <- 1L
    method <- template
    terra::values(method) <- method_vals
    names(method) <- sprintf("chn_irr_method_%d_1km", year)
    terra::writeRaster(
        method, file.path(output_dir, sprintf("chn_irr_method_%d_1km.tif", year)),
        overwrite = TRUE, datatype = "INT1U", NAflag = 255,
        gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=9", "TILED=YES")
    )

    consistency_vals <- rep(NA_real_, length(out_vals))
    drip_idx <- which(out_vals == 3L)
    sprinkler_idx <- which(out_vals == 2L)
    consistency_vals[drip_idx] <- drip_priority[drip_idx] / priority_max
    consistency_vals[sprinkler_idx] <- sprinkler_priority[sprinkler_idx] / priority_max
    consistency <- template
    terra::values(consistency) <- pmin(pmax(consistency_vals, 0), 1)
    names(consistency) <- sprintf("chn_consistency_%d_1km", year)
    terra::writeRaster(
        consistency, file.path(output_dir, sprintf("chn_consistency_%d_1km.tif", year)),
        overwrite = TRUE, datatype = "FLT4S", NAflag = -9999,
        gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=9", "TILED=YES")
    )
    bind_rows(summary_rows)
}

statistic_df <- prepare_statistic()
province_data <- prepare_province()
summary_df <- bind_rows(lapply(
    years, allocate_one_year,
    statistic_df = statistic_df,
    province_data = province_data
))
print(summary_df)
cat("Released irrigation-method and consistency products completed.\n")
