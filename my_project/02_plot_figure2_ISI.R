library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(sf)
library(readxl)
library(terra)

rm(list = ls(all = TRUE))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
package_dir <- dirname(dirname(script_file))
data_dir <- file.path(package_dir, "data")
input_dir <- file.path(package_dir, "input_data")
plot_dir <- Sys.getenv("FIGSHARE_OUTPUT_DIR", file.path(package_dir, "figures"))
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------- #
# Map
# ---------------------------- #
plotChn <- function(
    dfMap = NULL,
    value = "value",
    legend_pos = "right",
    legend_inside = NULL,
    legend_title = NULL,
    legend_limit = c(1, 3),
    legend_color = NULL,
    legend_breaks = waiver()
) {
    sf_use_s2(FALSE)

    chn_shp <- file.path(input_dir, "boundaries", "province", "province.shp")
    chn_map <- terra::vect(chn_shp) %>%
        sf::st_as_sf()
    chn_map$province_code <- as.integer(chn_map[["省代码"]])

    dfMap$province_code <- as.integer(dfMap$province_code)
    plot_map <- chn_map %>%
        left_join(dfMap, by = "province_code")

    sea_bbox <- st_bbox(
        c(xmin = 107, xmax = 124, ymin = 2, ymax = 26),
        crs = st_crs(4326)
    )
    sea_map <- st_crop(plot_map, sea_bbox)

    p <- ggplot() +
        theme(
            panel.grid = element_blank(),
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.border = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank()
        ) +
        geom_sf(
            data = plot_map,
            aes_string(fill = value),
            color = "black",
            linewidth = 0.2
        ) +
        coord_sf(xlim = c(70, 140), ylim = c(18, 55), expand = FALSE) +
        labs(title = NULL, x = NULL, y = NULL) +
        scale_fill_gradientn(
            name = legend_title,
            limits = legend_limit,
            colors = legend_color,
            breaks = legend_breaks,
            oob = scales::squish,
            na.value = "grey90"
        ) +
        theme(
            legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            legend.key.width = unit(1.0, "cm"),
            legend.key.height = unit(0.25, "cm"),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_blank(),
            plot.margin = margin(0, 0, 0, 0)
        ) +
        guides(
            fill = guide_colorbar(
                direction = "horizontal",
                title.position = "bottom",
                title.hjust = 0.5,
                label.position = "bottom"
            )
        )

    if (is.null(legend_inside)) {
        p <- p +
            theme(
                legend.position = legend_pos,
                legend.justification = c(0.5, 0.5)
            )
    } else {
        p <- p +
            theme(
                legend.position = legend_inside,
                legend.justification = c(0, 0)
            )
    }

    sea_p <- ggplot() +
        geom_sf(
            data = sea_map,
            aes_string(fill = value),
            color = "black",
            linewidth = 0.2
        ) +
        geom_rect(
            aes(xmin = 107, xmax = 124, ymin = 2, ymax = 26),
            fill = NA,
            color = "black",
            linewidth = 0.25
        ) +
        scale_fill_gradientn(
            colors = legend_color,
            limits = legend_limit,
            oob = scales::squish,
            na.value = "grey90",
            guide = "none"
        ) +
        theme_void() +
        theme(
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)
        ) +
        coord_sf(xlim = c(107, 124), ylim = c(2, 26), expand = FALSE)

    combine_p <- ggdraw() +
        draw_plot(p, x = 0, y = 0, width = 1, height = 1) +
        draw_plot(sea_p, x = 0.73, y = 0.1, width = 0.25, height = 0.28) +
        theme(
            plot.margin = margin(0, 0, 0, 0),
            plot.background = element_rect(fill = "transparent", color = NA)
        )

    list(combined_plot = combine_p, main_plot = p)
}


# ---------------------------- #
# Prepare WSI data
# ---------------------------- #
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

prepare_wsi_map_data <- function(year = 2020) {
    statistic_file <- file.path(input_dir, "statistics", "chn_water_save_irr_statistics_province_level_mean.xlsx")
    province_shp <- file.path(input_dir, "boundaries", "province", "province.shp")

    statistic_df <- as.data.frame(read_excel(statistic_file)) %>%
        mutate(
            mapped_year = as.integer(mapped_year),
            province_name = standardize_name(province_name)
        ) %>%
        filter(mapped_year == year) %>%
        transmute(
            province_name,
            drip_frac = as.numeric(drip_frac),
            sprinkler_frac = as.numeric(sprinkler_frac)
        ) %>%
        mutate(
            flood_frac = 100 - drip_frac - sprinkler_frac,
            value = (3 * drip_frac + 2 * sprinkler_frac + flood_frac) / 100
        )

    province_map <- terra::vect(province_shp)
    province_df <- as.data.frame(province_map) %>%
        transmute(
            province_code = as.integer(`省代码`),
            province_name = standardize_name(`省`)
        ) %>%
        distinct()

    province_df %>%
        left_join(statistic_df, by = "province_name") %>%
        filter(!province_code %in% c(710000L, 810000L, 820000L)) %>%
        arrange(province_code)
}


# ---------------------------- #
# Prepare correlation data
# ---------------------------- #
indicator_dir <- file.path(input_dir, "indicators")
province_shp <- file.path(input_dir, "boundaries", "province", "province.shp")
statistic_file <- file.path(input_dir, "statistics", "chn_water_save_irr_statistics_province_level_mean.xlsx")
years <- c(2000, 2005, 2010, 2015, 2020)

indicator_info <- tibble::tibble(
    indicator = c(
        "IWU",
        "dET",
        "dKNDVI_dET_byCrop",
        "GS_delta_ET_by_IWU_groupmean_byCrop"
    ),
    indicator_label = c(
        "IWW",
        "deltaET",
        "deltaKNDVI/deltaET_byCrop",
        "deltaET/IWW_byCrop"
    ),
    province_statistic = rep("median", 4),
    dir = c(
        file.path(indicator_dir, "iwu"),
        file.path(indicator_dir, "det"),
        file.path(indicator_dir, "dkndvi_by_det"),
        file.path(indicator_dir, "det_by_iwu")
    ),
    file_template = c(
        "GS_IWU_1km_%d_05_09.tif",
        "GS_mean_delta_ET_AEI65_priority_1km_%d_05_09_filled.tif",
        "GS_mean_ratio_dkNDVI_dET_AEI65_priority_1km_%d_05_09_filled_byCrop_Q10.tif",
        "GS_delta_ET_by_IWU_1km_%d_05_09_byCrop_Q90.tif"
    )
)

significance_label <- function(p_value) {
    dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
    )
}

indicator_file <- function(indicator_row, year) {
    file.path(indicator_row$dir, sprintf(indicator_row$file_template, year))
}

prepare_province <- function() {
    province <- terra::vect(province_shp)
    province$province_code <- as.integer(province[["\u7701\u4ee3\u7801"]][, 1])
    province$province_name <- standardize_name(province[["\u7701"]][, 1])

    province_feature <- as.data.frame(province)[, c("province_code", "province_name")]
    province_lookup <- unique(province_feature)
    province_lookup <- province_lookup[!is.na(province_lookup$province_code), ]
    province_lookup <- province_lookup[!province_lookup$province_code %in% c(710000L, 810000L, 820000L), ]
    province_lookup <- province_lookup[order(province_lookup$province_code), ]

    list(province = province, feature = province_feature, lookup = province_lookup)
}

prepare_wsi_data <- function() {
    as.data.frame(read_excel(statistic_file)) %>%
        mutate(
            year = as.integer(mapped_year),
            province_name = standardize_name(province_name),
            drip_frac = as.numeric(drip_frac),
            sprinkler_frac = as.numeric(sprinkler_frac),
            flood_frac = 100 - drip_frac - sprinkler_frac,
            WSI = (3 * drip_frac + 2 * sprinkler_frac + flood_frac) / 100
        ) %>%
        filter(year %in% years) %>%
        select(year, province_name, drip_frac, sprinkler_frac, flood_frac, WSI)
}

province_zone_cache <- new.env(parent = emptyenv())

make_grid_key <- function(r) {
    paste(
        terra::crs(r),
        paste(terra::ext(r), collapse = ","),
        paste(terra::res(r), collapse = ","),
        terra::nrow(r),
        terra::ncol(r),
        sep = "|"
    )
}

get_province_zone <- function(r, province) {
    grid_key <- make_grid_key(r)
    if (exists(grid_key, envir = province_zone_cache, inherits = FALSE)) {
        return(get(grid_key, envir = province_zone_cache, inherits = FALSE))
    }

    province_year <- province
    if (!identical(terra::crs(province_year), terra::crs(r))) {
        province_year <- terra::project(province_year, terra::crs(r))
    }

    province_r <- terra::rasterize(
        province_year,
        r,
        field = "province_code",
        touches = TRUE
    )
    assign(grid_key, province_r, envir = province_zone_cache)

    province_r
}

calc_province_statistic <- function(
    raster_file,
    province,
    province_lookup,
    province_statistic
) {
    r <- terra::rast(raster_file)[[1]]
    province_r <- get_province_zone(r, province)

    r_vals <- terra::values(r, mat = FALSE)
    province_vals <- terra::values(province_r, mat = FALSE)
    valid <- is.finite(r_vals) & !is.na(province_vals)

    statistic_fun <- switch(
        province_statistic,
        median = function(x) median(x, na.rm = TRUE),
        q90 = function(x) unname(quantile(x, probs = 0.9, na.rm = TRUE, type = 7)),
        stop("Unsupported province statistic: ", province_statistic)
    )
    statistic_values <- tapply(
        r_vals[valid],
        province_vals[valid],
        statistic_fun
    )
    zonal_df <- data.frame(
        province_code = as.integer(names(statistic_values)),
        indicator_value = as.numeric(statistic_values),
        stringsAsFactors = FALSE
    )

    province_lookup %>%
        left_join(zonal_df, by = "province_code") %>%
        mutate(indicator_value = as.numeric(indicator_value))
}

build_indicator_statistic_data <- function(province, province_lookup) {
    rows <- list()
    row_i <- 1L

    for (i in seq_len(nrow(indicator_info))) {
        indicator_row <- indicator_info[i, ]

        for (year in years) {
            f <- indicator_file(indicator_row, year)
            if (!file.exists(f)) {
                stop("Missing input raster: ", f)
            }

            message(
                "Processing ", indicator_row$indicator, " ", year,
                " (", indicator_row$province_statistic, ")"
            )
            statistic_df <- calc_province_statistic(
                f,
                province,
                province_lookup,
                indicator_row$province_statistic
            ) %>%
                mutate(
                    year = year,
                    indicator = indicator_row$indicator,
                    indicator_label = indicator_row$indicator_label,
                    province_statistic = indicator_row$province_statistic,
                    input_file = f
                ) %>%
                select(
                    year, indicator, indicator_label, province_statistic,
                    province_code, province_name,
                    indicator_value, input_file
                )

            rows[[row_i]] <- statistic_df
            row_i <- row_i + 1L
        }
    }

    bind_rows(rows)
}

calc_one_cor <- function(df) {
    valid_df <- df %>%
        filter(is.finite(indicator_value), is.finite(WSI))

    if (nrow(valid_df) < 3) {
        return(tibble::tibble(
            n = nrow(valid_df),
            cor = NA_real_,
            p_value = NA_real_
        ))
    }

    test <- suppressWarnings(cor.test(
        valid_df$indicator_value,
        valid_df$WSI,
        method = "pearson"
    ))

    tibble::tibble(
        n = nrow(valid_df),
        cor = unname(test$estimate),
        p_value = test$p.value
    )
}

build_correlation_data <- function(join_df) {
    join_df %>%
        group_by(indicator, indicator_label, province_statistic) %>%
        group_modify(function(.x, .y) {
            pearson <- calc_one_cor(.x)

            tibble::tibble(
                n = pearson$n,
                pearson_cor = pearson$cor,
                pearson_p_value = pearson$p_value
            )
        }) %>%
        ungroup() %>%
        arrange(desc(abs(pearson_cor)))
}

build_delta_data <- function(join_df) {
    join_df %>%
        arrange(indicator, province_code, year) %>%
        group_by(
            indicator, indicator_label, province_statistic,
            province_code, province_name
        ) %>%
        mutate(
            year_start = year,
            year_end = lead(year),
            indicator_value_start = indicator_value,
            indicator_value_end = lead(indicator_value),
            WSI_start = WSI,
            WSI_end = lead(WSI),
            delta_indicator = indicator_value_end - indicator_value_start,
            delta_WSI = WSI_end - WSI_start
        ) %>%
        ungroup() %>%
        filter(!is.na(year_end)) %>%
        mutate(delta_period = paste(year_start, year_end, sep = "_")) %>%
        select(
            indicator, indicator_label, province_statistic,
            province_code, province_name,
            delta_period, year_start, year_end,
            indicator_value_start, indicator_value_end, delta_indicator,
            WSI_start, WSI_end, delta_WSI
        )
}

build_delta_correlation_data <- function(delta_df) {
    delta_df %>%
        group_by(indicator, indicator_label, province_statistic) %>%
        group_modify(function(.x, .y) {
            cor_df <- .x %>%
                transmute(
                    indicator_value = delta_indicator,
                    WSI = delta_WSI
                )
            pearson <- calc_one_cor(cor_df)

            tibble::tibble(
                n = pearson$n,
                pearson_cor = pearson$cor,
                pearson_p_value = pearson$p_value
            )
        }) %>%
        ungroup() %>%
        arrange(desc(abs(pearson_cor)))
}

prepare_corr_data <- function() {
    province_data <- prepare_province()
    wsi_data <- prepare_wsi_data()

    indicator_statistic <- build_indicator_statistic_data(
        province = province_data$province,
        province_lookup = province_data$lookup
    )

    indicator_wsi <- indicator_statistic %>%
        left_join(wsi_data, by = c("year", "province_name"))

    correlation_df <- build_correlation_data(indicator_wsi)
    delta_df <- build_delta_data(indicator_wsi)
    delta_correlation_df <- build_delta_correlation_data(delta_df)

    corr_data <- correlation_df %>%
        select(
            indicator, indicator_label, province_statistic,
            raw_n = n,
            raw_R = pearson_cor,
            raw_p_value = pearson_p_value
        ) %>%
        left_join(
            delta_correlation_df %>%
                select(
                    indicator,
                    delta_n = n,
                    delta_R = pearson_cor,
                    delta_p_value = pearson_p_value
                ),
            by = "indicator"
        ) %>%
        arrange(desc(raw_R)) %>%
        mutate(
            indicator_name = recode(
                indicator,
                IWU = "IWW",
                dET = "\u0394ET",
                dKNDVI_dET_byCrop = "\u0394kNDVI/\u0394ET",
                GS_delta_ET_by_IWU_groupmean_byCrop = "\u0394ET/IWW"
            ),
            raw_label = paste0(sprintf("%.2f", raw_R), significance_label(raw_p_value)),
            delta_label = paste0(sprintf("%.2f", delta_R), significance_label(delta_p_value))
        )

    corr_data$indicator_name <- factor(
        corr_data$indicator_name,
        levels = rev(corr_data$indicator_name)
    )

    cat("Raw R and within-province delta R Pearson correlation with WSI:\n")
    print(as.data.frame(corr_data), row.names = FALSE)

    corr_data
}

plot_corr_bar <- function(
    corr_data,
    value_col,
    label_col,
    title,
    x_limit,
    fill_color,
    show_y_text = TRUE
) {
    ggplot(
        corr_data,
        aes(x = .data[[value_col]], y = indicator_name)
    ) +
        geom_vline(xintercept = 0, color = "grey35", linewidth = 0.35) +
        geom_col(width = 0.58, fill = fill_color) +
        geom_text(
            aes(
                x = ifelse(.data[[value_col]] >= 0, .data[[value_col]] + 0.015, .data[[value_col]] - 0.015),
                label = .data[[label_col]],
                hjust = ifelse(.data[[value_col]] >= 0, 0, 1)
            ),
            size = 4.6
        ) +
        scale_x_continuous(
            limits = x_limit,
            breaks = pretty(x_limit, n = 4),
            expand = expansion(mult = c(0.02, 0.1))
        ) +
        coord_cartesian(clip = "off") +
        labs(
            title = title,
            x = "Pearson R",
            y = NULL
        ) +
        theme_classic(base_size = 14) +
        theme(
            plot.title = element_text(size = 14, hjust = 0.5),
            axis.text.y = if (show_y_text) element_text(size = 14, color = "black") else element_blank(),
            axis.text.x = element_text(size = 14, color = "black"),
            axis.title.x = element_text(size = 14),
            axis.line.y = element_blank(),
            axis.ticks.y = element_blank(),
            plot.margin = margin(5, 18, 5, 5)
        )
}


# ---------------------------- #
# Build figure
# ---------------------------- #
plot_data <- prepare_wsi_map_data(2020)
corr_data <- prepare_corr_data()
map_colors <- brewer.pal(9, "YlGnBu")

map_p <- plotChn(
    dfMap = plot_data,
    value = "value",
    legend_inside = c(0.05, 0.045),
    legend_title = "ISI 2020",
    legend_limit = c(1, 3),
    legend_color = map_colors,
    legend_breaks = seq(1, 3, by = 0.5)
)$combined_plot

raw_p <- plot_corr_bar(
    corr_data = corr_data,
    value_col = "raw_R",
    label_col = "raw_label",
    title = "Raw R",
    x_limit = c(-0.28, 0.7),
    fill_color = map_colors[6],
    show_y_text = TRUE
)

raw_p <- ggdraw() +
    draw_plot(raw_p, x = 0, y = 0.1, width = 1, height = 0.8)

final_p <- plot_grid(
    map_p,
    raw_p,
    nrow = 1,
    rel_widths = c(1, 0.6),
    labels = c("a", "b"),
    label_x = c(0.055, -0.03),
    label_y = c(0.965, 0.965),
    label_size = 16
)

ggsave(
    file.path(plot_dir, "Figure2_ISI.jpg"),
    plot = final_p,
    height = 5,
    width = 10.2,
    dpi = 600
)
