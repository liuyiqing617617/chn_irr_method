library(ggplot2)
library(dplyr)
library(cowplot)
library(sf)
library(terra)

rm(list = ls(all = TRUE))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), mustWork = TRUE)
package_dir <- dirname(dirname(script_file))
data_dir <- file.path(package_dir, "data")
plot_dir <- Sys.getenv("FIGSHARE_OUTPUT_DIR", file.path(package_dir, "figures"))
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------- #
# Path
# ---------------------------- #
plot_year <- 2020
stage_years <- c(2000, 2005, 2010, 2015, 2020)

irr_file <- file.path(data_dir, sprintf("chn_irr_method_%d_1km.tif", plot_year))
stage_files <- file.path(
    data_dir,
    sprintf("chn_irr_method_%d_1km.tif", stage_years)
)
province_shp <- file.path(package_dir, "input_data", "boundaries", "province", "province.shp")


# ---------------------------- #
# Function
# ---------------------------- #
province_abbrev <- function(code) {
    province_lookup <- c(
        "110000" = "BJ", "120000" = "TJ", "130000" = "HE",
        "140000" = "SX", "150000" = "NM", "210000" = "LN",
        "220000" = "JL", "230000" = "HL", "310000" = "SH",
        "320000" = "JS", "330000" = "ZJ", "340000" = "AH",
        "350000" = "FJ", "360000" = "JX", "370000" = "SD",
        "410000" = "HA", "420000" = "HB", "430000" = "HN",
        "440000" = "GD", "450000" = "GX", "460000" = "HI",
        "500000" = "CQ", "510000" = "SC", "520000" = "GZ",
        "530000" = "YN", "540000" = "XZ", "610000" = "SN",
        "620000" = "GS", "630000" = "QH", "640000" = "NX",
        "650000" = "XJ"
    )

    code <- as.character(code)
    label <- unname(province_lookup[code])
    label[is.na(label)] <- code[is.na(label)]
    label
}

prepare_irr_type_data <- function(irr_file) {
    if (!file.exists(irr_file)) {
        stop("Missing irrigation type raster: ", irr_file)
    }

    irr_r <- terra::rast(irr_file)[[1]]
    tile_size <- terra::res(irr_r)
    irr_df <- as.data.frame(irr_r, xy = TRUE, na.rm = TRUE)
    names(irr_df) <- c("x", "y", "irr_type")

    irr_df <- irr_df %>%
        filter(irr_type %in% c(1, 2, 3)) %>%
        mutate(
            irr_type = factor(
                irr_type,
                levels = c(1, 2, 3),
                labels = c("Surface", "Sprinkler", "Micro")
            )
        )

    attr(irr_df, "tile_width") <- tile_size[1]
    attr(irr_df, "tile_height") <- tile_size[2]
    irr_df
}

prepare_wsi_stage_data <- function(stage_files, stage_years) {
    missing_files <- stage_files[!file.exists(stage_files)]
    if (length(missing_files) > 0) {
        stop("Missing irrigation type rasters: ", paste(missing_files, collapse = ", "))
    }

    irr_r <- do.call(c, lapply(stage_files, function(x) terra::rast(x)[[1]]))
    names(irr_r) <- paste0("irr_", stage_years)
    tile_size <- terra::res(irr_r)
    irr_df <- as.data.frame(irr_r, xy = TRUE, na.rm = FALSE)

    irr_df <- irr_df %>%
        filter(irr_2020 %in% c(2, 3)) %>%
        mutate(
            wsi_stage = case_when(
                irr_2000 %in% c(2, 3) ~ "2000",
                irr_2005 %in% c(2, 3) ~ "2005",
                irr_2010 %in% c(2, 3) ~ "2010",
                irr_2015 %in% c(2, 3) ~ "2015",
                irr_2020 %in% c(2, 3) ~ "2020",
                TRUE ~ NA_character_
            ),
            wsi_stage = factor(
                wsi_stage,
                levels = c("2000", "2005", "2010", "2015", "2020")
            )
        )

    attr(irr_df, "tile_width") <- tile_size[1]
    attr(irr_df, "tile_height") <- tile_size[2]
    irr_df
}

plot_irr_type_map <- function(irr_df, province_shp, fill_col, fill_colors, legend_nrow = 1) {
    sf_use_s2(FALSE)

    province_map <- terra::vect(province_shp) %>%
        sf::st_as_sf()

    sea_bbox <- st_bbox(
        c(xmin = 107, xmax = 124, ymin = 2, ymax = 26),
        crs = st_crs(4326)
    )
    sea_map <- st_crop(province_map, sea_bbox)
    sea_df <- irr_df %>%
        filter(x >= 107, x <= 124, y >= 2, y <= 26)
    tile_width <- attr(irr_df, "tile_width")
    tile_height <- attr(irr_df, "tile_height")

    p <- ggplot() +
        geom_tile(
            data = irr_df,
            aes(x = x, y = y, fill = .data[[fill_col]]),
            width = tile_width,
            height = tile_height
        ) +
        geom_sf(
            data = province_map,
            fill = NA,
            color = "black",
            linewidth = 0.18
        ) +
        coord_sf(xlim = c(70, 140), ylim = c(18, 55), expand = FALSE) +
        scale_fill_manual(
            name = NULL,
            values = fill_colors,
            drop = FALSE
        ) +
        guides(fill = guide_legend(nrow = legend_nrow, byrow = TRUE)) +
        labs(
            title = NULL,
            x = NULL,
            y = NULL
        ) +
        theme_void() +
        theme(
            legend.position = c(0.0, 0.1),
            legend.justification = c(0, 0),
            legend.text = element_text(size = 14),
            legend.key.width = unit(0.36, "cm"),
            legend.key.height = unit(0.34, "cm"),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            plot.margin = margin(0, 0, 0, 0)
        )

    sea_p <- ggplot() +
        geom_tile(
            data = sea_df,
            aes(x = x, y = y, fill = .data[[fill_col]]),
            width = tile_width,
            height = tile_height
        ) +
        geom_sf(
            data = sea_map,
            fill = NA,
            color = "black",
            linewidth = 0.18
        ) +
        geom_rect(
            aes(xmin = 107, xmax = 124, ymin = 2, ymax = 26),
            fill = NA,
            color = "black",
            linewidth = 0.25
        ) +
        coord_sf(xlim = c(107, 124), ylim = c(2, 26), expand = FALSE) +
        scale_fill_manual(values = fill_colors, drop = FALSE, guide = "none") +
        theme_void() +
        theme(
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)
        )

    ggdraw() +
        draw_plot(p, x = -0.3, y = -0.22, width = 1.7, height = 1.5) +
        draw_plot(sea_p, x = 0.92, y = -0.15, width = 0.3, height = 0.4) +
        theme(
            plot.margin = margin(0, 0, 0, 0),
            plot.background = element_rect(fill = "transparent", color = NA)
        )
}

prepare_national_fraction <- function(stage_files, stage_years) {
    bind_rows(lapply(seq_along(stage_files), function(i) {
        freq_df <- terra::freq(terra::rast(stage_files[i])[[1]]) %>%
            as.data.frame() %>%
            filter(value %in% c(1, 2, 3)) %>%
            transmute(
                year = stage_years[i],
                irr_code = as.integer(value),
                pixels = count
            )

        data.frame(irr_code = c(1L, 2L, 3L)) %>%
            left_join(freq_df, by = "irr_code") %>%
            mutate(
                year = stage_years[i],
                pixels = ifelse(is.na(pixels), 0, pixels),
                fraction = pixels / sum(pixels) * 100,
                irr_type = factor(
                    irr_code,
                    levels = c(1, 2, 3),
                    labels = c("Surface", "Sprinkler", "Micro")
                )
            )
    }))
}

prepare_province_wsi_timeseries <- function(
    stage_files,
    stage_years,
    province_shp,
    area_per_pixel_kha = 0.1
) {
    province_v <- terra::vect(province_shp)
    province_code <- as.character(terra::values(province_v)[[2]])
    province_v$zone_id <- seq_len(nrow(province_v))

    template_r <- terra::rast(stage_files[1])[[1]]
    province_zone_r <- terra::rasterize(
        province_v,
        template_r,
        field = "zone_id"
    )
    province_lookup <- data.frame(
        zone_id = seq_len(nrow(province_v)),
        province_label = province_abbrev(province_code)
    )

    province_year_df <- bind_rows(lapply(seq_along(stage_files), function(i) {
        irr_r <- terra::rast(stage_files[i])[[1]]
        if (!terra::compareGeom(irr_r, template_r, stopOnError = FALSE)) {
            stop("Irrigation rasters do not share the same geometry: ", stage_files[i])
        }

        wsi_r <- terra::ifel(irr_r == 2 | irr_r == 3, 1, NA)
        terra::zonal(wsi_r, province_zone_r, fun = "sum", na.rm = TRUE) %>%
            as.data.frame() %>%
            setNames(c("zone_id", "pixels")) %>%
            left_join(province_lookup, by = "zone_id") %>%
            transmute(
                year = stage_years[i],
                province_label,
                pixels = ifelse(is.finite(pixels), pixels, 0),
                area_mha = pixels * area_per_pixel_kha / 1000
            )
    }))

    top_provinces <- province_year_df %>%
        filter(year == max(stage_years)) %>%
        arrange(desc(area_mha)) %>%
        slice_head(n = 10) %>%
        pull(province_label)

    province_year_df %>%
        filter(province_label %in% top_provinces) %>%
        mutate(
            province_label = factor(province_label, levels = top_provinces),
            year = factor(year, levels = stage_years)
        )
}

plot_fraction_trend <- function(fraction_df, type_colors) {
    ggplot(fraction_df, aes(x = year, y = fraction, color = irr_type, group = irr_type)) +
        geom_line(linewidth = 1.2) +
        geom_point(size = 3.8) +
        scale_color_manual(values = type_colors, guide = "none") +
        scale_x_continuous(breaks = stage_years) +
        scale_y_continuous(
            limits = c(0, 100),
            breaks = seq(0, 100, 20),
            expand = expansion(mult = c(0, 0.02))
        ) +
        labs(x = "Year", y = "Fraction (%)") +
        theme_classic() +
        theme(
            axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 16, color = "black"),
            axis.line = element_line(linewidth = 0.35, color = "black"),
            axis.ticks = element_line(linewidth = 0.35, color = "black"),
            plot.margin = margin(8, 20, 8, 20)
        )
}

plot_province_wsi_stack <- function(province_wsi_df) {
    province_colors <- setNames(
        grDevices::colorRampPalette(
            c("#053061", "#4393c3", "#756bb1", "#ef8a62", "#67001f")
        )(10),
        levels(province_wsi_df$province_label)
    )

    ggplot(province_wsi_df, aes(x = year, y = area_mha, fill = province_label)) +
        geom_col(width = 0.52) +
        scale_fill_manual(
            name = NULL,
            values = province_colors,
            breaks = levels(province_wsi_df$province_label),
            drop = FALSE
        ) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
        labs(
            x = "Year",
            y = expression(Area~"("*"×"~10^6*","~ha*")")
        ) +
        guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
        theme_classic() +
        theme(
            axis.text = element_text(size = 16, color = "black"),
            axis.title = element_text(size = 16, color = "black"),
            axis.line = element_line(linewidth = 0.35, color = "black"),
            axis.ticks = element_line(linewidth = 0.35, color = "black"),
            legend.position = c(0.025, 1.05),
            legend.justification = c(0, 1),
            legend.text = element_text(size = 16),
            legend.key.width = unit(0.42, "cm"),
            legend.key.height = unit(0.36, "cm"),
            legend.background = element_rect(
                fill = scales::alpha("white", 0.82),
                color = NA
            ),
            plot.margin = margin(8, 20, 8, 20)
        )
}


# ---------------------------- #
# Build figure
# ---------------------------- #
type_colors <- c(
    Surface = "#bdbdbd",
    Sprinkler = "#d73027",
    Micro = "#2166ac"
)

stage_colors <- c(
    "2000" = "#053061",
    "2005" = "#4393c3",
    "2010" = "#756bb1",
    "2015" = "#ef6548",
    "2020" = "#99000d"
)

irr_df <- prepare_irr_type_data(irr_file)
stage_df <- prepare_wsi_stage_data(stage_files, stage_years)
fraction_df <- prepare_national_fraction(stage_files, stage_years)
province_wsi_df <- prepare_province_wsi_timeseries(
    stage_files,
    stage_years,
    province_shp
)

irr_p <- plot_irr_type_map(
    irr_df = irr_df,
    province_shp = province_shp,
    fill_col = "irr_type",
    fill_colors = type_colors,
    legend_nrow = 2
)
stage_p <- plot_irr_type_map(
    irr_df = stage_df,
    province_shp = province_shp,
    fill_col = "wsi_stage",
    fill_colors = stage_colors,
    legend_nrow = 2
)
fraction_p <- plot_fraction_trend(fraction_df, type_colors)
province_wsi_p <- plot_province_wsi_stack(province_wsi_df)

# Isolate the statistical plots from the map row. Their internal margins can
# now be adjusted without resizing the maps above.
fraction_panel <- ggdraw() +
    draw_plot(fraction_p, x = 0, y = 0, width = 1, height = 1)
province_wsi_panel <- ggdraw() +
    draw_plot(province_wsi_p, x = 0, y = 0, width = 1, height = 1)

map_row <- plot_grid(
    irr_p,
    stage_p,
    ncol = 2,
    rel_widths = c(1, 1),
    scale = c(0.62, 0.62)
)

stat_row <- plot_grid(
    fraction_panel,
    province_wsi_panel,
    ncol = 2,
    rel_widths = c(1, 1)
)

make_panel_header <- function(label, title_line_1, title_line_2) {
    ggdraw() +
        draw_label(
            label,
            x = 0.035,
            y = 0.72,
            hjust = 0,
            vjust = 0.5,
            size = 18,
            fontface = "bold"
        ) +
        draw_label(
            title_line_1,
            x = 0.5,
            y = 0.72,
            hjust = 0.5,
            vjust = 0.5,
            size = 14,
            fontface = "bold"
        ) +
        draw_label(
            title_line_2,
            x = 0.5,
            y = 0.16,
            hjust = 0.5,
            vjust = 0.5,
            size = 14,
            fontface = "bold"
        )
}

map_title_row <- plot_grid(
    make_panel_header(
        "a",
        "Distribution of irrigation systems",
        "in 2020"
    ),
    make_panel_header(
        "b",
        "Timing of the transition to",
        "water-saving irrigation"
    ),
    ncol = 2,
    rel_widths = c(1, 1)
)

stat_title_row <- plot_grid(
    make_panel_header(
        "c",
        "Temporal changes in irrigation",
        "system fractions"
    ),
    make_panel_header(
        "d",
        "Water-saving irrigation area trends",
        "in the top 10 provinces"
    ),
    ncol = 2,
    rel_widths = c(1, 1)
)

map_block <- plot_grid(
    map_title_row,
    map_row,
    ncol = 1,
    rel_heights = c(0.11, 0.89)
)

stat_block <- plot_grid(
    stat_title_row,
    stat_row,
    ncol = 1,
    rel_heights = c(0.15, 0.85)
)

fig_p <- plot_grid(
    map_block,
    stat_block,
    ncol = 1,
    rel_heights = c(0.9, 0.72)
)

out_file <- file.path(
    plot_dir,
    "Figure3_irrigation_methods.jpg"
)
ggsave(
    out_file,
    plot = fig_p,
    width = 11,
    height = 7.85,
    dpi = 600
)

cat("Irrigation type map written:\n")
cat("  ", out_file, "\n", sep = "")
