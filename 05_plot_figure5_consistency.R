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

irr_file <- file.path(data_dir, sprintf("chn_irr_method_%d_1km.tif", plot_year))
consistency_file <- file.path(data_dir, sprintf("chn_consistency_%d_1km.tif", plot_year))
province_shp <- file.path(package_dir, "input_data", "boundaries", "province", "province.shp")


# ---------------------------- #
# Function
# ---------------------------- #
prepare_consistency_data <- function(irr_file, consistency_file) {
    input_files <- c(irr_file, consistency_file)
    missing_files <- input_files[!file.exists(input_files)]
    if (length(missing_files) > 0) {
        stop("Missing input raster files: ", paste(missing_files, collapse = ", "))
    }

    r <- c(terra::rast(irr_file)[[1]], terra::rast(consistency_file)[[1]])
    names(r) <- c("irr_type", "consistency")
    tile_size <- terra::res(r)

    df <- as.data.frame(r, xy = TRUE, na.rm = TRUE) %>%
        filter(irr_type %in% c(1, 2, 3)) %>%
        mutate(
            irr_label = factor(
                irr_type,
                levels = c(1, 2, 3),
                labels = c("Flood", "Sprinkler", "Drip")
            ),
            consistency = ifelse(irr_type %in% c(2, 3), consistency, NA_real_),
            flood_background = factor(ifelse(irr_type == 1, "Flood", NA_character_))
        )

    attr(df, "tile_width") <- tile_size[1]
    attr(df, "tile_height") <- tile_size[2]
    df
}

plot_consistency_map <- function(consistency_df, province_shp) {
    sf_use_s2(FALSE)

    province_map <- terra::vect(province_shp) %>%
        sf::st_as_sf()

    sea_bbox <- st_bbox(
        c(xmin = 107, xmax = 124, ymin = 2, ymax = 26),
        crs = st_crs(4326)
    )
    sea_map <- st_crop(province_map, sea_bbox)
    sea_df <- consistency_df %>%
        filter(x >= 107, x <= 124, y >= 2, y <= 26)
    flood_df <- consistency_df %>%
        filter(irr_label == "Flood")
    sea_flood_df <- sea_df %>%
        filter(irr_label == "Flood")
    wsi_df <- consistency_df %>%
        filter(irr_label %in% c("Drip", "Sprinkler"), is.finite(consistency))
    sea_wsi_df <- sea_df %>%
        filter(irr_label %in% c("Drip", "Sprinkler"), is.finite(consistency))
    tile_width <- attr(consistency_df, "tile_width")
    tile_height <- attr(consistency_df, "tile_height")

    p <- ggplot() +
        geom_tile(
            data = flood_df,
            aes(x = x, y = y),
            width = tile_width,
            height = tile_height,
            fill = "#bdbdbd"
        ) +
        geom_tile(
            data = wsi_df,
            aes(x = x, y = y, fill = consistency),
            width = tile_width,
            height = tile_height
        ) +
        geom_sf(
            data = province_map,
            fill = NA,
            color = "black",
            linewidth = 0.18
        ) +
        coord_sf(xlim = c(70, 146), ylim = c(18, 55), expand = FALSE) +
        scale_fill_gradientn(
            name = "Consistency",
            colours = c("#ffffb2", "#fecc5c", "#fd8d3c", "#bd0026"),
            limits = c(0, 1),
            breaks = seq(0, 1, 0.25),
            na.value = "transparent",
            guide = guide_colorbar(
                title.position = "right",
                title.hjust = 0.5,
                barwidth = unit(0.3, "cm"),
                barheight = unit(6, "cm")
            )
        ) +
        labs(x = NULL, y = NULL) +
        theme_void() +
        theme(
            legend.position = c(0.955, 0.42),
            legend.justification = c(0.5, 0.5),
            legend.title = element_text(size = 17, angle = 270),
            legend.text = element_text(size = 18),
            legend.direction = "vertical",
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_blank(),
            plot.background = element_rect(fill = "transparent", color = NA),
            plot.margin = margin(0, 0, 0, 0)
        )

    sea_p <- ggplot() +
        geom_tile(
            data = sea_flood_df,
            aes(x = x, y = y),
            width = tile_width,
            height = tile_height,
            fill = "#bdbdbd"
        ) +
        geom_tile(
            data = sea_wsi_df,
            aes(x = x, y = y, fill = consistency),
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
        scale_fill_gradientn(
            colours = c("#ffffb2", "#fecc5c", "#fd8d3c", "#bd0026"),
            limits = c(0, 1),
            na.value = "transparent",
            guide = "none"
        ) +
        theme_void() +
        theme(
            panel.background = element_rect(fill = "transparent", color = NA),
            plot.background = element_rect(fill = "transparent", color = NA)
        )

    ggdraw() +
        draw_plot(p, x = 0, y = 0, width = 1, height = 1) +
        draw_plot(sea_p, x = 0.63, y = 0.03, width = 0.25, height = 0.28) +
        theme(
            plot.margin = margin(0, 0, 0, 0),
            plot.background = element_rect(fill = "transparent", color = NA)
        )
}

province_abbrev <- function(code) {
    code <- as.character(code)
    dplyr::case_when(
        code == "110000" ~ "BJ",
        code == "120000" ~ "TJ",
        code == "130000" ~ "HE",
        code == "140000" ~ "SX",
        code == "150000" ~ "NM",
        code == "210000" ~ "LN",
        code == "220000" ~ "JL",
        code == "230000" ~ "HL",
        code == "310000" ~ "SH",
        code == "320000" ~ "JS",
        code == "330000" ~ "ZJ",
        code == "340000" ~ "AH",
        code == "350000" ~ "FJ",
        code == "360000" ~ "JX",
        code == "370000" ~ "SD",
        code == "410000" ~ "HA",
        code == "420000" ~ "HB",
        code == "430000" ~ "HN",
        code == "440000" ~ "GD",
        code == "450000" ~ "GX",
        code == "460000" ~ "HI",
        code == "500000" ~ "CQ",
        code == "510000" ~ "SC",
        code == "520000" ~ "GZ",
        code == "530000" ~ "YN",
        code == "540000" ~ "XZ",
        code == "610000" ~ "SN",
        code == "620000" ~ "GS",
        code == "630000" ~ "QH",
        code == "640000" ~ "NX",
        code == "650000" ~ "XJ",
        TRUE ~ code
    )
}

add_province_label <- function(consistency_df, province_shp) {
    province_v <- terra::vect(province_shp)
    province_attr <- province_v[, 2]
    names(province_attr) <- "province_code"

    point_v <- terra::vect(
        consistency_df,
        geom = c("x", "y"),
        crs = terra::crs(province_v)
    )

    province_df <- terra::extract(province_attr, point_v) %>%
        as.data.frame() %>%
        transmute(
            province_label = province_abbrev(province_code)
        )

    bind_cols(consistency_df, province_df)
}

prepare_type_consistency_bar_data <- function(consistency_df, province_shp, type_name, n_province = 8) {
    province_df <- add_province_label(consistency_df, province_shp) %>%
        filter(!is.na(province_label), irr_label == type_name)

    top_province <- province_df %>%
        count(province_label, name = "wsi_pixels") %>%
        slice_max(wsi_pixels, n = n_province, with_ties = FALSE) %>%
        arrange(desc(wsi_pixels))

    province_df %>%
        inner_join(top_province, by = "province_label") %>%
        group_by(province_label) %>%
        summarise(
            consistency = mean(consistency, na.rm = TRUE),
            pixels = n(),
            wsi_pixels = first(wsi_pixels),
            .groups = "drop"
        ) %>%
        mutate(
            province_label = factor(
                province_label,
                levels = top_province$province_label
            )
        )
}

prepare_national_consistency_data <- function(consistency_df) {
    consistency_df %>%
        filter(irr_label %in% c("Drip", "Sprinkler"), is.finite(consistency)) %>%
        group_by(irr_label) %>%
        summarise(consistency = mean(consistency, na.rm = TRUE), .groups = "drop") %>%
        mutate(irr_label = factor(irr_label, levels = c("Drip", "Sprinkler")))
}

plot_type_consistency_bar <- function(bar_df, national_df, type_name, type_colors, show_y_axis = TRUE) {
    national_value <- national_df %>%
        filter(irr_label == type_name) %>%
        pull(consistency)
    display_type <- irr_display[[type_name]]
    bar_df$display_type <- display_type

    p <- ggplot(bar_df, aes(x = province_label, y = consistency, fill = display_type)) +
        geom_col(
            width = 0.6,
            color = "black",
            linewidth = 0.18
        ) +
        geom_hline(
            yintercept = national_value,
            color = type_colors[[type_name]],
            linewidth = 0.75,
            linetype = "dashed"
        ) +
        scale_fill_manual(
            name = NULL,
            values = c(setNames(type_colors[[type_name]], display_type)),
            labels = c(display_type),
            drop = FALSE
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = seq(0, 1, 0.2),
            expand = expansion(mult = c(0, 0.02))
        ) +
        labs(x = NULL, y = "Consistency") +
        theme_classic() +
        theme(
            axis.text = element_text(size = 18, color = "black"),
            axis.title.y = element_text(size = 19, color = "black"),
            axis.line = element_line(linewidth = 0.35, color = "black"),
            axis.ticks = element_line(linewidth = 0.35, color = "black"),
            legend.position = c(0.98, 0.98),
            legend.justification = c(1, 1),
            legend.text = element_text(size = 17),
            legend.key.width = unit(0.55, "cm"),
            legend.key.height = unit(0.45, "cm"),
            legend.background = element_rect(fill = "transparent", color = NA),
            plot.margin = margin(8, 10, 8, 10)
        )

    if (!show_y_axis) {
        p <- p +
            labs(y = NULL) +
            theme(
                axis.text.y = element_blank(),
                axis.title.y = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.y = element_blank()
            )
    }

    p
}


# ---------------------------- #
# Build figure
# ---------------------------- #
type_colors <- c(
    Drip = "#2c7fb8",
    Sprinkler = "#f03b20"
)

irr_display <- c(
    Flood = "Surface",
    Drip = "Micro",
    Sprinkler = "Sprinkler"
)

consistency_df <- prepare_consistency_data(
    irr_file = irr_file,
    consistency_file = consistency_file
)

map_p <- plot_consistency_map(
    consistency_df = consistency_df,
    province_shp = province_shp
)
drip_bar_df <- prepare_type_consistency_bar_data(
    consistency_df = consistency_df,
    province_shp = province_shp,
    type_name = "Drip",
    n_province = 8
)
sprinkler_bar_df <- prepare_type_consistency_bar_data(
    consistency_df = consistency_df,
    province_shp = province_shp,
    type_name = "Sprinkler",
    n_province = 8
)
national_consistency_df <- prepare_national_consistency_data(
    consistency_df = consistency_df
)
drip_bar_p <- plot_type_consistency_bar(
    bar_df = drip_bar_df,
    national_df = national_consistency_df,
    type_name = "Drip",
    type_colors = type_colors,
    show_y_axis = TRUE
)
sprinkler_bar_p <- plot_type_consistency_bar(
    bar_df = sprinkler_bar_df,
    national_df = national_consistency_df,
    type_name = "Sprinkler",
    type_colors = type_colors,
    show_y_axis = FALSE
)

bar_row_p <- plot_grid(
    drip_bar_p,
    sprinkler_bar_p,
    nrow = 1,
    labels = c("b", "c"),
    label_size = 22,
    label_fontface = "bold",
    label_x = c(0.01, 0.01)
)

fig_p <- plot_grid(
    map_p,
    bar_row_p,
    ncol = 1,
    rel_heights = c(1, 0.62),
    labels = c("a", ""),
    label_size = 22,
    label_fontface = "bold",
    label_x = c(0.035, 0.01),
    label_y = c(1, 1.03)
)

out_file <- file.path(
    plot_dir,
    "Figure5_consistency.jpg"
)
ggsave(
    out_file,
    plot = fig_p,
    width = 10.2,
    height = 9.6,
    dpi = 600
)

cat("Consistency figure written:\n")
cat("  ", out_file, "\n", sep = "")
