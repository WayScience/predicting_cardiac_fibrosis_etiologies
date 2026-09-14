suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(platetools))


platemap_dir <- file.path("platemaps")
platemap_files <- list.files(platemap_dir, pattern = "_platemap\\.csv$", full.names = TRUE)
print(platemap_files)

output_fig_dir <- file.path("platemap_figures")
if (!dir.exists(output_fig_dir)) {
    dir.create(output_fig_dir)
}
platemap_suffix <- "_platemap_figure.png"

# Define output figure paths
output_platemap_files <- list()
for (platemap_file in platemap_files) {
    # Extract plate name and remove suffix 
    plate <- basename(platemap_file)
    plate <- stringr::str_remove(plate, "_platemap.csv") 
    
    output_platemap_files[[plate]] <- file.path(output_fig_dir, paste0(plate, platemap_suffix))
}

print(output_platemap_files)


# Load in all platemap CSV files
platemap_dfs <- list()
for (plate in names(output_platemap_files)) {
    # Find the umap file associated with the plate
    platemap_file <- platemap_files[stringr::str_detect(platemap_files, plate)]
    
    # Load in the umap data
    df <- readr::read_csv(
    platemap_file,
    col_types = readr::cols(.default = "c")
)

    # for plotting replace NaN as "healthy" for heart failure type
    df$heart_failure_type[is.na(df$heart_failure_type)] <- "Healthy"


    platemap_dfs[[plate]] <- df 
}

print(platemap_dfs)


# Consistent color mapping for heart failure subtype across all platemaps
heart_failure_type_colors <- c(
    "Healthy" = "#009E73",  # green (Okabe-Ito palette)
    "DCM" = "#D55E00",      # vermillion
    "ICM" = "#0072B2",      # blue
    "HLHS" = "#CC79A7"      # pink-purple
)

for (plate in names(platemap_dfs)) {
 {
    # output for each plate
    output_file <- output_platemap_files[[plate]]
    output_file <- paste0(output_file)
    
    platemap <-
        platetools::raw_map(
            data = platemap_dfs[[plate]]$heart_failure_type,
            well = platemap_dfs[[plate]]$well_position,
            plate = 96,
            size = 8
        ) +
        ggtitle(paste("Platemap layout for plate", plate)) +
        theme(plot.title = element_text(hjust=0.5, size = 10, face = "bold", margin = margin(b = -5))) +
        ggplot2::geom_point(aes(shape = platemap_dfs[[plate]]$cell_type)) +
        ggplot2::scale_shape_discrete(name = "Heart condition") +
        ggplot2::scale_fill_manual(name = "Heart failure subtype", values = heart_failure_type_colors) +
        theme(
            legend.position = "right",
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8),
        ) +
        guides(
        shape = guide_legend(override.aes = list(size = 2)),
        fill = guide_legend(override.aes = list(size = 5))
    )

    ggsave(
        output_file,
        platemap,
        dpi = 500,
        height = 3.5,
        width = 6
    )
    }
}

