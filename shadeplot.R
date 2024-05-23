##### TITLE:      SHADE PLOT SOURCE CODE
#####
##### AUTHOR:     JILL THOMPSON-GRIM, M.S.
#####
##### OBJECTIVE:  Create shade plots to visualize species  
#####             composition and abundance data under
#####             mild data transformations (square-root,
#####             cube-root, and fourth-root) compared to
#####             raw untransformed data. Plots similar 
#####             to those created in Primer v7 (Clarke
#####             and Gorley 2014).
#####             
##### DATE:       7 APRIL 2021


shadeplot <- function(x) {
  
# INPUT -------------------------------------------------
# x = species data (as either dataframe or matrix) where,
#     columns = species 
#     rows    = numeric values (e.g., abundances)
#
# OUTPUT ------------------------------------------------
# Four panel plot of your data as:
#     (1) raw untransformed
#     (2) square-root transformed
#     (3) cube-root transformed
#     (4) fourth-root transformed data

### Load libraries
library(dplyr)              # for pipelines
library(tidyverse)          # for tidying data
library(ggplot2)            # for plotting
# library(plotly)             # for interactive plots
library(patchwork)

### Pivot the data for a heatmap 
  
# Convert to dataframe
df <- 
  as.data.frame(x) %>%
  # Retain the rownames
  rownames_to_column(var = "observation") %>%
  # Wrangle the data from a wide species matrix to a dataframe
  # with the 1st column as species and the 2nd column as counts
  pivot_longer(
    !observation, 
    names_to  = "species",
    values_to = "Abundance"
  ) %>%     
  # Set order for columns 
  select(order(as.character(colnames(.))))                          

### Extract factor order 
splevels <- unique(df$species)
oblevels <- unique(df$observation)

### Set base plot 
baseplot <- 
  ggplot(
    data = df %>% select(order(as.character(colnames(.)))),
    aes(
      y = fct_rev(factor(species,     levels = splevels)),
      x = fct_rev(factor(observation, levels = oblevels))
    )
  ) +
  labs(
    x = "Sample",
    y = "Taxon ID"
  ) +
  scale_fill_gradient(
    low  = "white", 
    high = "black"
  ) +
  theme_classic() +
  theme(
    plot.title         = element_text(hjust = 0.5),
    axis.text.x        = element_blank(),
    # turn off major interior lines
    panel.grid.major   = element_blank(),                                 
    # turn off minor interior lines
    panel.grid.minor   = element_blank(),                                 
    axis.line          = element_line(colour = "black"),
    legend.position    = "top"
  ) +
  scale_y_discrete(
    guide = guide_axis(n.dodge = 2)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5)
  )

### Make a heatmap of the data with RAW untransformed data 
sp_raw <- baseplot +
  geom_tile(
    aes(fill = Abundance)
  ) +
  labs(
    title = "Raw density",
    fill  = NULL
  )

### Make a heatmap of the data with square-root transformation 
sp_sqrt <- baseplot +
  geom_tile(
    aes(fill = Abundance^(1/2))
  ) +
  labs(
    title = "Square-root\ntransformed\ndensity",
    fill  = NULL
  ) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

### Make a heatmap of the data with cube-root transformation 
sp_cbrt <- baseplot +
  geom_tile(
    aes(fill = Abundance^(1/3))
  ) +
  labs(
    title = "Cube-root\ntransformed\ndensity",
    fill  = NULL
  ) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )

### Make a heatmap of the data with fourth-root transformation 
sp_4rt <- baseplot +
  geom_tile(
    aes(fill = Abundance^(1/4))
  ) +
  labs(
    title = "Fourth-root\ntransformed\ndensity",
    fill  = NULL
  ) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank()
  )


### Facet plots together to compare 
sp_options <- 
   sp_raw | sp_sqrt | sp_cbrt | sp_4rt

### View all the transformation shade plots 
sp_options
}