plot_basin_map <- function(
    basins_sf = NULL,
    rivers_sf = NULL,
    points_sf = NULL,
    main_stem_order = 2,
    add_scale = FALSE,
    add_north = FALSE
) {
  
  # ---------------------------------------------------------
  # 1. Determine reference CRS
  # ---------------------------------------------------------
  ref_crs <- NULL
  if (!is.null(basins_sf)) {
    ref_crs <- sf::st_crs(basins_sf)
  } else if (!is.null(rivers_sf)) {
    ref_crs <- sf::st_crs(rivers_sf)
  } else if (!is.null(points_sf)) {
    ref_crs <- sf::st_crs(points_sf)
  }
  
  # ---------------------------------------------------------
  # 2. CRS harmonization
  # ---------------------------------------------------------
  if (!is.null(basins_sf) && !identical(sf::st_crs(basins_sf), ref_crs)) {
    basins_sf <- sf::st_transform(basins_sf, ref_crs)
  }
  
  if (!is.null(rivers_sf) && !identical(sf::st_crs(rivers_sf), ref_crs)) {
    rivers_sf <- sf::st_transform(rivers_sf, ref_crs)
  }
  
  if (!is.null(points_sf) && !identical(sf::st_crs(points_sf), ref_crs)) {
    points_sf <- sf::st_transform(points_sf, ref_crs)
  }
  
  # Extract coordinates if points exist
  if (!is.null(points_sf)) {
    coords <- sf::st_coordinates(points_sf)
    points_df <- cbind(points_sf, coords)
    
    # Compute automatic nudges based on map extent
    dy <- diff(range(points_df$Y)) * 0.07   # vertical offset ≈ 2.5% of height
    dx <- diff(range(points_df$X)) * 0.05   # slight horizontal offset
  }
  
  # ---------------------------------------------------------
  # 3. Build plot
  # ---------------------------------------------------------
  p <- ggplot()
  
  # Basins
  if (!is.null(basins_sf)) {
    p <- p +
      geom_sf(
        data = basins_sf,
        fill = "grey95",
        color = "grey70",
        linewidth = 0.2
      )
  }
  
  # Rivers
  if (!is.null(rivers_sf)) {
    rivers_filt <- dplyr::filter(rivers_sf, strord > main_stem_order)
    
    p <- p +
      geom_sf(
        data = rivers_filt,
        aes(linewidth = strord),
        color = "steelblue4",
        alpha = 0.9
      ) +
      scale_linewidth(range = c(0.2, 1.5), guide = "none")
  }
  
  
  # Points + labels
  if (!is.null(points_sf)) {
    
    # point style: smaller, filled, dark
    p <- p +
      geom_label(
        data = points_df,
        aes(x = X, y = Y, label = name),
        size = 3,
        label.size = 0,                 # removes border
        fill = scales::alpha("white", 0.6),   # transparent background
        color = "black",
        label.r = unit(0.1, "lines"),
        nudge_y = dy,
        nudge_x = dx
      ) +
      geom_point(
        data = points_df,
        aes(x = X, y = Y),
        shape = 21,          # filled circle
        fill = "black",
        color = "white",
        size = 2.5,          # smaller
        stroke = 0.3
      ) 
   
  }
  
  # ---------------------------------------------------------
  # 4. Optional scale + north arrow
  # ---------------------------------------------------------
  if (add_scale) {
    p <- p +
      ggspatial::annotation_scale(
        location = "bl",
        width_hint = 0.3,
        text_cex = 0.75
      )
  }
  
  if (add_north) {
    p <- p +
      ggspatial::annotation_north_arrow(
        location = "tl",
        which_north = "true",
        style = ggspatial::north_arrow_fancy_orienteering,
        height = unit(1, "cm"),
        width  = unit(1, "cm")
      )
  }
  
  # ---------------------------------------------------------
  # 5. Final layout
  # ---------------------------------------------------------
  p +
    theme_bw() +
    labs(x = NULL, y = NULL) +
    coord_sf(expand = FALSE)
}
