
plot.world.distr <- function(sp_distr, sp_name, gift_shapes, odir,
                             world, worldcountries, eckertIV, bb, equator){
  
  # Function to plot the world distribution of the selected CWRs
  # based on the GIFT database
  
  sp_statuses <-
    sp_distr %>%
    mutate(native = ifelse(native == 1, "native", "non-native"),
           naturalized = ifelse(naturalized == 1, "naturalized",
                                "non-naturalized"),
           endemic_list = ifelse(endemic_list == 1, "endemic_list",
                                 "non-endemic_list")) %>%
    select(entity_ID, native, naturalized, endemic_list)
  
  # We rename the statuses based on the distinct combinations
  sp_statuses <-
    sp_statuses %>%
    mutate(Status = case_when(
      native == "native" & naturalized == "non-naturalized" ~ "native",
      native == "native" & is.na(naturalized) ~ "native",
      native == "non-native" & is.na(naturalized) ~ "non-native",
      native == "non-native" & naturalized == "naturalized" ~ "naturalized",
      native == "non-native" & naturalized == "non-naturalized" ~ "non-native",
      is.na(native) & is.na(naturalized) ~ "unknown"
    ))
  
  # Merge with the shapes
  sp_shape <- gift_shapes[which(gift_shapes$entity_ID %in% 
                                       unique(sp_distr$entity_ID)), ]
  sp_map <- left_join(sp_shape, sp_statuses, by = "entity_ID")
  
  sp_plot <-
    ggplot(world) +
    geom_sf(data  = bb, fill = "aliceblue") +
    geom_sf(data  = equator, color = "gray50", linetype = "dashed", linewidth = 0.1) +
    geom_sf(color = "gray70") +
    geom_sf(data  = world_countries, fill = "antiquewhite1", color = NA) +
    geom_sf(color = "gray50", linewidth = 0.1) +
    geom_sf(data  = bb, fill = NA) +
    geom_sf(data  = sp_map, color = "black", aes(fill = as.factor(Status))) +
    scale_fill_brewer("Status", palette = "Set2") +
    labs(title = expr(paste("Distribution map of ",
                            italic(!!sp_name),
                            " from GIFT v", !!GIFT_v)),
         subtitle = "Projection EckertIV") +
    coord_sf(crs = eckertIV) +
    theme_void(base_size = 7)
  
  ggsave(sp_plot,
         filename = paste0(sub(" ", "_", sp_name), ".png"),
         path = odir,
         bg = "white",
         width = 15,
         height = 10,
         units = "cm",
         dpi = 300)
  
}