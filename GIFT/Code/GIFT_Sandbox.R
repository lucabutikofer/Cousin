
#==============#
# Libraries ####
#==============#

library(GIFT)
library(sf)
library(fuzzyjoin)
library(fst)
library(data.table)
library(rnaturalearth)
library(rnaturalearthdata)
library(ape)             # pylogenetic tree operations
library(ggplot2)
library(BiocManager)   # needed to run the following packages?
library(ggtree)        # plotting phylogenetic trees
library(tidytree)      # plotting phylogenetic trees
library(ggtreeExtra)   # plotting phylogenetic trees
library(dplyr)
library(magrittr)

source("../Functions/Plot_World_Distr.R")


#=======================#
# Set-up environment ####
#=======================#

od <- # output directory
  "../Outputs/GIFT_Sandbox" %T>%
  dir.create(recursive = T, showWarnings = F)

# gift_shapes <- GIFT_shapes() # retrieves all shapefiles by default
# dir.create(file.path(od, "GIFT_Polygons"), showWarnings = F)
# st_write(gift_shapes, file.path(od, "GIFT_Polygons", "GIFT_Polygons.shp"))
gift_shapes <- st_read(file.path(od, "GIFT_Polygons", "GIFT_Polygons.shp"))

#==========================================#
# Retrieve CWR harmonised species names ####
#==========================================#

gsp <- function(x){ # function to retain only genus and species
  strsplit(x, " ") %>% first %>% `[`(1:2) %>% paste(collapse = " ")
}

cwrl <- # crop wild relatives list
  read.csv("../Inputs/20240226_EuropeanCWRchecklist.csv") %>%
  tibble %>%
  mutate(taxon = sapply(FullAcceptedName, gsp))

sp_list <- cwrl$taxon %>% `names<-`(NULL) %>% unique
gift_sp <- GIFT_species() %>% tibble

# Hard matches
foo <- sapply(sp_list, function(x) grep(x, gift_sp$work_species)) %>% unlist
bar <- gift_sp[as.numeric(foo),] %>% distinct

# With fuzzy matching with fuzzyjoin package
sp_list <- tibble(work_species = sp_list)

system.time({
  fuzz <-
    stringdist_join(sp_list, gift_sp,
                    by = "work_species",
                    mode = "left",
                    ignore_case = FALSE,
                    method = "jw",
                    max_dist = 99,
                    distance_col = "dist")
}) # 24 min.

cwr_hnms <- # CWR harmonised names
  fuzz %>%
  group_by(work_species.x) %>%
  slice_min(order_by = dist, n = 1)

write_fst(fuzz, file.path(od, "CWR_GIFT_fuzzy_matches.fst"))
write_fst(cwr_hnms, file.path(od, "CWR_GIFT_fuzzy_Spp_List.fst"))
cwr_hnms <- read_fst(file.path(od, "CWR_GIFT_fuzzy_Spp_List.fst")) %>% tibble


#===============================#
# Plot species distributions ####
#===============================#

GIFT_v <- "2.2"
# odp <- # output directory for plots
#   file.path("../Outputs/CWR_World_Maps") %T>%
#   dir.create(showWarnings = F)
# 
# # > Plotting aids ####
# 
# world <-
#   ne_coastline(scale = "medium", returnclass = "sf") %>%
#   st_wrap_dateline
# 
# world_countries <-
#   ne_countries(scale = "medium", returnclass = "sf") %>%
#   st_wrap_dateline
# 
# eckertIV <-
#   "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
# 
# bb <- # blue background
#   st_bbox(c(xmin = st_bbox(world)[["xmin"]],
#             xmax = st_bbox(world)[["xmax"]],
#             ymax = st_bbox(world)[["ymax"]],
#             ymin = st_bbox(world)[["ymin"]]),
#           crs = st_crs(4326)) %>%
#   st_make_grid(n = 100) %>%
#   st_union
# 
# equator <-
#   matrix(c(-180, 0, 180, 0), ncol = 2, byrow = TRUE) %>%
#   st_linestring %>%
#   st_sfc(crs = st_crs(world))
# 
# 
# # > Plot maps ####
# 
# pb <- txtProgressBar(max = nrow(cwr_hnms), style = 3) 
# 
# for (sp in 1:nrow(cwr_hnms)){ # for each CWR
#   
#   sp_name <-
#     cwr_hnms %>%
#     slice(sp) %>%
#     pull(work_species.x)
#   
#   sp_distr <-
#     GIFT_species_distribution(genus = strsplit(sp_name, " ")[[1]][1],
#                               epithet = strsplit(sp_name, " ")[[1]][2],
#                               aggregation = TRUE,
#                               GIFT_version = GIFT_v) %>%
#     select(entity_ID, native, naturalized, endemic_list)
#   
#   plot.world.distr(sp_distr, sp_name, gift_shapes, odir = odp,
#                    world, worldcountries, eckertIV, bb, equator)
#   
#   setTxtProgressBar(pb, sp)
#   
# }


#===================#
# Plot phylogeny ####
#===================#

odph <- # output directory for plots
  file.path("../Outputs/CWR_Phylogeny") %T>%
  dir.create(showWarnings = F)


# > Download phylogeny ####

phy <- GIFT_phylogeny(clade = "Tracheophyta", GIFT_version = "beta")
# Problem with hybrid species on the tip labels of the phylo tree
phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"] <-
  substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"],
            3,
            nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))

phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"] <-
  substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"],
            3,
            nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))

tax <- GIFT_taxonomy(GIFT_version = "beta")
gift_sp <- GIFT_species(GIFT_version = "beta")

gf <- GIFT_traits(trait_IDs = "1.2.1", agreement = 0.66, bias_ref = FALSE,
                  bias_deriv = FALSE, GIFT_version = "beta")

# Replacing space with _ for the species names
gf$work_species <- gsub(" ", "_", gf$work_species, fixed = TRUE)

# Retrieving family of each species
sp_fam <- GIFT_taxgroup(work_ID = unique(cwr_hnms$work_ID),
                        taxon_lvl = "family", GIFT_version = "beta")

sp_genus_fam <- # add family
  tibble(work_ID = unique(cwr_hnms$work_ID),
         work_species = unique(cwr_hnms$work_species.x),
         family = sp_fam)

sp_genus_fam <- # add genus
  left_join(sp_genus_fam,
            gift_sp[, c("work_ID", "work_genus")],
            by = "work_ID") %>%
  rename(genus = work_genus)


# > Add functional trait prevalence ####

sp_genus_fam <- # add functional trait
  left_join(sp_genus_fam,
            gf[, c("work_ID", "trait_value_1.2.1")],
            by = "work_ID")

genus_gf <- # compute prevalence of functional trait for each genus
  sp_genus_fam %>%
  group_by(genus) %>%
  mutate(prop_gf = round(100*sum(is.na(trait_value_1.2.1))/n(), 2)) %>%
  ungroup() %>%
  dplyr::select(-work_ID, -work_species, -family, -trait_value_1.2.1) %>%
  distinct(.keep_all = TRUE)

fam_gf <- # compute prevalence of functional trait for each family
  sp_genus_fam %>%
  group_by(family) %>%
  mutate(prop_gf = round(100*sum(is.na(trait_value_1.2.1))/n(), 2)) %>%
  ungroup() %>%
  dplyr::select(-work_ID, -work_species, -genus, -trait_value_1.2.1) %>%
  distinct(.keep_all = TRUE)

sp_genus_fam$species <- gsub("([[:punct:]])|\\s+", "_",
                             sp_genus_fam$work_species)

# Keeping one species with available phylogeny for each genus
one_sp_per_gen <- data.frame()
for(i in 1:n_distinct(sp_genus_fam$genus)){ # loop over genera
  
  focal_gen <- # focal genus
    unique(sp_genus_fam$genus)[i]
  
  gen_sp_i <- # All species in that genus
    sp_genus_fam[which(sp_genus_fam$genus == focal_gen),
                 "species"]
  
  gen_sp_i <- # Species from the genus available in the phylogeny
    gen_sp_i[gen_sp_i$species %in% phy$tip.label,]
  
  # Taking the first one (if at least one is available)
  if(ncol(gen_sp_i) > 0 & nrow(gen_sp_i) > 0){
    one_sp_per_gen <- rbind(one_sp_per_gen,
                            data.frame(species = gen_sp_i,
                                       genus = focal_gen))
  }
} 

one_sp_per_gen %<>% na.omit

one_sp_per_gen <- # Adding the trait coverage per genus
  left_join(one_sp_per_gen, genus_gf, by = "genus")

one_sp_per_gen <- # Adding the family
  left_join(one_sp_per_gen,
            sp_genus_fam[!duplicated(sp_genus_fam$genus),
                         c("genus", "family")],
            by = "genus") %>%
  rename(prop_gf_gen = prop_gf)

one_sp_per_gen <- # adding the trait coverage per family
  left_join(one_sp_per_gen, fam_gf, by = "family") %>%
  rename(prop_gf_fam = prop_gf) %>%
  tibble



# > Plot phylogenetic tree with functional trait ####

phy_gen <- # prune the tree at the genus level with ape package
  ape::keep.tip(phy = phy,
                tip = one_sp_per_gen$species)

phy_plot <- # phylogenetic plot
  ggtree(phy_gen,
         layout = "circular",
         color = "grey70") %<+%
  one_sp_per_gen +
  geom_fruit(geom = geom_tile,
             mapping = aes(fill = prop_gf_gen),
             width = 50,
             offset = 0.1) +
  geom_fruit(geom = geom_tile,
             mapping = aes(color = prop_gf_fam, fill = prop_gf_fam),
             width = 50,
             offset = 0.1,
             show.legend = FALSE) +
  scale_color_viridis_c() +
  scale_fill_viridis_c("Growth form availability per genus (%)") +
  # geom_text(aes(label = genus), color = 'red', size = 2, adj = "left") +
  geom_tiplab(size = 2, offset = 60) +
  theme(legend.text = element_text(size = 5),
        legend.title = element_text(size = 7),
        legend.key.size = unit(4, 'mm'),
        legend.position = "bottom",
        legend.box.spacing = unit(1, "cm"),
        plot.margin = unit(c(1,1,1,1),"cm"))

ggsave(phy_plot,
       filename = "CRW_Tree_All.png",
       path = odph,
       bg = "white",
       width = 30,
       height = 30,
       units = "cm",
       dpi = 300)







#___________________####
# FIRST ADAPTED CODE FROM GIFT TUTORIAL, TO BE MODIFIED ABOVE.
# 
# # > Download phylogeny ####
# 
# phy <- GIFT_phylogeny(clade = "Tracheophyta", GIFT_version = "beta")
# # Problem with hybrid species on the tip labels of the phylo tree
# phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"] <-
#   substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"],
#             3,
#             nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))
# 
# phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"] <-
#   substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"],
#             3,
#             nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))
# 
# tax <- GIFT_taxonomy(GIFT_version = "beta")
# gift_sp <- GIFT_species(GIFT_version = "beta")
# 
# gf <- GIFT_traits(trait_IDs = "1.2.1", agreement = 0.66, bias_ref = FALSE,
#                   bias_deriv = FALSE, GIFT_version = "beta")
# 
# # Replacing space with _ for the species names
# gf$work_species <- gsub(" ", "_", gf$work_species, fixed = TRUE)
# 
# # Retrieving family of each species
# sp_fam <- GIFT_taxgroup(work_ID = unique(cwr_hnms$work_ID),
#                         taxon_lvl = "family", GIFT_version = "beta")
# 
# sp_genus_fam <- # add family
#   tibble(work_ID = unique(cwr_hnms$work_ID),
#          work_species = unique(cwr_hnms$work_species.x),
#          family = sp_fam)
# 
# sp_genus_fam <- # add genus
#   left_join(sp_genus_fam,
#             gift_sp[, c("work_ID", "work_genus")],
#             by = "work_ID") %>%
#   rename(genus = work_genus)
# 
# 
# # > Add functional trait prevalence ####
# 
# sp_genus_fam <- # add functional trait
#   left_join(sp_genus_fam,
#             gf[, c("work_ID", "trait_value_1.2.1")],
#             by = "work_ID")
# 
# genus_gf <- # compute prevalence of functional trait for each genus
#   sp_genus_fam %>%
#   group_by(genus) %>%
#   mutate(prop_gf = round(100*sum(is.na(trait_value_1.2.1))/n(), 2)) %>%
#   ungroup() %>%
#   dplyr::select(-work_ID, -work_species, -family, -trait_value_1.2.1) %>%
#   distinct(.keep_all = TRUE)
# 
# fam_gf <- # compute prevalence of functional trait for each family
#   sp_genus_fam %>%
#   group_by(family) %>%
#   mutate(prop_gf = round(100*sum(is.na(trait_value_1.2.1))/n(), 2)) %>%
#   ungroup() %>%
#   dplyr::select(-work_ID, -work_species, -genus, -trait_value_1.2.1) %>%
#   distinct(.keep_all = TRUE)
# 
# sp_genus_fam$species <- gsub("([[:punct:]])|\\s+", "_",
#                              sp_genus_fam$work_species)
# 
# # Keeping one species with available phylogeny for each genus
# one_sp_per_gen <- data.frame()
# for(i in 1:n_distinct(sp_genus_fam$genus)){ # loop over genera
#   
#   focal_gen <- # focal genus
#     unique(sp_genus_fam$genus)[i]
#   
#   gen_sp_i <- # All species in that genus
#     sp_genus_fam[which(sp_genus_fam$genus == focal_gen),
#                            "species"]
#   
#   gen_sp_i <- # Species from the genus available in the phylogeny
#     gen_sp_i[gen_sp_i$species %in% phy$tip.label,]
#   
#   # Taking the first one (if at least one is available)
#   if(ncol(gen_sp_i) > 0){
#     gen_sp_i <- gen_sp_i[1,]
#     one_sp_per_gen <- rbind(one_sp_per_gen,
#                             data.frame(species = gen_sp_i,
#                                        genus = focal_gen))
#     }
#   } 
# 
# one_sp_per_gen %<>% na.omit
# 
# one_sp_per_gen <- # Adding the trait coverage per genus
#   left_join(one_sp_per_gen, genus_gf, by = "genus")
# 
# one_sp_per_gen <- # Adding the family
#   left_join(one_sp_per_gen,
#             sp_genus_fam[!duplicated(sp_genus_fam$genus),
#                          c("genus", "family")],
#             by = "genus") %>%
#   rename(prop_gf_gen = prop_gf)
# 
# one_sp_per_gen <- # adding the trait coverage per family
#   left_join(one_sp_per_gen, fam_gf, by = "family") %>%
#   rename(prop_gf_fam = prop_gf) %>%
#   tibble
# 
# 
# 
# # > Plot phylogenetic tree with functional trait ####
# 
# phy_gen <- # prune the tree at the genus level with ape package
#   ape::keep.tip(phy = phy,
#            tip = one_sp_per_gen$species)
# 
# phy_plot <- # phylogenetic plot
#   ggtree(phy_gen,
#        layout = "circular",
#        color = "grey70") %<+%
#   one_sp_per_gen +
#   geom_fruit(geom = geom_tile,
#              mapping = aes(fill = prop_gf_gen),
#              width = 50,
#              offset = 0.1) +
#   geom_fruit(geom = geom_tile,
#              mapping = aes(color = prop_gf_fam, fill = prop_gf_fam),
#              width = 50,
#              offset = 0.1,
#              show.legend = FALSE) +
#   scale_color_viridis_c() +
#   scale_fill_viridis_c("Growth form availability per genus (%)") +
#   # geom_text(aes(label = genus), color = 'red', size = 2, adj = "left") +
#   geom_tiplab(size = 2, offset = 60) +
#   theme(legend.text = element_text(size = 5),
#         legend.title = element_text(size = 7),
#         legend.key.size = unit(4, 'mm'),
#         legend.position = "bottom",
#         legend.box.spacing = unit(1, "cm"),
#         plot.margin = unit(c(1,1,1,1),"cm"))
# 
# ggsave(phy_plot,
#        filename = "CRW_Tree.png",
#        path = odph,
#        bg = "white",
#        width = 15,
#        height = 15,
#        units = "cm",
#        dpi = 300)
#___________________####  


#======================#
# Find related crop ####
#======================#

crps <- # crops list
  "../../Crops/FAO/Crop List_PlantTreaty_2022_9_5.csv" %>%
  read.csv %>%
  pull(Taxa) %>%
  sapply(. %>% strsplit(", ") %>% first) %>%
  unlist %>%
  `names<-`(NULL)
  
gift_sp <- GIFT_species() %>% tibble
sp_list <- tibble(work_species = crps)

system.time({
  fuzz_crop <-
    stringdist_join(sp_list, gift_sp,
                    by = "work_species",
                    mode = "left",
                    ignore_case = FALSE,
                    method = "jw",
                    max_dist = 99,
                    distance_col = "dist")
}) # 21 min.

crop_hnms <- # CWR harmonised names
  fuzz_crop %>%
  group_by(work_species.x) %>%
  slice_min(order_by = dist, n = 1)

fwrite(fuzz_crop, file.path(od, "Crops_GIFT_fuzzy_matches.csv"), nThread = 8)
fwrite(crop_hnms, file.path(od, "Crops_GIFT_fuzzy_Spp_List.csv"), nThread = 8)

