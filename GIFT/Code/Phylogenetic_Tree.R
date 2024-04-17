
#==============#
# Libraries ####
#==============#

library(GIFT)
library(data.table)
library(fst)
library(ape)           # pylogenetic tree operations
library(ggplot2)
library(BiocManager)   # needed to run the following packages?
library(ggtree)        # plotting phylogenetic trees
library(tidytree)      # plotting phylogenetic trees
library(ggtreeExtra)   # plotting phylogenetic trees
library(ape)           # pairwise distances from a phylogenetic tree
library(ggnewscale)
library(dplyr)
library(magrittr)


od <- # output directory
  "../Outputs/GIFT_Phylogeny" %T>%
  dir.create(recursive = T, showWarnings = F)


#=======================#
# Download phylogeny ####
#=======================#

# phy <- # phylogeny of vascular plants
#   GIFT_phylogeny(clade = "Tracheophyta",
#                  GIFT_version = "3.1") %T>%
#   saveRDS(file.path(od, "GIFT_Phylogeny.rds"))

phy <- # phylogeny of vascular plants
  readRDS(file.path(od, "GIFT_Phylogeny.rds"))

# Problem with hybrid species on the tip labels of the phylo tree
phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"] <-
  substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "x_"],
            3,
            nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))

phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"] <-
  substring(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"],
            3,
            nchar(phy$tip.label[substring(phy$tip.label, 1, 2) == "×_"]))


# Load crops and CWRs ####

cwr <- # crop wild relatives
  read.fst("../Outputs/GIFT_Harmonised_Names/CWR_GIFT_fuzzy_Spp_List.fst") %>%
  tibble

crops <- # crop species
  fread("../Outputs/GIFT_Harmonised_Names/Crops_GIFT_fuzzy_Spp_List.csv") %>%
  tibble

conflict <- # species listed both as crops and CWRs
  cwr[cwr$work_ID %in% crops$work_ID,] %>%
  arrange(work_ID) %>%
  select(work_species.y,
         work_genus,
         work_author,
         work_ID) %>%
  rename(work_species = work_species.y) %T>%
  fwrite(file.path(od, "Spp_both_crop_and_CWR.csv"))

cwr %<>% # remove species listed as crops from the CWR list
  filter(!(work_ID %in% conflict$work_ID))

cwr_sel <- 
  select(cwr, work_species.y, work_genus, work_ID) %>%
  rename(work_species = work_species.y)

crps_sel <-
  select(crops, harmonised_name, work_genus, work_ID) %>%
  rename(work_species = harmonised_name)

sp <- # list of species of interest
  bind_rows(cwr_sel,
            crps_sel) %>%
  mutate(is_crop = c(rep(F, nrow(cwr)), rep(T, nrow(crops)))) %>%
  filter(!duplicated(.[["work_ID"]])) # remove duplicates (same work_ID)

# sp_fam <- # family of each species
#   GIFT_taxgroup(work_ID = sp$work_ID,
#                 taxon_lvl = "family",
#                 GIFT_version = "3.1") %T>%
#   saveRDS(file.path(od, "GIFT_families.rds"))

sp_fam <- # family of each species
  readRDS(file.path(od, "GIFT_families.rds"))

sp %<>% # add family
  mutate(family = sp_fam) %>% 
  mutate(species = gsub("([[:punct:]])|\\s+",
                        "_",
                        work_species))

sp_phy <- # species with philogeny
  sp %>%
  filter(species %in% phy$tip.label)

phy <- # keep only crops and CWRs
  ape::keep.tip(phy = phy, # ape:: is needed as same function in tidytree{}
                tip = phy$tip.label[phy$tip.label %in% sp_phy$species])


#==================#
# Nearest crops ####
#==================#

dm <- # distance matrix
  cophenetic(phy)

sp_phy$nearest_crop <- # nearest crop
  sapply(seq_along(sp_phy$species), function(n){
    
    cr_candidates <- # crop candidates
      dm[sp_phy$species[n], sp_phy %>% filter(is_crop == T) %>% pull(species)]
    
    crop_dist <- # genetic distance
      cr_candidates %>%
      min
    
    crop_names <- # crop names (there may be multiple equidistant crops)
      cr_candidates[cr_candidates == crop_dist] %>%
      names
    
    return(sort(crop_names[1])) # return first in alfabetical order
    
  })

nok <- # next of keen crops for CWRs
  sp_phy %>%
  filter(is_crop == F) %>%
  pull(nearest_crop) %>%
  unique

sp_phy %<>% # keep only CWRs and next of keens
  filter(is_crop == F | species %in% nok)

sp_phy$crop_common_name <- # add crop's common name
  crops %>%
  filter(!duplicated(.[["harmonised_name"]])) %>% # remove duplicates (same harmonised_name)
  mutate(nearest_crop = gsub("([[:punct:]])|\\s+",
                             "_",
                             harmonised_name)) %>%
  left_join(sp_phy, ., by = "nearest_crop") %>%
  pull(common_name)


#===========#
# Plot ####
#===========#

phy <- # keep only crops and next of keens
  ape::keep.tip(phy = phy, # ape:: is needed as same function in tidytree{}
                tip = phy$tip.label[phy$tip.label %in% sp_phy$species])

sp_phy <- # sort species in the same order as phylogeny
  left_join(tibble(species = phy$tip.label),
            sp_phy,
            by = "species")

phy_plot <- # phylogenetic plot
  ggtree(phy, layout = "circular", colour = "grey") +
    geom_fruit(data = sp_phy,
               geom = geom_tile,
               mapping = aes(y = species,
                             fill = is_crop),
               width = 3,
               offset = 0.01) +
    geom_tiplab(size = 1.5,
                offset = 4) +
    geom_fruit(data = sp_phy,
               geom = geom_text,
               mapping = aes(y = species,
                             label = crop_common_name,
                             color = crop_common_name),
               size = 1.5,
               hjust = "outward",
               offset = 0.2,
               show.legend = F) +
  labs(title = "Phylogeny of European CWR and their nearest crops",
       subtitle = paste("The inner text circle contains the list of CWRs and wild crop varieties,",
                        "(these two are distiguished by the differently coloured flags).",
                        "The outer text circle indicates the nearest crop."))
  
ggsave(phy_plot,
       filename = "CRW_Tree_Next_Of_Keen.png",
       path = od,
       width = 40,
       height = 40,
       units = "cm",
       dpi = 300)
