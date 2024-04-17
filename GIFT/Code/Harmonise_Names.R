
#==============#
# Libraries ####
#==============#

library(GIFT)
library(fuzzyjoin)
library(data.table)
library(dplyr)
library(magrittr)

od <- # output directory
  "../Outputs/GIFT_Harmonised_Names" %T>%
  dir.create(recursive = T, showWarnings = F)


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

# # Hard matches
# foo <- sapply(sp_list, function(x) grep(x, gift_sp$work_species)) %>% unlist
# bar <- gift_sp[as.numeric(foo),] %>% distinct

# With fuzzy matching with fuzzyjoin package
sp_list <- tibble(work_species = sp_list)

system.time({
  fuzz_cwr <-
    stringdist_join(sp_list, gift_sp,
                    by = "work_species",
                    mode = "left",
                    ignore_case = FALSE,
                    method = "jw",
                    max_dist = 99,
                    distance_col = "dist")
}) # 24 min.

cwr_hnms <- # CWR harmonised names
  fuzz_cwr %>%
  group_by(work_species.x) %>%
  slice_min(order_by = dist, n = 1)

write_fst(fuzz_cwr, file.path(od, "CWR_GIFT_fuzzy_matches.fst"))
write_fst(cwr_hnms, file.path(od, "CWR_GIFT_fuzzy_Spp_List.fst"))


#============================================#
# Retrieve crops harmonised species names ####
#============================================#

crps_raw <- # crops list
  "../../Crops/FAO/Crop List_PlantTreaty_2022_9_5.csv" %>%
  read.csv %>%
  tibble %>%
  select(Common.name...main, Taxa) %>%
  rename(Common_name = Common.name...main) %>%
  na.omit

crps <- # crops list
  crps_raw %>%
  pull(Taxa) %>%
  sapply(. %>% strsplit(", ") %>% first) %>%
  unlist %>%
  `names<-`(NULL) %>%
  unique

crps_name <-
  lapply(seq_along(crps), function(n){
    crps_raw %>%
      slice(grep(crps[n], crps_raw$Taxa)[1]) %>%
      pull(Common_name)
  }) %>%
  unlist


gift_sp <- GIFT_species(GIFT_version = "3.1") %>% tibble
sp_list <- tibble(work_species = crps,
                  common_name = crps_name)

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

crop_hrms <- # CWR harmonised names
  fuzz_crop %>%
  group_by(work_species.x) %>%
  slice_min(order_by = dist, n = 1) %>%
  rename(harmonised_name = work_species.y) %>%
  select(-work_species.x)

# fwrite(fuzz_crop, file.path(od, "Crops_GIFT_fuzzy_matches.csv"), nThread = 8)
fwrite(crop_hrms, file.path(od, "Crops_GIFT_fuzzy_Spp_List.csv"), nThread = 8)

