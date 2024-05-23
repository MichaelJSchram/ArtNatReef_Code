##### Dynamically extract fishbase caudal fin aspect ratio #####
## Author: M. Schram
## Date: 2021-07-28
## Orig: 1.0 (2021-07-28)
## Vers: 1.1 (2021-09-03)
##
## This code iteratively extracts the caudal fin aspect ratio estimate of each
## SpCode in tbl_Taxa of the FishSurveys database. When SpCodes do not have
## species level idenfication (e.g. HASP) or a caudal fin aspect ratio estimates
## is not available for that species, estimates are based on representative
## species## found in the Atlantic, Western Central FAO. The step-wise
## broadening of representative species continues through subfamily and family.
## 
## 
library(rfishbase)
library(tidyverse)
library(dplyr)
library(readxl)
library(EnvStats)

##### Data import and correction for subsequent loops #####
Sp_Tbl <- read_excel("Data/Surveys_TaxaTable.xlsx")
Sp_Tbl <- mutate(Sp_Tbl, 
                 Genus       = replace(Genus,       is.na(Genus), 
                                       "Unknown"),
                 Subfamily   = replace(Subfamily,   is.na(Subfamily), 
                                       "Unknown"),
                 Family      = replace(Family,      is.na(Family), 
                                       "Unknown"),
                 Superfamily = replace(Superfamily, is.na(Superfamily), 
                                       "Unknown"),
                 Suborder    = replace(Suborder,    is.na(Suborder), 
                                       "Unknown"),
                 Order       = replace(Order,       is.na(Order), 
                                       "Unknown"),
                 Class       = replace(Class,       is.na(Class), 
                                       "Unknown"))

#If caudal fin aspect ratio already included, removes all values to re-query
Sp_Tbl <- select(Sp_Tbl, c(SpCode:Class))

#### Species level caudal fin aspect ratio ####
sp_caud <- NULL
for (i in 1:length(Sp_Tbl$SpName)){
  temp <- morphometrics(Sp_Tbl$SpName[i], 
                        fields = c("Species","AspectRatio"))
  temp <- group_by(temp, Species) %>%
    summarize(CaudalAspect = mean(AspectRatio),
              .groups      = "drop_last") %>%
    mutate(SpCode  = Sp_Tbl$SpCode[i],
           Common  = Sp_Tbl$Common[i],
           .before = "Species")
  sp_caud <- bind_rows(sp_caud, temp)
}
rm(temp, i)

Sp_Tbl <- Sp_Tbl %>%
  mutate(CaudalAR = NaN, 
         .after = "Class") %>%
  mutate(CaudalAR = replace(CaudalAR, 
                            SpCode %in% sp_caud$SpCode, 
                            sp_caud$CaudalAspect)) %>%
  arrange(SpCode)

Spec  <- filter(Sp_Tbl, !is.na(CaudalAR))
Sub   <- filter(Sp_Tbl, is.na(CaudalAR))

#### Genus level caudal fin aspect ratio ####
gn_caud <- NULL
for (i in 1:length(Sub$Genus)){
  temp_gn <- filter(distribution(species_list(Genus = local(Sub$Genus[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_gn$Species)){
    temp <- morphometrics(temp_gn$Species[j], 
                          fields = c("Species","AspectRatio"))
    temp <- group_by(temp, Species) %>%
      summarize(tempCaudAR = mean(AspectRatio, na.rm = TRUE),
                .groups    = "drop_last")
    temp <- separate(temp, 
                     Species, 
                     c("Genus","Species"), 
                     sep = " ")
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- group_by(temp_sp, Genus) %>%
    summarize(CaudalAspect = mean(tempCaudAR, na.rm = TRUE),
              .groups      = "drop_last") %>%
    mutate(SpCode  = Sub$SpCode[i],
           Common  = Sub$Common[i],
           Species = Sub$SpName[i],
           .before = "Genus")
  gn_caud <- bind_rows(gn_caud, temp_sp)
}
rm(temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(CaudalAR = replace(CaudalAR,
                            SpCode %in% gn_caud$SpCode,
                            gn_caud$CaudalAspect)) %>%
  arrange(SpCode)

Spec  <- filter(Sp_Tbl, !is.na(CaudalAR))
Sub   <- filter(Sp_Tbl, is.na(CaudalAR))

#### Subfamily level caudal fin aspect ratio ####
sf_caud <- NULL
for (i in 1:length(Sub$Subfamily)){
  temp_sf <- filter(distribution(species_list(Subfamily = local(Sub$Subfamily[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_sf$Species)){
    temp <- morphometrics(temp_sf$Species[j], 
                          fields = c("Species","AspectRatio"))
    temp <- group_by(temp, Species) %>%
      summarize(tempCaudAR = mean(AspectRatio, na.rm = TRUE),
                .groups    = "drop_last")
    temp <- separate(temp, 
                     Species, 
                     c("Genus","Species"), 
                     sep = " ")
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- mutate(temp_sp,
                    SpCode  = Sub$SpCode[i],
                    Common  = Sub$Common[i],
                    .before = "Genus")
  temp_sp <- group_by(temp_sp, SpCode) %>%
    summarize(CaudalAspect = mean(tempCaudAR, na.rm = TRUE),
              .groups      = "drop_last")
  sf_caud <- bind_rows(sf_caud, temp_sp)
}
rm(temp, i, j)


Sp_Tbl <- Sp_Tbl %>%
  mutate(CaudalAR = replace(CaudalAR, 
                            SpCode %in% sf_caud$SpCode, 
                            sf_caud$CaudalAspect)) %>%
  arrange(SpCode)

Spec  <- filter(Sp_Tbl, !is.na(CaudalAR))
Sub   <- filter(Sp_Tbl, is.na(CaudalAR))

#### Family level caudal fin aspect ratio ####
fa_caud <- NULL
for (i in 1:length(Sub$Family)){
  temp_fa <- filter(distribution(species_list(Family = local(Sub$Family[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_fa$Species)){
    temp <- morphometrics(temp_fa$Species[j], 
                          fields = c("Species","AspectRatio"))
    temp <- group_by(temp, Species) %>%
      summarize(tempCaudAR = mean(AspectRatio, na.rm = TRUE),
                .groups = "drop_last")
    temp <- separate(temp, 
                     Species, 
                     c("Genus","Species"), 
                     sep = " ")
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- mutate(temp_sp,
                    SpCode  = Sub$SpCode[i],
                    Common  = Sub$Common[i],
                    .before = "Genus")
  temp_sp <- group_by(temp_sp, SpCode) %>%
    summarize(CaudalAspect = mean(tempCaudAR, na.rm = TRUE),
              .groups      = "drop_last")
  fa_caud <- bind_rows(fa_caud, temp_sp)
}
rm(temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(CaudalAR = replace(CaudalAR,
                            SpCode %in% fa_caud$SpCode, 
                            fa_caud$CaudalAspect)) %>%
  arrange(SpCode)

Spec  <- filter(Sp_Tbl, !is.na(CaudalAR))
Sub   <- filter(Sp_Tbl, is.na(CaudalAR))

Sp_Tbl$CaudalAR[is.nan(Sp_Tbl$CaudalAR)] <- NA

#### Corrections for specific SpCodes ####

# BLSP - Blenny species
# Assigns the average Caudal Fin Aspect Ratio for other blennys in the database
idx <- which(Sp_Tbl$SpCode == "BLSP")
temp <- filter(Sp_Tbl, grepl("Blenny",Common)) %>%
  filter(SpCode != "BLSP") %>%
  summarize(CaudalAR = mean(CaudalAR, na.rm = TRUE))
Sp_Tbl$CaudalAR[idx] <- temp$CaudalAR
rm(temp, idx)

# PLFM - Flounder species
# Assigns the average Caudal Fin Aspect Ratio for other flounders in the database
idx <- which(Sp_Tbl$SpCode == "PLFM")
temp <- filter(Sp_Tbl, grepl("Flounder",Common)) %>%
  filter(SpCode != "PLFM") %>%
  summarize(CaudalAR = mean(CaudalAR, na.rm = TRUE))
Sp_Tbl$CaudalAR[idx] <- temp$CaudalAR
rm(temp, idx)

NO_CaudFinAR <- filter(Sp_Tbl, is.na(CaudalAR))

write.csv(Sp_Tbl, paste0("Outputs/CaudalAspect_", Sys.Date(), ".csv"))

# #### Genus level caudal fin aspect ratio w/o location restriction ####
# gn_caud <- NULL
# for (i in 1:length(Sub$Genus)){
#   temp_gn <- distribution(species_list(Genus = local(Sub$Genus[i])))
#   temp_gn <- distinct(temp_gn,)
#   temp_sp <- NULL
#   for (j in 1:length(temp_gn$Species)){
#     temp <- morphometrics(temp_gn$Species[j], fields = c("Species","AspectRatio"))
#     temp <- group_by(temp, Species) %>%
#       summarize(tempCaudAR = mean(AspectRatio, na.rm = TRUE),
#                 .groups = "drop_last")
#     temp <- separate(temp, Species, c("Genus","Species"), sep = " ")
#     temp_sp <- bind_rows (temp_sp, temp)
#   }
#   temp_sp <- group_by(temp_sp, Genus) %>%
#     summarize(CaudalAspect = mean(tempCaudAR, na.rm = TRUE),
#               .groups = "drop_last") %>%
#     mutate(SpCode = Sub$SpCode[i],
#            Common = Sub$Common[i],
#            Species = Sub$SpName[i],
#            .before = "Genus")
#   gn_caud <- bind_rows(gn_caud, temp_sp)
# }
# rm(temp, i, j)
# 
# Sp_Tbl <- Sp_Tbl %>%
#   mutate(CaudalAR = replace(CaudalAR, SpCode %in% gn_caud$SpCode, gn_caud$CaudalAspect)) %>%
#   arrange(SpCode)
# 
# Spec  <- filter(Sp_Tbl, !is.na(CaudalAR))
# Sub   <- filter(Sp_Tbl, is.na(CaudalAR))
# 
# #### Family level caudal fin aspect ratio w/o location restriction ####
# fa_caud <- NULL
# for (i in 1:length(Sub$Family)){
#   temp_gn <- distribution(species_list(Family = local(Sub$Family[i])))
#   temp_gn <- distinct(temp_gn,)
#   temp_sp <- NULL
#   for (j in 1:length(temp_gn$Species)){
#     temp <- morphometrics(temp_gn$Species[j], fields = c("Species","AspectRatio"))
#     temp <- group_by(temp, Species) %>%
#       summarize(tempCaudAR = mean(AspectRatio, na.rm = TRUE),
#                 .groups = "drop_last")
#     temp <- separate(temp, Species, c("Genus","Species"), sep = " ")
#     temp_sp <- bind_rows (temp_sp, temp)
#   }
#   temp_sp <- mutate(temp_sp,
#                     SpCode = Sub$SpCode[i],
#                     Common = Sub$Common[i],
#                     .before = "Genus")
#   temp_sp <- group_by(temp_sp, SpCode) %>%
#     summarize(CaudalAspect = mean(tempCaudAR, na.rm = TRUE),
#               .groups = "drop_last")
#   fa_caud <- bind_rows(fa_caud, temp_sp)
# }
# 
# Sp_Tbl <- Sp_Tbl %>%
#   mutate(CaudalAR = replace(CaudalAR, SpCode %in% fa_caud$SpCode, fa_caud$CaudalAspect)) %>%
#   arrange(SpCode)
# 
# Spec  <- filter(Sp_Tbl, !is.na(CaudalAR))
# Sub   <- filter(Sp_Tbl, is.na(CaudalAR))
