##### Dynamically extract fishbase Trophic Level #####
## Author: M. Schram
## Date: 2021-06-08
## Vers: 1.1 (2021-09-03)
##
## This code iteratively extracts the Trophic Level estimate of each SpCode in
## tbl_Taxa of the FishSurveys database. When SpCodes do not have species level
## idenfication (e.g. HASP) or a Trophic Level estimates is not available for
## that species, estimates are based on representative species## found in the
## Atlantic, Western Central FAO. The step-wise broadening of representative
## species continues through subfamily and family.
## 
library(rfishbase)
library(tidyverse)
library(dplyr)
library(readxl)
library(EnvStats)

# DietTroph = diet composition data estimation
# DietTLu   = diet composition data in unfished system
# FoodTroph = prey item (as opposed to % composition) using a randomization
#             method
# Chose FoodTroph because it allows equitable comparison across species instead
# of limited diet data

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

#If trophic level already included, this will remove
Sp_Tbl <- select(Sp_Tbl, c(SpCode:Class))

#### Species level Trophic Level estimate ####
sp_TrphL <- NULL
for (i in 1:length(Sp_Tbl$SpName)){
  temp <- ecology(Sp_Tbl$SpName[i]) %>%
    select(matches(c("Diet","Food")))
  temp <- add_column (temp, 
                      Species = Sp_Tbl$SpName[i], 
                      SpCode  = Sp_Tbl$SpCode[i], 
                      Common  = Sp_Tbl$Common[i])
  sp_TrphL <- bind_rows(sp_TrphL, temp)
}
rm(temp, i)

# Remove unreliastic estimates of trophic level, 'inferred' from 1 item, or
# 'derived' from basal resource + 1

sp_TrphL <- sp_TrphL[grep("inferred", sp_TrphL$FoodRemark, invert = TRUE),]
sp_TrphL <- sp_TrphL[grep("derived",  sp_TrphL$FoodRemark, invert = TRUE),]

Sp_Tbl <- Sp_Tbl %>%
  mutate(Trophic_Level = NaN, 
         .after = "Class") %>%
  mutate(Trophic_Level = replace(Trophic_Level, 
                                 SpCode %in% sp_TrphL$SpCode, 
                                 sp_TrphL$FoodTroph)) %>%
  arrange(SpCode)

Sp_Tbl$Trophic_Level[is.nan(Sp_Tbl$Trophic_Level)] <- NA

Spec <- filter(Sp_Tbl, !is.na(Trophic_Level))
Sub  <- filter(Sp_Tbl, is.na(Trophic_Level))

#### Genus level Trophic Level estimate ####
gn_TrphL <- NULL
for (i in 1:length(Sub$Genus)){
  temp_gn <- filter(distribution(species_list(Genus = local(Sub$Genus[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_gn$Species)){
    temp <- ecology(temp_gn$Species[j])%>%
      select(matches(c("Diet","Food")))
    temp <- bind_cols(temp, "Species" = temp_gn$Species[j])
    temp <- separate(temp, Species, c("Genus","Species"), sep = " ")
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- temp_sp[grep("inferred", temp_sp$FoodRemark, invert = TRUE),]
  temp_sp <- temp_sp[grep("derived",  temp_sp$FoodRemark, invert = TRUE),]
  temp_sp <- group_by(temp_sp, Genus) %>%
    summarize(Mean_TrphL = mean(temp_sp$FoodTroph, na.rm = TRUE),
              .groups    = "drop_last") %>%
    add_column(SpCode  = Sub$SpCode[i],
               Common  = Sub$Common[i],
               .before = "Genus") %>%
    add_column(Species = Sub$SpName[i],
               .before = "Mean_TrphL")
  gn_TrphL <- bind_rows(gn_TrphL, temp_sp)
}
rm(temp_gn, temp_sp, temp, i, j)

gn_TrphL$Mean_TrphL[is.nan(gn_TrphL$Mean_TrphL)] <- NA

Sp_Tbl <- Sp_Tbl %>%
  mutate(Trophic_Level = replace(Trophic_Level, 
                                 SpCode %in% gn_TrphL$SpCode, 
                                 gn_TrphL$Mean_TrphL)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Trophic_Level))
Sub  <- filter(Sp_Tbl, is.na(Trophic_Level))

#### Family level Trophic Level estimates ####
fa_TrphL <- NULL
for (i in 1:length(Sub$Family)){
  temp_fa <- filter(distribution(species_list(Family = local(Sub$Family[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_fa$Species)){
    temp <- ecology(temp_fa$Species[j])%>%
      select(matches(c("Diet","Food")))
    temp <- bind_cols(temp, "Species" = temp_fa$Species[j])
    temp <- separate(temp, 
                     Species, 
                     c("Genus","Species"), 
                     sep = " ")
    temp <- add_column(temp, SpCode = Sub$SpCode[i])
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- temp_sp[grep("inferred", temp_sp$FoodRemark, invert = TRUE),]
  temp_sp <- temp_sp[grep("derived",  temp_sp$FoodRemark, invert = TRUE),]
  temp_sp <- group_by(temp_sp, SpCode) %>%
    summarize(Mean_TrphL = mean(temp_sp$FoodTroph, na.rm = TRUE),
              .groups    = "drop_last") %>%
    add_column(Common  = Sub$Common[i],
               .before = "SpCode") %>%
    add_column(Species = Sub$SpName[i],
               .before = "Mean_TrphL")
  fa_TrphL <- bind_rows(fa_TrphL, temp_sp)
}

fa_TrphL$Mean_TrphL[is.nan(fa_TrphL$Mean_TrphL)] <- NA

rm(temp_fa, temp_sp, temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Trophic_Level = replace(Trophic_Level, 
                                 SpCode %in% fa_TrphL$SpCode, 
                                 fa_TrphL$Mean_TrphL)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Trophic_Level))
Sub  <- filter(Sp_Tbl, is.na(Trophic_Level))


#### Corrections for specific SpCodes ####

# BLSP - Blenny species
# Assigns the average Trophic Level for other blennys in the database
idx <- which(Sp_Tbl$SpCode == "BLSP")
temp <- filter(Sp_Tbl, grepl("Blenny",Common)) %>%
  filter(SpCode != "BLSP") %>%
  summarize(Trophic_Level = mean(Trophic_Level, na.rm = TRUE))
Sp_Tbl$Trophic_Level[idx] <- temp$Trophic_Level
rm(temp, idx)

# PLFM - Flounder species
# Assigns the average Trophic Level for other flounders in the database
idx <- which(Sp_Tbl$SpCode == "PLFM")
temp <- filter(Sp_Tbl, grepl("Flounder",Common)) %>%
  filter(SpCode != "PLFM") %>%
  summarize(Trophic_Level = mean(Trophic_Level, na.rm = TRUE))
Sp_Tbl$Trophic_Level[idx] <- temp$Trophic_Level
rm(temp, idx)

NO_coefs  <- filter(Sp_Tbl, is.na(Trophic_Level))

Spec <- filter(Sp_Tbl, !is.na(Trophic_Level))
Sub  <- filter(Sp_Tbl, is.na(Trophic_Level))


write.csv(Sp_Tbl, paste0("Outputs/TrophicLevel_", Sys.Date(), ".csv"))

#### Testing secondary run including derived/inferred data
sp_DietL <- NULL
for (i in 1:length(Sub$SpName)){
  temp <- ecology(Sub$SpName[i]) %>%
    select(matches(c("Diet","Food")))
  temp <- add_column (temp, 
                      Species = Sub$SpName[i], 
                      SpCode  = Sub$SpCode[i], 
                      Common  = Sub$Common[i])
  sp_DietL <- bind_rows(sp_DietL, temp)
}
rm(temp, i)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Trophic_Level = NaN, .after = "Class") %>%
  mutate(Trophic_Level = replace(Trophic_Level, 
                                 SpCode %in% sp_DietL$SpCode, 
                                 sp_DietL$FoodTroph)) %>%
  arrange(SpCode)

write.csv(Sp_Tbl, paste0("Outputs/TrophicLevel_DerInf_", Sys.Date(), ".csv"))
