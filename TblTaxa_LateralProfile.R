##### Dynamically extract fishbase Body Shape I #####
## Author: M. Schram
## Date: 2021-07-29
## Vers: 1.1 (2021-09-08)
##
## This code iteratively extracts the lateral profile values of each SpCode in
## tbl_Taxa of the FishSurveys database. When SpCodes do not have species level
## idenfication (e.g. HASP) or a lateral profile values is not available for
## that species, estimates are based on representative species## found in the
## Atlantic, Western Central FAO. The step-wise broadening of representative
## species continues through subfamily and family.
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

# Removes parse to taxonomic hierarchy only
Sp_Tbl <- select(Sp_Tbl, c(SpCode:Class))

#### Species level MaxLength estimate ####
sp_Body <- NULL
for (i in 1:length(Sp_Tbl$SpName)){
  temp <- morphology(Sp_Tbl$SpName[i], field = c("Species","BodyShapeI"))
  if (nrow(temp) > 1) {
    temp <- filter(temp, !is.na(BodyShapeI))
  }
  temp    <- add_column (temp, 
                         SpCode = Sp_Tbl$SpCode[i], 
                         Common = Sp_Tbl$Common[i])
  sp_Body <- bind_rows(sp_Body, temp)
}
rm(temp, i)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Lateral_Profile = NA, 
         .after = "Class") %>%
  mutate(Lateral_Profile = replace(Lateral_Profile, 
                              SpCode %in% sp_Body$SpCode, 
                              sp_Body$BodyShapeI)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Lateral_Profile))
Sub  <- filter(Sp_Tbl, is.na(Lateral_Profile))

#### Genus level MaxLength estimate ####
gn_Body <- NULL
for (i in 1:length(Sub$Genus)){
  temp_gn <- filter(distribution(species_list(Genus = local(Sub$Genus[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_gn$Species)){
    temp    <- morphology(temp_gn$Species[j], 
                       field = c("Species","BodyShapeI"))
    if (nrow(temp) > 1) {
      temp <- filter(temp, !is.na(BodyShapeI))
    }
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- count(temp_sp, BodyShapeI)
  temp_sp <- select(filter(temp_sp, n == max(n)), BodyShapeI)
  if (nrow(temp_sp) > 1) {
    temp_sp <- filter(temp_sp, !is.na(BodyShapeI))
  }
  temp_sp <- add_column(temp_sp, 
                        SpCode = Sub$SpCode[i],
                        Common = Sub$Common[i],)
  gn_Body <- bind_rows(gn_Body, temp_sp)
}
rm(temp_gn, temp_sp, temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Lateral_Profile = replace(Lateral_Profile, 
                              SpCode %in% gn_Body$SpCode, 
                              gn_Body$BodyShapeI)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Lateral_Profile))
Sub  <- filter(Sp_Tbl, is.na(Lateral_Profile))

#### Family level MaxLength estimates ####
fa_Body <- NULL
for (i in 1:length(Sub$Genus)){
  temp_fa <- filter(distribution(species_list(Family = local(Sub$Family[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_fa$Species)){
    temp <- morphology(temp_fa$Species[j],
                       field = c("Species","BodyShapeI"))
    if (nrow(temp) > 1) {
      temp <- filter(temp, !is.na(BodyShapeI))
    }
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- count(drop_na(temp_sp), BodyShapeI)
  temp_sp <- select(filter(temp_sp, n == max(n)), BodyShapeI)
  temp_sp <- add_column(temp_sp, 
                        SpCode = Sub$SpCode[i],
                        Common = Sub$Common[i],)
  fa_Body <- bind_rows(fa_Body, temp_sp)
}
rm(temp_fa, temp_sp, temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Lateral_Profile = replace(Lateral_Profile, 
                              SpCode %in% fa_Body$SpCode, 
                              fa_Body$BodyShapeI)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Lateral_Profile))
Sub  <- filter(Sp_Tbl, is.na(Lateral_Profile))


#### Corrections for specific SpCodes ####

# BLSP - Blenny species
# Assigns the average Lateral Profile for other blennys in the database
idx  <- which(Sp_Tbl$SpCode == "BLSP")
temp <- filter(Sp_Tbl, grepl("Blenny",Common)) %>%
  filter(SpCode != "BLSP") %>%
  count(Lateral_Profile)
temp <- select(filter(temp, n == max(n)), Lateral_Profile)
Sp_Tbl$Lateral_Profile[idx] <- temp$Lateral_Profile
rm(temp, idx)

# PLFM - Flounder species
# Assigns the average Lateral Profile for other flounders in the database
idx  <- which(Sp_Tbl$SpCode == "PLFM")
temp <- filter(Sp_Tbl, grepl("Flounder",Common)) %>%
  filter(SpCode != "PSFM") %>%
  count(Lateral_Profile)
temp <- select(filter(temp, n == max(n)), Lateral_Profile)
temp <- filter(temp, !is.na(Lateral_Profile))
Sp_Tbl$Lateral_Profile[idx] <- temp$Lateral_Profile
rm(temp, idx)

# Corrects two classifications that convey duplicate information for analyses
#
# "fusiform/normal"   --> fusiform
# "short and/or deep" --> short
Sp_Tbl$Lateral_Profile[grep("fusiform", Sp_Tbl$Lateral_Profile)] <- "fusiform"
Sp_Tbl$Lateral_Profile[grep("short"   , Sp_Tbl$Lateral_Profile)] <- "deep"

NO_LatPro  <- filter(Sp_Tbl, is.na(Lateral_Profile))

write.csv(Sp_Tbl, paste0("Outputs/LateralProfile_", Sys.Date(), ".csv"))
