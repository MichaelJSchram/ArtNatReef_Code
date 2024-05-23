##### Dynamically extract fishbase substrate type #####
## Author: M. Schram
## Date: 2021-08-13
## Vers: 1.1 (2021-09-08)
##
## This code iteratively extracts the Substrate Type estimate of each SpCode in
## tbl_Taxa of the FishSurveys database. When SpCodes do not have species level
## idenfication (e.g. HASP) or a Substrate estimates is not available for that
## species, estimates are based on representative species## found in the
## Atlantic, Western Central FAO. The step-wise broadening of representative
## species continues through subfamily and family. at which time all remaining
## species lacking coefs must be manually examined.
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

#Removes Substrate Type if already included
Sp_Tbl <- select(Sp_Tbl, c(SpCode:Class))


#### Species level Substrate_Type estimate ####
sp_Subs <- NULL
for (i in 1:length(Sp_Tbl$SpName)){
  temp    <- ecology(Sp_Tbl$SpName[i], 
                     field = c("Species","SoftBottom","HardBottom"))
  temp    <- add_column(temp, 
                        SpCode = Sp_Tbl$SpCode[i], 
                        Common = Sp_Tbl$Common[i])
  sp_Subs <- bind_rows(sp_Subs, temp)
}
rm(temp, i)

sp_Subs <- sp_Subs %>%
  mutate(Substrate_Type = NA) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  sp_Subs$SoftBottom < 0,
                                  "Soft Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  sp_Subs$HardBottom < 0,
                                  "Hard Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  sp_Subs$SoftBottom & sp_Subs$HardBottom < 0,
                                  "Mixed"))
  
Sp_Tbl <- Sp_Tbl %>%
  mutate(Substrate_Type = NA, 
         .after = "Class") %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  SpCode %in% sp_Subs$SpCode, 
                                  sp_Subs$Substrate_Type)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Substrate_Type))
Sub  <- filter(Sp_Tbl, is.na(Substrate_Type))

#### Genus level Substrate type ####
gn_Subs <- NULL
for (i in 1:length(Sub$Genus)){
  temp_gn <- filter(distribution(species_list(Genus = local(Sub$Genus[i]))),
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_gn$Species)){
    temp    <- ecology(temp_gn$Species[j], 
                       field = c("Species","SoftBottom","HardBottom"))
    temp_sp <- bind_rows (temp_sp, temp)
  }
  if (nrow(temp_sp > 0)) {
    temp_sp <- summarize(temp_sp, 
                         SoftBottom = mean(SoftBottom, na.rm = TRUE), 
                         HardBottom = mean(HardBottom, na.rm = TRUE)) %>%
      add_column(SpCode = Sub$SpCode[i],
                 Common = Sub$Common[i])
  } else {
    temp_sp <- cbind(temp_sp,
                     SpCode = Sub$SpCode[i],
                     Common = Sub$Common[i])
  }
  temp_sp$SoftBottom[is.nan(temp_sp$SoftBottom)] <- NA
  temp_sp$HardBottom[is.nan(temp_sp$HardBottom)] <- NA
  gn_Subs <- bind_rows(gn_Subs, temp_sp)
}
rm(temp_gn, temp_sp, temp, i, j)

gn_Subs <- gn_Subs %>%
  mutate(Substrate_Type = NA) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  gn_Subs$SoftBottom < 0,
                                  "Soft Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  gn_Subs$HardBottom < 0,
                                  "Hard Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  gn_Subs$SoftBottom & gn_Subs$HardBottom < 0,
                                  "Mixed"))

Sp_Tbl <- Sp_Tbl %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  SpCode %in% gn_Subs$SpCode, 
                                  gn_Subs$Substrate_Type)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Substrate_Type))
Sub  <- filter(Sp_Tbl, is.na(Substrate_Type))

#### Family level Substrate type ####
fa_Subs <- NULL
for (i in 1:length(Sub$Family)){
  temp_fa <- filter(distribution(species_list(Family = local(Sub$Family[i]))),
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_fa$Species)){
    temp    <- ecology(temp_fa$Species[j], 
                       field = c("Species","SoftBottom","HardBottom"))
    temp_sp <- bind_rows (temp_sp, temp)
  }
  if (nrow(temp_sp > 0)) {
    temp_sp <- summarize(temp_sp, 
                         SoftBottom = mean(SoftBottom, na.rm = TRUE), 
                         HardBottom = mean(HardBottom, na.rm = TRUE)) %>%
      add_column(SpCode = Sub$SpCode[i],
                 Common = Sub$Common[i])
  } else {
    temp_sp <- cbind(temp_sp,
                     SpCode = Sub$SpCode[i],
                     Common = Sub$Common[i])
  }
  temp_sp$SoftBottom[is.nan(temp_sp$SoftBottom)] <- NA
  temp_sp$HardBottom[is.nan(temp_sp$HardBottom)] <- NA
  fa_Subs <- bind_rows(fa_Subs, temp_sp)
}
rm(temp_fa, temp_sp, temp, i, j)

fa_Subs <- fa_Subs %>%
  mutate(Substrate_Type = NA) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  fa_Subs$SoftBottom < 0,
                                  "Soft Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  fa_Subs$HardBottom < 0,
                                  "Hard Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  fa_Subs$SoftBottom & fa_Subs$HardBottom < 0,
                                  "Mixed"))

Sp_Tbl <- Sp_Tbl %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  SpCode %in% fa_Subs$SpCode, 
                                  fa_Subs$Substrate_Type)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Substrate_Type))
Sub  <- filter(Sp_Tbl, is.na(Substrate_Type))

#### Corrections for specific SpCodes ####

# BLSP - Blenny species
# Assigns the average Body Shape II for other blennys in the database
idx <- which(Sp_Tbl$SpCode == "BLSP")
temp <- filter(Sp_Tbl, grepl("Blenny",Common)) %>%
  filter(SpCode != "BLSP") %>%
  count(Substrate_Type)
temp <- select(filter(temp, n == max(n)), Substrate_Type)
Sp_Tbl$Substrate_Type[idx] <- temp$Substrate_Type
rm(temp, idx)

# PLFM - Flounder species
# Assigns the average Body Shape II for other flounders in the database
idx <- which(Sp_Tbl$SpCode == "PLFM")
temp <- filter(Sp_Tbl, grepl("Flounder",Common)) %>%
  filter(SpCode != "PLFM") %>%
  count(Substrate_Type)
temp <- select(filter(temp, n == max(n)), Substrate_Type)
Sp_Tbl$Substrate_Type[idx] <- temp$Substrate_Type
rm(temp, idx)

NO_Subs  <- filter(Sp_Tbl, is.na(Substrate_Type))

write.csv(Sp_Tbl, paste0("Outputs/Substrate_Type_", Sys.Date(), ".csv"))

fa_Subs_ALL <- NULL
for (i in 1:length(NO_Subs$Family)){
  temp_fa <- species_list(Family = local(NO_Subs$Family[i]))
  temp_sp <- NULL
  for (j in 1:length(temp_fa)){
    temp    <- ecology(temp_fa[j], 
                       field = c("Species","SoftBottom","HardBottom"))
    temp_sp <- bind_rows (temp_sp, temp)
  }
  if (nrow(temp_sp > 0)) {
    temp_sp <- summarize(temp_sp, 
                         SoftBottom = mean(SoftBottom, na.rm = TRUE), 
                         HardBottom = mean(HardBottom, na.rm = TRUE)) %>%
      add_column(SpCode = NO_Subs$SpCode[i],
                 Common = NO_Subs$Common[i])
  } else {
    temp_sp <- cbind(temp_sp,
                     SpCode = NO_Subs$SpCode[i],
                     Common = NO_Subs$Common[i])
  }
  temp_sp$SoftBottom[is.nan(temp_sp$SoftBottom)] <- NA
  temp_sp$HardBottom[is.nan(temp_sp$HardBottom)] <- NA
  fa_Subs_ALL <- bind_rows(fa_Subs_ALL, temp_sp)
}
rm(temp_fa, temp_sp, temp, i, j)

fa_Subs_ALL <- fa_Subs_ALL %>%
  mutate(Substrate_Type = NA) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  fa_Subs_ALL$SoftBottom < 0,
                                  "Soft Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  fa_Subs_ALL$HardBottom < 0,
                                  "Hard Bottom")) %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  fa_Subs_ALL$SoftBottom & fa_Subs_ALL$HardBottom < 0,
                                  "Mixed"))

Sp_Tbl <- Sp_Tbl %>%
  mutate(Substrate_Type = replace(Substrate_Type, 
                                  SpCode %in% fa_Subs_ALL$SpCode, 
                                  fa_Subs_ALL$Substrate_Type)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Substrate_Type))
Sub  <- filter(Sp_Tbl, is.na(Substrate_Type))


write.csv(Sp_Tbl, paste0("Outputs/Substrate_Type_ALL_", Sys.Date(), ".csv"))
