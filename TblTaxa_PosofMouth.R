##### Dynamically extract fishbase Position of Mouth #####
## Author: M. Schram
## Date: 2021-07-29
## Vers: 1.1 (2021-09-03)
##
## This code iteratively extracts the Position of Mouth values of each SpCode in
## tbl_Taxa of the FishSurveys database. When SpCodes do not have species level
## idenfication (e.g. HASP) or a Position of Mouth values is not available for
## that species, estimates are based on representative species## found in the
## Atlantic, Western Central FAO. The step-wise broadening of representative
## species continues through subfamily and family.
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

#Removes Body_Shape_II if already included
Sp_Tbl <- select(Sp_Tbl, c(SpCode:Class))

# CABA causing problems, PosofMouth has 2 entries.
# Isolate and process individually, then append back onto table. 
Sp_CABA <- filter(Sp_Tbl, SpCode == "CABA")
Sp_Tbl  <- filter(Sp_Tbl, SpCode != "CABA")

#### Species level PosofMouth estimate ####
sp_Body <- NULL
for (i in 1:length(Sp_Tbl$SpName)){
  temp <- morphology(Sp_Tbl$SpName[i], field = c("Species","PosofMouth"))
  if (nrow(temp) > 1) {
    temp <- filter(temp, !is.na(PosofMouth))
  }
  temp    <- add_column(temp, 
                        SpCode = Sp_Tbl$SpCode[i], 
                        Common = Sp_Tbl$Common[i])
  sp_Body <- bind_rows(sp_Body, temp)
}
rm(temp, i)

Sp_Tbl <- Sp_Tbl %>%
  mutate(PosofMouth = NaN, 
         .after    = "Class") %>%
  mutate(PosofMouth = replace(PosofMouth, 
                              SpCode %in% sp_Body$SpCode,
                              sp_Body$PosofMouth)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(PosofMouth))
Sub  <- filter(Sp_Tbl, is.na(PosofMouth))

#### Genus level PosofMouth estimate ####
gn_Body <- NULL
for (i in 1:length(Sub$Genus)){
  temp_gn <- filter(distribution(species_list(Genus = local(Sub$Genus[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_gn$Species)){
    temp <- morphology(temp_gn$Species[j], 
                       field = c("Species","PosofMouth"))
    if (nrow(temp) > 1) {
      temp <- filter(temp, !is.na(PosofMouth))
    }
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- count(drop_na(temp_sp),PosofMouth)
  temp_sp <- select(filter(temp_sp, n == max(n)), PosofMouth)
  temp_sp <- add_column(temp_sp, 
                        SpCode = Sub$SpCode[i],
                        Common = Sub$Common[i],)
  gn_Body <- bind_rows(gn_Body, temp_sp)
}
rm(temp_gn, temp_sp, temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(PosofMouth = replace(PosofMouth, 
                              SpCode %in% gn_Body$SpCode, 
                              gn_Body$PosofMouth)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(PosofMouth))
Sub  <- filter(Sp_Tbl, is.na(PosofMouth))

#### Family level PosofMouth estimates ####
fa_Body <- NULL
for (i in 1:length(Sub$Family)){
  temp_fa <- filter(distribution(species_list(Family = local(Sub$Family[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_fa$Species)){
    temp <- morphology(temp_fa$Species[j], 
                       field = c("Species","PosofMouth"))
    if (nrow(temp) > 1){
      temp <- filter(temp, !is.na(PosofMouth))
    }
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- count(drop_na(temp_sp),PosofMouth)
  temp_sp <- select(filter(temp_sp, n == max(n)), PosofMouth)
  temp_sp <- add_column(temp_sp, 
                        SpCode = Sub$SpCode[i],
                        Common = Sub$Common[i],)
  fa_Body <- bind_rows(fa_Body, temp_sp)
}
rm(temp_fa, temp_sp, temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(PosofMouth = replace(PosofMouth, 
                              SpCode %in% fa_Body$SpCode, 
                              fa_Body$PosofMouth)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(PosofMouth))
Sub  <- filter(Sp_Tbl, is.na(PosofMouth))

#### Testing manual entry at Family Level ####
#### 2021-09-03: Desktop runs out of space using Family level query, will need
#### to revise approach or go with more storage space (laptop upgrade?)
fam_list <- read.csv("fam_list.csv")
fam_list <- "Echeneidae"
fa_TEST <- NULL
for (i in 1:length(fam_list)) {
  temp_fa <- filter(distribution(species_list(Family = local(fam_list[i]))))
  temp_sp <- NULL
  for (j in 1:length(temp_fa$Species)){
    temp <- morphology(temp_fa$Species[j], field = c("Species","PosofMouth"))
    if (nrow(temp) > 1) {
      temp <- filter(temp, !is.na(PosofMouth))
    }
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- distinct(temp_sp)
  temp_sp <- bind_cols(temp_sp, Family = fam_list[i])
  temp_sp <- count(temp_sp, PosofMouth, Family)
  fa_TEST <- bind_rows(fa_TEST, temp_sp)
}

#### Corrections for specific SpCodes ####

# BLSP - Blenny species
# Assigns the average Mouth Position for other blennys in the database
idx <- which(Sp_Tbl$SpCode == "BLSP")
temp <- filter(Sp_Tbl, grepl("Blenny",Common)) %>%
  filter(SpCode != "BLSP") %>%
  count(PosofMouth)
temp <- select(filter(temp, n == max(n)), PosofMouth)
Sp_Tbl$PosofMouth[idx] <- temp$PosofMouth
rm(temp, idx)

# PLFM - Flounder species
# Assigns the average Mouth Position for other flounders in the database
idx <- which(Sp_Tbl$SpCode == "PLFM")
temp <- filter(Sp_Tbl, grepl("Flounder",Common)) %>%
  filter(SpCode != "PLFM") %>%
  count(PosofMouth)
temp <- select(filter(temp, n == max(n)), PosofMouth)
Sp_Tbl$PosofMouth[idx] <- temp$PosofMouth
rm(temp, idx)

NO_PoM  <- filter(Sp_Tbl, is.na(PosofMouth))

write.csv(Sp_Tbl, paste0("Outputs/PosofMouth_", Sys.Date(), ".csv"))
