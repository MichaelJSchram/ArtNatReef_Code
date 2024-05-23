##### Dynamically extract fishbase Cross Section #####
## Author: M. Schram
## Date: 2021-07-29
## Vers: 1.0 (2021-09-08)
##
## This code iteratively extracts the Cross Section values of each SpCode in
## tbl_Taxa of the FishSurveys database. When SpCodes do not have species level
## idenfication (e.g. HASP) or a Cross Section values is not available for that
## species, estimates are based on representative species## found in the
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

#### Species level BodyShapeII estimate ####
sp_Body <- NULL
for (i in 1:length(Sp_Tbl$SpName)){
  temp <- morphology(Sp_Tbl$SpName[i], field = c("Species","BodyShapeII"))
  if (nrow(temp) > 1) {
    temp <- filter(temp, !is.na(BodyShapeII))
  }
  temp <- add_column (temp, 
                      SpCode = Sp_Tbl$SpCode[i], 
                      Common = Sp_Tbl$Common[i])
  sp_Body <- bind_rows(sp_Body, temp)
}
rm(temp, i)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Cross_Section= NA, 
         .after = "Class") %>%
  mutate(Cross_Section = replace(Cross_Section, 
                                 SpCode %in% sp_Body$SpCode, 
                                 sp_Body$BodyShapeII)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Cross_Section))
Sub  <- filter(Sp_Tbl, is.na(Cross_Section))

#### Genus level BodyShapeII estimate ####
gn_Body <- NULL
for (i in 1:length(Sub$Genus)){
  temp_gn <- filter(distribution(species_list(Genus = local(Sub$Genus[i]))),
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_gn$Species)){
    temp    <- morphology(temp_gn$Species[j], 
                          field = c("Species","BodyShapeII"))
    if (nrow(temp) > 1) {
      temp <- filter(temp, !is.na(BodyShapeII))
    }
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- count(drop_na(temp_sp),BodyShapeII)
  temp_sp <- select(filter(temp_sp, n == max(n)), BodyShapeII)
  temp_sp <- add_column(temp_sp, 
                        SpCode = Sub$SpCode[i],
                        Common = Sub$Common[i],)
  gn_Body <- bind_rows(gn_Body, temp_sp)
}
rm(temp_gn, temp_sp, temp, i, j)

gn_SYPP <- filter(gn_Body, SpCode == "SYPP")
gn_TIMA <- filter(gn_Body, SpCode == "TIMA")
gn_SPHS <- filter(gn_Body, SpCode == "SPHS")
gn_Body <- filter(gn_Body, SpCode != "SYPP")
gn_Body <- filter(gn_Body, SpCode != "TIMA")
gn_Body <- filter(gn_Body, SpCode != "SPHS")

Sp_Tbl <- Sp_Tbl %>%
  mutate(Cross_Section = replace(Cross_Section, 
                                 SpCode %in% gn_Body$SpCode, 
                                 gn_Body$BodyShapeII)) %>%
  arrange(SpCode)

# Combines and assigns genera w/ two body shape classifications.
# 2021-09-08 needs further correction to self-reference likely representative
# from observations
    # gn_SYPP                   <- cbind(BodyShapeII = paste(gn_SYPP[1,1],
    #                                                        gn_SYPP[2,1],
    #                                                        sep = "/"),
    #                                    SpCode = gn_SYPP[1,2],
    #                                    Common = gn_SYPP[1,3])
    # idx                       <- which(Sp_Tbl$SpCode == "SYPP")
    # Sp_Tbl$Cross_Section[idx] <- gn_SYPP$BodyShapeII
    # rm(idx)
    # 
    # gn_TIMA                   <- cbind(BodyShapeII = paste(gn_TIMA[1,1],
    #                                                        gn_TIMA[2,1],
    #                                                        sep = "/"),
    #                                    SpCode = gn_TIMA[1,2],
    #                                    Common = gn_TIMA[1,3])
    # idx                       <- which(Sp_Tbl$SpCode == "TIMA")
    # Sp_Tbl$Cross_Section[idx] <- gn_TIMA$BodyShapeII
    # rm(idx)
    # 
    # gn_SPHS                   <- cbind(BodyShapeII = paste(gn_SPHS[1,1],
    #                                                        gn_SPHS[2,1],
    #                                                        sep = "/"),
    #                                    SpCode = gn_SPHS[1,2],
    #                                    Common = gn_SPHS[1,3])
    # idx                       <- which(Sp_Tbl$SpCode == "SPHS")
    # Sp_Tbl$Cross_Section[idx] <- gn_SPHS$BodyShapeII
    # rm(idx)

Spec <- filter(Sp_Tbl, !is.na(Cross_Section))
Sub  <- filter(Sp_Tbl, is.na(Cross_Section))

#### Family level BodyShapeII estimates ####
fa_Body <- NULL
for (i in 1:length(Sub$Genus)){
  temp_fa <- filter(distribution(species_list(Family = local(Sub$Family[i]))), 
                    FAO == "Atlantic, Western Central")
  temp_sp <- NULL
  for (j in 1:length(temp_fa$Species)){
    temp    <- morphology(temp_fa$Species[j], 
                          field = c("Species","BodyShapeII"))
    temp_sp <- bind_rows (temp_sp, temp)
  }
  temp_sp <- count(drop_na(temp_sp),BodyShapeII)
  temp_sp <- select(filter(temp_sp, n == max(n)), BodyShapeII)
  temp_sp <- add_column(temp_sp, 
                        SpCode = Sub$SpCode[i],
                        Common = Sub$Common[i],)
  fa_Body <- bind_rows(fa_Body, temp_sp)
}
rm(temp_fa, temp_sp, temp, i, j)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Cross_Section = replace(Cross_Section, 
                                 SpCode %in% fa_Body$SpCode, 
                                 fa_Body$BodyShapeII)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Cross_Section))
Sub  <- filter(Sp_Tbl, is.na(Cross_Section))


#### Corrections for specific SpCodes w/ ambiguity (e.g., Family = Unknown) ####

# BLSP - Blenny species
# Assigns the average Body Shape II for other blennys in the database
idx  <- which(Sp_Tbl$SpCode == "BLSP")
temp <- filter(Sp_Tbl, grepl("Blenny",Common)) %>%
  filter(SpCode != "BLSP") %>%
  count(Cross_Section)
temp <- select(filter(temp, n == max(n)), Cross_Section)
Sp_Tbl$Cross_Section[idx] <- temp$Cross_Section
rm(temp, idx)

# PLFM - Flounder species
# Assigns the average Body Shape II for other flounders in the database
idx  <- which(Sp_Tbl$SpCode == "PLFM")
temp <- filter(Sp_Tbl, grepl("Flounder",Common)) %>%
  filter(SpCode != "PLFM") %>%
  count(Cross_Section)
temp <- select(filter(temp, n == max(n)), Cross_Section)
Sp_Tbl$Cross_Section[idx] <- temp$Cross_Section
rm(temp, idx)

NO_BSII  <- filter(Sp_Tbl, is.na(Cross_Section))

write.csv(Sp_Tbl, paste0("Outputs/CrossSection_", Sys.Date(), ".csv"))