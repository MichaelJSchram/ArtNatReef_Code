##### Dynamically extract fishbase length-weight coefficients #####
## Author: M. Schram
## Date: 2021-05-04
## Original: (2021-05-04)
## Vers: 1.1 (2021-09-03)
##
## This code iteratively extracts the Max Length estimate of each SpCode in
## tbl_Taxa of the FishSurveys database. When SpCodes do not have species level
## idenfication (e.g. HASP) or a Max Length estimates is not available for that
## species, estimates are based on representative species## found in the
## Atlantic, Western Central FAO. The step-wise broadening of representative
## species continues through subfamily and family. at which time all remaining
## species lacking coefs must be manually examined.
## 
## 1.1 Update: cleaned up code, adjusted output format to data stamp

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

#If a-b coefs already included, this will remove
Sp_Tbl <- select(Sp_Tbl, c(SpCode:Class))

#### Species level MaxLength estimate ####
sp_MaxL <- NULL
for (i in 1:length(Sp_Tbl$SpName)){
    temp <- species(Sp_Tbl$SpName[i], field = c("Length","LengthFemale"))
    temp <- mutate(temp, 
                   Max_Length = as.numeric(pmax(Length, 
                                                LengthFemale, 
                                                na.rm = TRUE)))
    temp <- add_column(temp, 
                       Species = Sp_Tbl$SpName[i], 
                       SpCode  = Sp_Tbl$SpCode[i], 
                       Common  = Sp_Tbl$Common[i])
  sp_MaxL <- bind_rows(sp_MaxL, temp)
}
rm(temp, i)

Sp_Tbl <- Sp_Tbl %>%
  mutate(Max_Length = NaN, 
         .after    = "Class") %>%
  mutate(Max_Length = replace(Max_Length, 
                              SpCode %in% sp_MaxL$SpCode, 
                              sp_MaxL$Max_Length)) %>%
  arrange(SpCode)

Spec <- filter(Sp_Tbl, !is.na(Max_Length))
Sub  <- filter(Sp_Tbl, is.na(Max_Length))

#### Genus level MaxLength estimate ####
gn_MaxL <- NULL
for (i in 1:length(Sub$Genus)){
    temp_gn <- filter(distribution(species_list(Genus = local(Sub$Genus[i]))), 
                                   FAO == "Atlantic, Western Central")
    temp_sp <- NULL
    for (j in 1:length(temp_gn$Species)){
        temp <- species(temp_gn$Species[j], 
                        field = c("Length","LengthFemale","Genus"))
        temp <- mutate(temp, 
                       Max_Length = as.numeric(pmax(Length, 
                                                    LengthFemale, 
                                                    na.rm = TRUE)))
        temp_sp <- bind_rows (temp_sp, temp)
    }
    temp_sp <- group_by(temp_sp, Genus) %>%
        summarize(G.Mean_MaxL = geoMean(temp_sp$Max_Length, na.rm = TRUE),
                  A.Mean_MaxL = mean(temp_sp$Max_Length,    na.rm = TRUE),
                  .groups     = "drop_last") %>%
        add_column(SpCode  = Sub$SpCode[i],
                   Common  = Sub$Common[i],
                   .before = "Genus") %>%
        add_column(Species = Sub$SpName[i],
                   .before = "G.Mean_MaxL")
    gn_MaxL <- bind_rows(gn_MaxL, temp_sp)
}
gn_MaxL$G.Mean_MaxL[is.nan(gn_MaxL$G.Mean_MaxL)] <- NA
gn_MaxL$A.Mean_MaxL[is.nan(gn_MaxL$A.Mean_MaxL)] <- NA

rm(temp_gn, temp_sp, temp)

Sub <- Sub %>%
  mutate(Max_Length = replace(Max_Length, 
                              SpCode %in% gn_MaxL$SpCode, 
                              gn_MaxL$A.Mean_MaxL))
Sp_Tbl <- bind_rows(Spec, Sub) %>%
    arrange(SpCode)
rm(Spec, Sub, i, j)

Spec <- filter(Sp_Tbl, !is.na(Max_Length))
Sub  <- filter(Sp_Tbl, is.na(Max_Length))

#### Family level MaxLength estimates ####
fa_MaxL <- NULL
for (i in 1:length(Sub$Family)){
    temp_fa <- filter(distribution(species_list(Family = local(Sub$Family[i]))), 
                      FAO == "Atlantic, Western Central")
    temp_sp <- NULL
    for (j in 1:length(temp_fa$Species)){
        temp    <- species(temp_fa$Species[j], 
                           field = c("Length",
                                     "LengthFemale",
                                     "Genus",
                                     "Species"))
        temp    <- mutate(temp, 
                          Max_Length = as.numeric(pmax(Length, 
                                                       LengthFemale, 
                                                       na.rm = TRUE)))
        temp    <- mutate(temp, 
                          Family = Sub$Family[i])
        temp_sp <- bind_rows (temp_sp, temp)
    }
    temp_sp <- group_by(temp_sp, Family) %>%
      summarize(G.Mean_MaxL = geoMean(temp_sp$Max_Length, na.rm = TRUE),
                A.Mean_MaxL = mean(temp_sp$Max_Length,    na.rm = TRUE),
                .groups     = "drop_last") %>%
      add_column(SpCode  = Sub$SpCode[i],
                 Common  = Sub$Common[i],
                 .before = "Family") %>%
      add_column(Species = Sub$SpName[i],
                 .before = "G.Mean_MaxL")
    fa_MaxL <- bind_rows(fa_MaxL, temp_sp)
}
fa_MaxL$G.Mean_MaxL[is.nan(fa_MaxL$G.Mean_MaxL)] <- NA
fa_MaxL$A.Mean_MaxL[is.nan(fa_MaxL$A.Mean_MaxL)] <- NA
rm(temp_fa, temp_sp, temp)

Sub <- Sub %>%
  mutate(Max_Length = replace(Max_Length, 
                              SpCode %in% fa_MaxL$SpCode, 
                              fa_MaxL$A.Mean_MaxL))
Sp_Tbl <- bind_rows(Spec, Sub) %>%
  arrange(SpCode)
rm(Spec, Sub, i, j)

Spec <- filter(Sp_Tbl, !is.na(Max_Length))
Sub  <- filter(Sp_Tbl, is.na(Max_Length))

#### Corrections for specific SpCodes ####

# BLSP - Blenny species
# Assigns the average Max Length for other blennys in the database
idx <- which(Sp_Tbl$SpCode == "BLSP")
temp <- filter(Sp_Tbl, grepl("Blenny",Common)) %>%
    filter(SpCode != "BLSP") %>%
    summarize(Max_Length = mean(Max_Length, na.rm = TRUE))
Sp_Tbl$Max_Length[idx] <- temp$Max_Length
rm(temp, idx)

# PLFM - Flounder species
# Assigns the average Max Length for other flounders in the database
idx <- which(Sp_Tbl$SpCode == "PLFM")
temp <- filter(Sp_Tbl, grepl("Flounder",Common)) %>%
    filter(SpCode != "PLFM") %>%
    summarize(Max_Length = mean(Max_Length, na.rm = TRUE))
Sp_Tbl$Max_Length[idx] <- temp$Max_Length
rm(temp, idx)

NO_MaxLength  <- filter(Sp_Tbl, is.na(Max_Length))

write.csv(Sp_Tbl, paste0("Outputs/MaxLength_", Sys.Date(), ".csv"))
