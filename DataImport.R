# Preloading and parameters -----

## Load in  packages -----

# Tidyverse
library(tidyverse)
library(readxl)
# Multivariate analytics
library(vegan)
library(BiodiversityR)
library(labdsv)
library(goeveg)
# Functional diversity package
library(FD)
# Package for pulling data directly from Fishbase
library(rfishbase)
# Beta diversity package for partitioning turnover vs. nestedness
library(betapart)
library(ape)
# Assist with plotting/visualizationsfor creating multi-panel plots
library(ggpubr)
library(patchwork)
library(ggallin)
library(ggrepel)
# Fast dummy coding (great for traits)
library(fastDummies)
# Misc
library(scales)

## Misc operators and parameters ----

# Establishes "not in" operator
`%!in%` <- Negate(`%in%`)
#Disables 'summarize()' regrouping... output
options(dplyr.summarise.inform = FALSE)

# Hyper-abundant, broad groups (e.g. genus level & larvae) with no novel terms
# of identification.
HypAbn  <- c("BAIT","LARV")
# One-off species codes that include single-count observations at family level
# (e.g. CAFM - Carangidae) when species level representation is far greater
# throughout the data, thus not constituting a unique observation.
Drop <- c("CAFM","CALF","CLFM","EPSP","GOFM","GYSP","HOSP",
          "LAFM","LUSP","OPSP","OSFM","SCPP","SPPP","TEFM","UNKN")

## Import site parameters -----
Sites <- 
  read_excel(
    "Data/Surveys_SiteDetails.xlsx"
  ) %>%
  mutate(
    Distance  = factor(Depth, levels = c("Deep", "Intermediate", "Shallow")),
    Depth     = NULL,
    Proximity = factor(Position, levels = c("North", "Middle", "South")),
    Position  = NULL,
    Type      = factor(Type, levels = c("Artificial", "Natural")),
    Series    = factor(Series, "levels" = c("OG8", "N6"))
  )

# Import and calculate density data -----
  
# Reads in the Excel file as exported from "Schram_SpeciesCount-by-Site" cross
# tab query in FishSurveys Access DB. NOTE: File must be renamed on export to
# match below, or change below to match exported file name
dat <- 
  as.data.frame(
    read_excel(
      "Data/Surveys_SpeciesCount-by-Site.xlsx"
    )
  ) %>%
  # Converts densities from character strings as exported from Access into
  # numeric values
  mutate_at(
    vars(
      -(SampleUnit:Survey_Radius)
    ), 
    as.numeric
  ) %>%
  # Singular survey performed out of season - should not be included in data. 
  filter(
    SampleUnit != "TI2_WI19_5_BB"
  ) %>%
  # Converts survey date to decimal value 
  mutate(
    J_Date = round(decimal_date(as.Date(parse_date_time(Cal_Date, "Ymd"))), 3),
    .after = Cal_Date
  ) %>%
  # NOTE: Some surveys on the same site/season were performed on different days
  # in early years. This only applies to 2 sites (AC5 and CWR) for a single
  # season (Spring 2016) and the differences in surveys dates were 4 days (AC5;
  # 2 x 2016/05/08 and 4 x 2016/05/012) and 5 days (CWR, 2 x 2016/05/09 + 2 x
  # 2016/05/14). The short time spans are not enough to cause concern in terms
  # of temporal similarity, thus finding a 'mean date' was deemed adequate to
  # compile site x season x year averages. J_Date and Cal_Date are equivalent in
  # ALL OTHER INSTANCES.
  group_by(
    Year, 
    Season, 
    Site
  ) %>%
  mutate(
    J_Date = 
      mean(J_Date)
  ) %>%
  ungroup() %>%
  # Calculates survey radius and subsequent taxon densities for each individual
  # survey
  mutate(
    Survey_Area = pi * (Survey_Radius^2), 
    .after      = Survey_Radius
  ) %>%
  # Expands from a Species by Site matrix for density estimates
  pivot_longer(
    -c(SampleUnit:Survey_Area),
    names_to  = "SpCode",
    values_to = "Count"
  ) %>%
  mutate(
    Density = Count/Survey_Area
  ) %>%
  # Estimates mean density based on all surveys performed on each site, each
  # season, each year.
  group_by(
    Year, 
    Season, 
    J_Date, 
    Site, SpCode
  ) %>%
  summarize(
    Density       = mean(Density),
    Survey_Radius = mean(Survey_Radius),
    Survey_Area   = mean(Survey_Area),
    Survey_Count  = n()
  ) %>%
  ungroup() %>%
  # Expands back into a Species by Site matrix
  pivot_wider(
    names_from  = SpCode, 
    values_from = Density
  ) %>%
  # Pulls in site details
  merge(
    Sites,
    by.x = "Site",
    by.y = "SiteName"
  ) %>%
  # Relocates specific columns
  relocate(
    c(Site, Survey_Count, Survey_Radius, Survey_Area, 
      Distance, Proximity, Series, Type), 
    .after = J_Date
  ) %>%
  # Sorts order of data
  arrange(
    Year,
    match(
      Season, 
      c("Winter", "Spring", "Summer", "Fall")
    ),
    match(
      Series, 
      c("OG8", "N6")
    )
  )

## Establish "Reefs" variable -----
Reefs <- 
  filter(
    dat,
    Series == "OG8"
  ) %>%
  # Assigns factor to order data by collection year & season
  # 2, 5, 8, & 11 used as midpoint of each season (month)
  mutate(
    Month = 
      case_when(
        Season == "Winter" ~ "02",
        Season == "Spring" ~ "05",
        Season == "Summer" ~ "08",
        Season == "Fall"   ~ "11"
      ),
    Year = 
      factor(Year),
    Season =
      factor(Season, levels = c("Winter", "Spring", "Summer", "Fall"))
  ) %>%
  # Establishes a YYYY-MM character vector from Year and Month
  unite(
    Date,
    c(Year, Month),
    sep    = "-",
    remove = FALSE
  ) %>%
  mutate(
  # Converts YYYY-MM character vector into a YM-date value for temporal analyses
    Date =
      decimal_date(
        as.Date(
          parse_date_time(
            Date,
            orders = "ym"
          )
        )
      )
  ) %>%
  relocate(
    c(Month, Date), 
    .before = Site
  ) %>%
  # Exclude specific taxa as defined above
  select(
    -all_of(c(HypAbn, Drop))
  ) %>%
  # Consolidates all observed hamlet species into single SpCode. 2 NEW SPECIES
  # identified ~2013 when surveys began, but updated field guide not obtained
  # for several years. Chances were that many "Florida Hamlet" (one of 2 new
  # species) were misidentified as "Barred Hamlet" during that time. 
  mutate(
    HYSP = rowSums(select(., c(HYSP, HYPU, HYFL))),
    HYPU = NULL,
    HYFL = NULL
  ) %>%
  # Drops taxa with ZERO observations
  select(
    where(
      ~ any(. != 0)  
    )
  ) %>%
  as.data.frame()

rm(dat)

# Import taxonomic details for each species code -----

# Reads in the Excel file as exported from "tbl_Taxa" query in FishSurveys
# Access DB. NOTE: File must be renamed on export to match below, or change
# below to match exported file name

Taxa_table <- 
  data.frame(
    read_excel(
      "Data/Surveys_TaxaTable.xlsx"
    )
  ) %>%
  mutate(
    Genus = 
      replace(
        Genus,
        is.na(Genus),
        "Unknown"
      ),
    Subfamily = 
      replace(
        Subfamily,
        is.na(Subfamily),
        "Unknown"
      ),
    Family =
      replace(
        Family,
          is.na(Family),
          "Unknown"
      ),
    Superfamily = 
      replace(
        Superfamily, 
        is.na(Superfamily), 
        "Unknown"
      ),
    Suborder = 
      replace(
        Suborder,
        is.na(Suborder), 
        "Unknown"
      ),
    Order =
      replace(
        Order,
        is.na(Order), 
        "Unknown"
      ),
    Class = 
      replace(
        Class,
        is.na(Class),
        "Unknown"
      )
  ) %>%
  add_column(
    GSpecies = "NA", 
    .before  = "SpName"
  )

Taxa_table <- 
  mutate(
    Taxa_table,
    GSpecies = 
      ifelse(Taxa_table$Genus != "Unknown",
             paste(paste(substr(Taxa_table$Genus,
                                start = 1, 
                                stop  = 1),
                         ".",
                         sep = ""),
                   sub("^\\S+\\s+", '', 
                       Taxa_table$SpName), 
                   sep = " "),
             paste(Taxa_table$Common))
  ) %>%
  mutate(
    Lateral_Profile     = factor(Lateral_Profile),
    Cross_Section       = factor(Cross_Section),
    Mouth_Position      = factor(Mouth_Position),
    Distribution_Center = factor(Distribution_Center),
    Substrate_Type      = factor(Substrate_Type),
    Feeding_Habit       = factor(Feeding_Habit),
    Water_Column        = factor(Water_Column),
    Gregariousness      = factor(Gregariousness, ordered = TRUE))

# Additional function(s) -----------

cap_pctvars <- function(cap_out){
  # A function that takes the output of the CAPdiscrim function and 
  # produces the percent of among-group variability explained by each
  # canonical axis.
  #
  # Inputs: 
  # cap_out = 'CAPDiscrim' function list output
  # 
  # Output: 
  # pct_var = A numerical array with dimensions 1 x p, where p is the
  #           number of canonical axes, corresponding to the percent
  #           of among-group variability explained by each axis.
  #         
  
  # Extract eigenvalues from manova sublist
  eig     <- cap_out[["manova"]][["Eigenvalues"]]
  
  # Convert eigenvalues to analog of canonical variability
  vars    <- eig*cap_out$tot*(cap_out$varm/100)/(eig+1)
  
  # Calculate percent of variability explained by each axis
  pct_var <- 100*vars/sum(vars)
  
  # Trim columns not corresponding to canonical axes
  pct_var <- pct_var[1:ncol(cap_out$x)]
}

# clears console (toggle on/off as necessary, on initial import good to keep off
# for error tracking)
# cat("\014")