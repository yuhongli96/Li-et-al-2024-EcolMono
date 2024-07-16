library(renv)
renv::restore() # restore the r libraries that are used in this script development
library(tidyverse)
library(lme4)
library(sjPlot)
library(here)

# # read the plant_species_list info
# list <- read.csv(file="https://docs.google.com/spreadsheets/d/e/2PACX-1vTBd3rlvUtVHKrbUn80JK4LF-uDMDBD9VBFkYKh29XHzBxdwL2ti7MQjMzS8bmyjMpgwXxSfWnSGZVB/pub?gid=705823615&single=true&output=csv") %>%
#   dplyr::tibble() %>%
#   dplyr::select(SpCode6, Sci_Name) %>%
#   dplyr::rename(spcode6=SpCode6,
#                 sci_name=Sci_Name)
# site_id2 <- data.frame(site_id = c("SC0001", "SC0006", "SC0042", "SC0048",
#                                    "SC0058", "SC0059", "SC0060", "SC0128", "SC0129"), 
#                        site_id2 = c("5", "4", "1", "3", "6", "8", "7", "2", "9"))
# 
# # add environmental data to the all dataset
# envir<-read.csv(file = "https://docs.google.com/spreadsheets/d/e/2PACX-1vTBd3rlvUtVHKrbUn80JK4LF-uDMDBD9VBFkYKh29XHzBxdwL2ti7MQjMzS8bmyjMpgwXxSfWnSGZVB/pub?gid=1921808020&single=true&output=csv") %>%
#   dplyr::tibble() %>%
#   dplyr::select(site_id, rainfall_2000.2021_mm) 
# names(envir)
# 
# nutri <- read.csv("./data/sample_nutrient.csv", header = T) %>%
#   dplyr::left_join(list, by="spcode6") %>%
#   dplyr::left_join(site_id2, by="site_id") %>%
#   dplyr::left_join(envir, by="site_id")
# 
# names(nutri)
# nutri <- nutri[, c(2, 11, 4, 12, 13, 3, 5, 6, 10, 7, 9, 8)]
# nutri <- nutri %>%
#   dplyr::rename(site_id=site_id2,
#                 annual_rainfall_mm=rainfall_2000.2021_mm)
# dat <- nutri %>%
#   dplyr::arrange(., site_id)
#write.csv(dat, "D:/data.csv", row.names=F)
dat <- read_csv("./data/data.csv") %>%
  dplyr::mutate(NP = nitrogen_per/phosphorus_per, # calculate the nutrient ratios
                CaP = calcium_per/phosphorus_per,
                KNa = potassium_per/sodium_per)

#### estimate how site-average elemental contents and ratios varied with rainfall, associated with Figure 3 and Figure 6------------------
# calculate the weighted average elemental contents and ratios for each site
trasct_nutri <- dat %>% 
  dplyr::group_by(site_id) %>%
  dplyr::mutate(layer_total=sum(layer_abundance)) %>%
  dplyr::mutate(occr_proprt=layer_abundance/layer_total) %>% # calculate the relative layer abundance for each species at each site
  dplyr::reframe(nitrogen.sp=occr_proprt*nitrogen_per,
                 phosphorus.sp=occr_proprt*phosphorus_per,
                 potassium.sp=occr_proprt*potassium_per,
                 calcium.sp=occr_proprt*calcium_per,
                 magnesium.sp=occr_proprt*magnesium_per,
                 sodium.sp=occr_proprt*sodium_per) %>%
  dplyr::group_by(site_id) %>%
  dplyr::reframe(nitrogen=sum(nitrogen.sp, na.rm = T), # calculate the weighted average elemental contents for each site
                 phosphorus=sum(phosphorus.sp, na.rm = T),
                 potassium=sum(potassium.sp, na.rm = T),
                 calcium=sum(calcium.sp, na.rm = T),
                 magnesium=sum(magnesium.sp, na.rm = T),
                 sodium=sum(sodium.sp, na.rm = T)) %>% 
  dplyr::left_join(unique(dat[, c("site_id", "annual_rainfall_mm")]), by="site_id") %>%
  dplyr::mutate(ratioNP = nitrogen/phosphorus,
                ratioCaP = calcium/phosphorus,
                ratioKNa = potassium/sodium) 

# linear models of the site average element contents along the rainfall gradient
c1<-c("nitrogen", "phosphorus", "potassium", 
      "log10(calcium)", "log10(magnesium)", "log10(sodium)")
for (i in 1:length(c1)) {
  lm <- lm(paste(c1[i], "~annual_rainfall_mm"), 
           data = trasct_nutri)
  print(c1[i])
  print(summary(lm))
  cat("\n")
}


# linear models of the site element ratios along the rainfall gradient
c1 <- c("ratioNP", "ratioCaP", "ratioKNa")
for (i in 1:length(c1)) {
  lm <- lm(paste(c1[i], "~annual_rainfall_mm"), 
           data = trasct_nutri)
  print(c1[i])
  print(summary(lm))
  cat("\n")
}


#### estimate how intraspecific leaf nutrient contents vary along the rainfall gradient, associated with Figure 4 and Figure 6------------------
# select the plant species sampled at >=3 sites
dat2 <- dat %>%
  dplyr::group_by(spcode6) %>%
  dplyr::reframe(n=n()) %>%
  dplyr::filter(n>=3)
plant_sp <- dat2$spcode6
dat3 <- dat[dat$spcode6 %in% plant_sp,] %>%
  dplyr::mutate(NP = nitrogen_per/phosphorus_per,
                CaP = calcium_per/phosphorus_per,
                KNa = potassium_per/sodium_per)

# general linear mixed model of nutrient content variation within species along the rainfall gradient
c1 <- c("nitrogen_per", "phosphorus_per", "potassium_per",
        "log10(calcium_per)", "log10(magnesium_per)", "log10(sodium_per)")
for (i in 1:length(c1)) {
  glmm <- lme4::lmer(paste(c1[i], "~ annual_rainfall_mm + (1 | spcode6)"), 
                     data = dat3)
  print(c1[i])
  print(summary(glmm))
  print(confint(glmm))
  print(sjPlot::tab_model(glmm))
  cat("\n")
}

# general linear mixed model of nutrient ratio variation within species along the rainfall gradient
c1 <- c("NP", "log10(CaP)", "KNa")
for (i in 1:length(c1)) {
  glmm <- lme4::lmer(paste(c1[i], "~ annual_rainfall_mm + (1 | spcode6)"), 
                     data = dat3)
  print(c1[i])
  print(summary(glmm))
  print(confint(glmm))
  print(sjPlot::tab_model(glmm))
  cat("\n")
}


#### decompose the total element variation across sites into intraspecific variation and species turnover, associated with Figure 5 and Figure S3----------
# reference:
# chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://digital.csic.es/bitstream/10261/221270/3/R_Material_traits.pdf
# refer to 6.4 ITV between communities
# the calculation used an external R project saved in this Git project
source(here::here("./scripts/trait-flex.R"))

# assign the plant species id and site id
spID <- dat$spcode6 
plotID <- dat$site_id

# create a matrix where the header is species name, and rows show the relative abundance of species in each plot
commMatrix <- dat %>%
  dplyr::mutate(layer_total=sum(layer_abundance)) %>%
  dplyr::mutate(occr_proprt=layer_abundance/layer_total) %>% 
  dplyr::select(site_id, spcode6, occr_proprt) %>%
  tidyr::pivot_wider(names_from = spcode6, values_from = occr_proprt, values_fill = 0) %>%
  dplyr::arrange(., site_id) %>%
  tibble::column_to_rownames(., var = "site_id")%>%
  data.matrix(., rownames.force = NA)

# prepare the environmental gradient
gradient <- dat %>%
  dplyr::select(site_id, annual_rainfall_mm) %>%
  dplyr::arrange(., site_id) %>%
  distinct(.) %>%
  dplyr::mutate(annual_rainfall_mm=ifelse(site_id==4, 797.1, annual_rainfall_mm)) %>% # because the model doesn't allow two sites with identical environmental conditions, we modified the rainfall level of site 4 from 797 to 797.1
  tibble::column_to_rownames(., var = "annual_rainfall_mm")
gradient<-as.numeric(rownames(gradient))

# decompose the total variation to species turnover and intraspecific variation
c1 <- c("nitrogen_per", "phosphorus_per", "potassium_per", 
        "calcium_per", "magnesium_per", "sodium_per")

# to run the analysis for nutrient ratios by running the code of next line
# c1 <- c("NP", "CaP", "KNa")

for (i in 1:length(c1)) {
  traitVector <- dat[[c1[i]]]
  traitsForCWM <- trait.transform(spIDs = as.factor(spID),
                                  traitVals = traitVector,
                                  sampleID = plotID)
  
  CWMSpecFix<-trait.CWM(matLocal = traitsForCWM$matLocal,
                        matFixed = traitsForCWM$matFixed,
                        commMatrix = commMatrix)
  CWMSpecFix # show the community-weighted mean calculated by local mean and grand mean
  
  bITV<-trait.flex.anova(formula = ~1, specif.avg = CWMSpecFix$specific.avg,
                         fixed.avg = CWMSpecFix$fixed.avg)
  print(c1[i])
  print(bITV)
  
  # check how much variation of each type is explained by rainfall
  bITV2<-trait.flex.anova(formula = ~gradient,
                          specif.avg = CWMSpecFix$specific.avg,
                          fixed.avg = CWMSpecFix$fixed.avg)
  print(bITV2)
  cat("\n")
}
