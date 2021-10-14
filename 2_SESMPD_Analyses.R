
# Notes: MalAvi database and FASTA downloaded 10 Aug 2021

# load packages for analysis
library(dplyr)
library(picante)
library(ape)


#### prepare final dataset for analysis
data_final <- read.delim("Global_parasite_database.txt")
data_final$species <- as.character(data_final$species)

# load deaths data file 
malaria_death_data <- read.delim("Avian_Malaria_Deaths_For_Analysis.txt", stringsAsFactors = F)
malaria_death_data <- filter(malaria_death_data, Dead_Host_Order != "")
malaria_death_lineages <- unique(malaria_death_data$Parasite_Cytb_Lineage)

# load tree distribution of all host species in malavi database
# Ericson all species
all_host_tree <- read.nexus("output_MCC.tre")

# make sure all the species are in the tree
tree_list <- gsub("_", " ", all_host_tree$tip.label)
all_mal_list <- gsub("_", " ", data_final$species)

setdiff(all_mal_list, tree_list) 

# create parasiteXhost table for phylogenetic diversity analyses
mal_x_host_table <- table(data_final$Lineage_Name, data_final$species)

##########
# SES.MPD analyses
##########

# get SES.MPD values for specific regions based on where the haemosporidian lineages have been found

# Asia, Europe (Lineages: COCOR09)
europe_asia <- filter(data_final, continent == "Asia" | continent == "Europe")
europe_asia$Lineage_Name <- as.character(europe_asia$Lineage_Name)
europe_asia_para_host_table <- table(europe_asia$Lineage_Name, europe_asia$species)
europe_asia_para_host_table <- europe_asia_para_host_table[row.names(europe_asia_para_host_table) %in% malaria_death_lineages,]
ses.mpd_europe_asia_all_hosts <- ses.mpd(europe_asia_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Asia, Europe, Oceania (Lineages: FANTAIL01)
europe_asia_oc <- filter(data_final, continent == "Asia" | continent == "Europe" | continent == "Oceania")
europe_asia_oc$Lineage_Name <- as.character(europe_asia_oc$Lineage_Name)
europe_asia_oc_para_host_table <- table(europe_asia_oc$Lineage_Name, europe_asia_oc$species)
europe_asia_oc_para_host_table <- europe_asia_oc_para_host_table[row.names(europe_asia_oc_para_host_table) %in% malaria_death_lineages,]
ses.mpd_europe_asia_oc_all_hosts <- ses.mpd(europe_asia_oc_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# North America, South America (Lineages: CATUST05, PADOM11, LAIRI01, DENPET03, PADOM09, PHPAT01)
americas <- filter(data_final, continent == "North America" | continent == "South America")
americas$Lineage_Name <- as.character(americas$Lineage_Name)
americas_para_host_table <- table(americas$Lineage_Name, americas$species)
americas_para_host_table <- americas_para_host_table[row.names(americas_para_host_table) %in% malaria_death_lineages,]
ses.mpd_americas <- ses.mpd(americas_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# South America (Lineages: TASCH01, SPMAG06, SPMAG04, TURALB01, SPMAG08)
south_america <- filter(data_final, continent == "South America")
south_america$Lineage_Name <- as.character(south_america$Lineage_Name)
sa_para_host_table <- table(south_america$Lineage_Name, south_america$species)
sa_para_host_table <- sa_para_host_table[row.names(sa_para_host_table) %in% malaria_death_lineages,]
ses.mpd_sa <- ses.mpd(sa_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Oceania (Lineages: MEGANT01)
# Note that SES MPD cannot be calculated for this lineage because it has only been found in one host species
oceania <- filter(data_final, continent == "Oceania")
oceania$Lineage_Name <- as.character(oceania$Lineage_Name)
oceania_para_host_table <- table(oceania$Lineage_Name, oceania$species)
oceania_para_host_table <- oceania_para_host_table[row.names(oceania_para_host_table) %in% malaria_death_lineages,]
ses.mpd_oceania <- ses.mpd(oceania_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# North America (Lineages: SIAMEX01, GASAN01, STOCC16, BUVIR06, MELSTR01, SEIAUR01)
north_america <- filter(data_final, continent == "North America")
north_america$Lineage_Name <- as.character(north_america$Lineage_Name)
na_para_host_table <- table(north_america$Lineage_Name, north_america$species)
na_para_host_table <- na_para_host_table[row.names(na_para_host_table) %in% malaria_death_lineages,]
ses.mpd_na <- ses.mpd(na_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Africa, North America (Lineages: STVAR01)
north_america_afr <- filter(data_final, continent == "North America" | continent == "Africa")
north_america_afr$Lineage_Name <- as.character(north_america_afr$Lineage_Name)
na_afr_para_host_table <- table(north_america_afr$Lineage_Name, north_america_afr$species)
na_afr_para_host_table <- na_afr_para_host_table[row.names(na_afr_para_host_table) %in% malaria_death_lineages,]
ses.mpd_na_afr <- ses.mpd(na_afr_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Africa, Asia, Europe (Lineages: GRW11)
eur_asia_afr <- filter(data_final, continent == "Europe" | continent == "Africa" | continent == "Asia")
eur_asia_afr$Lineage_Name <- as.character(eur_asia_afr$Lineage_Name)
eur_asia_afr_host_table <- table(eur_asia_afr$Lineage_Name, eur_asia_afr$species)
eur_asia_afr_host_table <- eur_asia_afr_host_table[row.names(eur_asia_afr_host_table) %in% malaria_death_lineages,]
ses.mpd_eur_asia_afr <- ses.mpd(eur_asia_afr_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Africa, Asia, Europe, North America, South America, (Lineages: TURDUS2, PYERY01)
eur_na_sa_asia_afr <- filter(data_final, continent == "Europe" | continent == "North America" | continent == "South America" | continent == "Asia" | continent == "Africa")
eur_na_sa_asia_afr$Lineage_Name <- as.character(eur_na_sa_asia_afr$Lineage_Name)
eur_na_sa_asia_afr_host_table <- table(eur_na_sa_asia_afr$Lineage_Name, eur_na_sa_asia_afr$species)
eur_na_sa_asia_afr_host_table <- eur_na_sa_asia_afr_host_table[row.names(eur_na_sa_asia_afr_host_table) %in% malaria_death_lineages,]
ses.mpd_eur_na_sa_asia_afr <- ses.mpd(eur_na_sa_asia_afr_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Africa, Europe, North America, South America (Lineages: TUPHI01, WW3)
eur_na_sa_afr <- filter(data_final, continent == "Europe" | continent == "North America" | continent == "South America" | continent == "Africa")
eur_na_sa_afr$Lineage_Name <- as.character(eur_na_sa_afr$Lineage_Name)
eur_na_sa_afr_host_table <- table(eur_na_sa_afr$Lineage_Name, eur_na_sa_afr$species)
eur_na_sa_afr_host_table <- eur_na_sa_afr_host_table[row.names(eur_na_sa_afr_host_table) %in% malaria_death_lineages,]
ses.mpd_eur_na_sa_afr <- ses.mpd(eur_na_sa_afr_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")
 
# Africa, Asia, Europe, North America, Oceania (Lineages: LINN1)
eur_na_oceania_asia_afr <- filter(data_final, continent == "Europe" | continent == "North America" | continent == "Oceania" | continent == "Asia" | continent == "Africa")
eur_na_oceania_asia_afr$Lineage_Name <- as.character(eur_na_oceania_asia_afr$Lineage_Name)
eur_na_oceania_asia_afr_host_table <- table(eur_na_oceania_asia_afr$Lineage_Name, eur_na_oceania_asia_afr$species)
eur_na_oceania_asia_afr_host_table <- eur_na_oceania_asia_afr_host_table[row.names(eur_na_oceania_asia_afr_host_table) %in% malaria_death_lineages,]
ses.mpd_eur_na_oceania_asia_afr <- ses.mpd(eur_na_oceania_asia_afr_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Africa, Asia, Europe, North America, South America, Oceania (Lineages: GRW04, GRW06, SGS1, SYAT05)
eur_na_sa_asia_afr_oceania <- filter(data_final, continent == "Europe" | continent == "North America" | continent == "South America" | continent == "Asia" | continent == "Africa" | continent == "Oceania")
eur_na_sa_asia_afr_oceania$Lineage_Name <- as.character(eur_na_sa_asia_afr_oceania$Lineage_Name)
eur_na_sa_asia_afr_oceania_host_table <- table(eur_na_sa_asia_afr_oceania$Lineage_Name, eur_na_sa_asia_afr_oceania$species)
eur_na_sa_asia_afr_oceania_host_table <- eur_na_sa_asia_afr_oceania_host_table[row.names(eur_na_sa_asia_afr_oceania_host_table) %in% malaria_death_lineages,]
ses.mpd_eur_na_sa_asia_afr_oceania <- ses.mpd(eur_na_sa_asia_afr_oceania_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Asia, South America (Lineages: DENVID02)
asia_sa <- filter(data_final, continent == "Asia" | continent == "South America")
asia_sa$Lineage_Name <- as.character(asia_sa$Lineage_Name)
asia_sa_para_host_table <- table(asia_sa$Lineage_Name, asia_sa$species)
asia_sa_para_host_table <- asia_sa_para_host_table[row.names(asia_sa_para_host_table) %in% malaria_death_lineages,]
ses.mpd_asia_sa <- ses.mpd(asia_sa_para_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Africa, Europe (Lineages: TUMER01)
af_eur <- filter(data_final, continent == "Africa" | continent == "Europe")
af_eur$Lineage_Name <- as.character(af_eur$Lineage_Name)
af_eur_host_table <- table(af_eur$Lineage_Name, af_eur$species)
af_eur_host_table <- af_eur_host_table[row.names(af_eur_host_table) %in% malaria_death_lineages,]
ses.mpd_af_eur <- ses.mpd(af_eur_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Oceania, South America (Lineages: EUMIN01)
oceania_SA <- filter(data_final, continent == "Oceania" | continent == "South America")
oceania_SA$Lineage_Name <- as.character(oceania_SA$Lineage_Name)
oceania_SA_host_table <- table(oceania_SA$Lineage_Name, oceania_SA$species)
oceania_SA_host_table <- oceania_SA_host_table[row.names(oceania_SA_host_table) %in% malaria_death_lineages,]
ses.mpd_oceania_SA_all_hosts <- ses.mpd(oceania_SA_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")

# Africa, Asia, Europe, North America (Lineages: AFTRU5)
asia_af_NA_eur <- filter(data_final, continent == "Europe" | continent == "North America" | continent == "Asia" | continent == "Africa")
asia_af_NA_eur$Lineage_Name <- as.character(asia_af_NA_eur$Lineage_Name)
asia_af_NA_eur_host_table <- table(asia_af_NA_eur$Lineage_Name, asia_af_NA_eur$species)
asia_af_NA_eur_host_table <- asia_af_NA_eur_host_table[row.names(asia_af_NA_eur_host_table) %in% malaria_death_lineages,]
ses.mpd_asia_af_NA_eur <- ses.mpd(asia_af_NA_eur_host_table, cophenetic(all_host_tree), abundance.weighted = TRUE, null.model = "independentswap")
