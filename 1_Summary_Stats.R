"""
Summary stats for malaria death dataset

"""
library(dplyr)

# Load the bird deaths database
malaria_death_data <- read.delim("Avian_Malaria_Deaths_For_Analysis.txt", stringsAsFactors = F)
malaria_death_data <- filter(malaria_death_data, Dead_Host_Order != "")

# Load the database of global parasite occurrences
data_final <- read.delim("Global_parasite_database.txt")
data_final$species <- as.character(data_final$species)

###
### Summary of the deaths database
###

# number of studies
length(unique(malaria_death_data$Paper_Title))

# number of continents
length(unique(malaria_death_data$Continent))

# number of countries
length(unique(malaria_death_data$Country))

# number of records
nrow(malaria_death_data)

# number of unique haemosporidian lineages
length(unique(malaria_death_data$Parasite_Cytb_Lineage))

# number in the database for each of the haemosporidian lineages
sort(table(malaria_death_data$Parasite_Cytb_Lineage))

# number of records of each haemosporidian genus
table(malaria_death_data$Parasite_Genus)

# number of hosts that were naive or not to the parasite associated with their death
table(malaria_death_data$Naive_To_Parasite)

# number of hosts in the database that were in captivity or were wild
table(malaria_death_data$Host_Captive_Or_Wild)

# haemosporidian lineages that were not kept for analysis because there is only one record
not_kept <- malaria_death_data %>% filter(Parasite_Distribution == "Unknown_Single_Record") %>% distinct(Parasite_Cytb_Lineage)
not_kept_lins <- not_kept$Parasite_Cytb_Lineage


###
### Summary of the global infection records database
###


# total records of the focal lineages that were associated with host death
malaria_death_lineages <- unique(malaria_death_data$Parasite_Cytb_Lineage)
all_data_filtered_to_death_lineages <- filter(data_final, Lineage_Name %in% malaria_death_lineages)
all_data_filtered_to_death_lineages$Lineage_Name <- as.character(all_data_filtered_to_death_lineages$Lineage_Name)
sum(table(all_data_filtered_to_death_lineages$Lineage_Name)) # this is the number of total global records of the haemosporidians that were associated with the death of a host

# number of records of each focal lineage that was associated with host death
table(all_data_filtered_to_death_lineages$Lineage_Name)


# calculate the number of hosts that were tested in total to generate the global database

# just hosts tested from MalAvi
malavi <- read.delim("MalAvi_10Aug2021.txt", stringsAsFactors = F)
malavi <- filter(malavi, Host_Environment != "Captivity")
malavi <- filter(malavi, !grepl("Olias", Reference_Name))
malavi <- filter(malavi, !grepl("Donovan", Reference_Name))
malavi <- filter(malavi, !grepl("Martinsen et al 2017", Reference_Name))
malavi <- filter(malavi, !grepl("Schrenzel", Reference_Name))
malavi <- filter(malavi, !grepl("Bueno et al 2010", Reference_Name))
malavi <- filter(malavi, !grepl("Ferrell", Reference_Name))
malavi <- filter(malavi, !grepl("Vanstreels et al 2015", Reference_Name))
malavi <- filter(malavi, !grepl("Hill et al 2010", Reference_Name))
malavi <- filter(malavi, !grepl("Silveira et al 2013", Reference_Name))
malavi <- filter(malavi, !grepl("Jia et al 2018", Reference_Name))
malavi <- filter(malavi, !grepl("Verwey et al 2018", Reference_Name))
malavi <- filter(malavi, !grepl("Chagas et al 2013", Reference_Name))
malavi <- filter(malavi, !grepl("Groff et al 2019", Reference_Name))
malavi <- filter(malavi, !grepl("Argilla et al 2013", Reference_Name))
malavi <- filter(malavi, !grepl("Howe et al 2012", Reference_Name))
malavi <- filter(malavi, !grepl("Cannell et al 2013", Reference_Name))
malavi <- filter(malavi, !grepl("Vanstreels et al 2019", Reference_Name))
malavi <- filter(malavi, !grepl("Vanstreels et al 2020", Reference_Name))
malavi <- filter(malavi, !grepl("Meisteret al unpubl", Reference_Name))
malavi$tested[is.na(malavi$tested)] <- 1
malavi_count <- malavi %>% group_by(Reference_Name, continent, order, species, site, tested) %>% summarize() %>% as.data.frame()
sum(malavi_count$tested) # this is the minimum total tested from malavi

# hosts from Galen unpublished data from New York
galen <- read.delim("Galen_NY_malaria_data.txt", stringsAsFactors = F)
galen_ny <- filter(galen, Location == "New York")
galen_edited <- data.frame(Lineage_Name=c(galen_ny$Lineage.1.Name, galen_ny$Lineage.2.Name, galen_ny$Lineage.3.Name, galen_ny$Lineage.4.Name), parasiteGenus=c(galen_ny$Lineage.1.Genus, galen_ny$Lineage.2.Genus, galen_ny$Lineage.3.Genus, galen_ny$Lineage.4.Genus), order=c(galen_ny$BirdTree.Order, galen_ny$BirdTree.Order, galen_ny$BirdTree.Order, galen_ny$BirdTree.Order), family=c(galen_ny$BirdTree.Family, galen_ny$BirdTree.Family, galen_ny$BirdTree.Family, galen_ny$BirdTree.Family), species=c(galen_ny$Host.Latin.Name, galen_ny$Host.Latin.Name, galen_ny$Host.Latin.Name, galen_ny$Host.Latin.Name), continent=rep("North America", nrow(galen_ny)*4), stringsAsFactors = F)
galen_edited_remove_blanks <- filter(galen_edited, Lineage_Name != " " & Lineage_Name != "")
galen_edited_remove_blanks$parasiteGenus[galen_edited_remove_blanks$parasiteGenus == "PL"] <- "Plasmodium"
galen_edited_remove_blanks$parasiteGenus[galen_edited_remove_blanks$parasiteGenus == "PA"] <- "Parahaemoproteus"
galen_edited_remove_blanks$parasiteGenus[galen_edited_remove_blanks$parasiteGenus == "LE"] <- "Leucocytozoon"
galen_edited_remove_blanks$parasiteGenus[galen_edited_remove_blanks$parasiteGenus == "HA"] <- "Haemoproteus"

# hosts tested from Pennsylvania
rushton <- read.delim("PA_malaria_data.txt", stringsAsFactors = F)
rushton_edited <- data.frame(Lineage_Name=c(rushton$Lineage.Name.1, rushton$Lineage.Name.2, rushton$Lineage.Name.3, rushton$Lineage.Name.4), parasiteGenus=c(rushton$Parasite.1, rushton$Parasite.2, rushton$Parasite.3, rushton$Parasite.4), order=c(rushton$Taxon_order, rushton$Taxon_order, rushton$Taxon_order, rushton$Taxon_order), family=c(rushton$Taxon_family, rushton$Taxon_family, rushton$Taxon_family, rushton$Taxon_family), species=c(rushton$Taxon, rushton$Taxon, rushton$Taxon, rushton$Taxon), continent=rep("North America", nrow(rushton)*4), stringsAsFactors = F)
rushton_edited_remove_blanks <- filter(rushton_edited, Lineage_Name != " " & Lineage_Name != "")
rushton_edited_remove_blanks$parasiteGenus[rushton_edited_remove_blanks$parasiteGenus == "Haemoproteus"] <- "Parahaemoproteus"

# total number of hosts that were tested to produce the global database:
sum(malavi_count$tested) + nrow(galen_ny) + nrow(rushton) 


###
### Summary of the SES-MPD analysis
###

# total records in database that were from lineages with negative SESMPD
neg_lins <- c("BUVIR06", "COCOR09", "STOCC16", "TUMER01", "MELSTR01", "PYERY01", "SIAMEX01", "STVAR01", "TASCH01", "TUPHI01", "TURDUS2", "AFTRU5", "CATUST05", "DENVID02", "GASAN01", "GRW04", "GRW11", "LINN1", "PADOM11", "SEIAUR01", "SPMAG04", "SPMAG06", "SYAT05", "TURALB01", "WW3")
neg_records <- filter(malaria_death_data, Parasite_Cytb_Lineage %in% neg_lins) %>% nrow()

pos_lins <- c("EUMIN01", "DENPET03", "FANTAIL01", "GRW06", "LAIRI01", "PADOM09", "PHPAT01", "SGS1", "SPMAG08")
pos_records <- filter(malaria_death_data, Parasite_Cytb_Lineage %in% pos_lins) %>% nrow()

# proportion of total records that were associated with haemosporidians that have negative SESMPD values
neg_records/(neg_records+pos_records)

# MPD summary stats
my_mpd_vals <- data.frame(lins = malaria_death_data$Parasite_Cytb_Lineage, genus = malaria_death_data$Parasite_Genus, MPD = malaria_death_data$SESMPD) %>% distinct()

# Average MPD values per genus
my_mpd_vals %>% group_by(genus) %>% summarize(mean(MPD, na.rm=T))

# Kruskal-Wallis test for difference in MPD values between Parahaem and Plas
my_mpd_vals_reduced <- filter(my_mpd_vals, genus == "Haemoproteus" | genus == "Plasmodium" | genus == "Leucocytozoon")

my_test <- kruskal.test(c(filter(my_mpd_vals_reduced, genus=="Haemoproteus")$MPD, filter(my_mpd_vals_reduced, genus=="Plasmodium")$MPD, filter(my_mpd_vals_reduced, genus=="Leucocytozoon")$MPD), c(rep("Haemoproteus", 10), rep("Plasmodium", 34), rep("Leucocytozoon", 5)))
             