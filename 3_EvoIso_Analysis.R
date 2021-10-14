"""
Notes:
MalAvi database and FASTA downloaded Aug 2021
"""
# load packages for analysis
library(dplyr)
library(ggplot2)
library(picante)
library(ape)

#### prepare final dataset for analysis
data_final <- read.delim("Global_parasite_database.txt")
data_final$species <- as.character(data_final$species)

# load main deaths data file 
malaria_death_data <- read.delim("Avian_Malaria_Deaths_For_Analysis.txt", stringsAsFactors = F)
malaria_death_data <- filter(malaria_death_data, Dead_Host_Order != "")

# create list of malaria lineages that were associated with the death of a host
malaria_death_lineages <- unique(malaria_death_data$Parasite_Cytb_Lineage)
all_data_filtered_to_death_lineages <- filter(data_final, Lineage_Name %in% malaria_death_lineages)
all_data_filtered_to_death_lineages$Lineage_Name <- as.character(all_data_filtered_to_death_lineages$Lineage_Name)
table(all_data_filtered_to_death_lineages$Lineage_Name) # check that everything looks correct

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
# EvoIso analyses
# Get evolutionary isolation estimates for each host that died relative to the typical
# hosts for that parasite
##########

# get list of unique parasite-dead host combinations
paraXdeadHost <- data.frame(Lineage = malaria_death_data$Parasite_Cytb_Lineage, Host = malaria_death_data$BirdTree_Dead_Host_Species, captive_or_wild = malaria_death_data$Host_Captive_Or_Wild, native_or_not = malaria_death_data$Host_Native_Or_Not, naive = malaria_death_data$Naive_To_Parasite, COD = malaria_death_data$Authors_State_Caused_Or_Associated_W_Death, SESMPD = malaria_death_data$SESMPD)
paraXdeadHost <- distinct(paraXdeadHost) # just keep unique rows
paraXdeadHost <- filter(paraXdeadHost, Host != "") # remove blank rows that may have remained
paraXdeadHost <- filter(paraXdeadHost, Lineage != "NIDAV01", Lineage != "PADOM10", Lineage != "COREG01", Lineage != "ANTPAR01", Lineage != "SPMAG07", Lineage != "SPMAG01", Lineage != "SPMAG02", Lineage != "SPMAG09", Lineage != "SPMAG11", Lineage != "EUMIN02", Lineage != "SPMAG10", Lineage != "TUMER12", Lineage != "TUPHI09", Lineage != "TUMER16", Lineage != "BUBSCA02") # removed lineages that are not known from any other hosts other than the one that died
paraXdeadHost <- filter(paraXdeadHost, Lineage != "MEGANT01") # remove because it only infects one host species
paraXdeadHost <- filter(paraXdeadHost, Host != "Cyanoramphus_sp.") # remove because no species level info is available
paraXdeadHost$Lineage <- as.character(paraXdeadHost$Lineage)
paraXdeadHost$Host <- as.character(paraXdeadHost$Host)

# list of unique parasite-normal host associations
paraXnormalHost <- data.frame(Lineage = all_data_filtered_to_death_lineages$Lineage_Name, Host = all_data_filtered_to_death_lineages$species, captive_or_wild = rep("Wild", nrow(all_data_filtered_to_death_lineages)), native_or_not = rep("Native", nrow(all_data_filtered_to_death_lineages)), naive = rep("Not_Naive", nrow(all_data_filtered_to_death_lineages)),  COD = rep("Unknown", nrow(all_data_filtered_to_death_lineages)), SESMPD = rep("Unknown", nrow(all_data_filtered_to_death_lineages)))
paraXnormalHost <- distinct(paraXnormalHost) # just keep unique rows
paraXnormalHost <- filter(paraXnormalHost, Host != "") # remove blank rows that may have remained
paraXnormalHost$Lineage <- as.character(paraXnormalHost$Lineage)
paraXnormalHost$Host <- as.character(paraXnormalHost$Host)
paraXnormalHost <- filter(paraXnormalHost, Host != "Passer_domesticus_x_hispaniolensis")

# make function to get EvoIso metric (weighted by abundance)
# my_df is a dataframe containing separate columns of the hosts that died (the focal host) and the parasite associated with its death (parasite first column, host second column)
# my_mat is a host by parasite matrix with parasites as rows and hosts as columns
# my_host_tree is host phylogeny
get_EvoIso <- function(my_df, my_mat, my_host_tree) {
  for(i in 1:nrow(my_df)){
    parasite <- my_df[i,1] # parasite has to be first column
    dead_host <- my_df[i,2] # focal host has to be second column
    para_all_hosts <- colnames(my_mat)[which(my_mat[parasite,]>=1)] # get the names of the hosts that parasite has been found to infect from MalAvi and other databases
    
    # save the number of infections for each host species from mal_x_host_table
    # e.g. even if a host was infected by the lineage multiple times, each record is listed in infected_hosts
    infected_hosts <- c()
    for (j in para_all_hosts) {
      inf_num <- my_mat[parasite, j]
      infected_hosts <- c(infected_hosts, rep(j, inf_num))
    }
    
    # combine those hosts with the dead host from the database
    para_all_hosts <- unique(c(dead_host, para_all_hosts)) 
    
    if(length(para_all_hosts)==1) {
      my_df$EvoIso[i] <- 0
    } else {
      host_tree <- drop.tip(my_host_tree, my_host_tree$tip.label[-match(para_all_hosts, my_host_tree$tip.label)])
      host_tree_dist_mat <- cophenetic(host_tree)
      host_tree_dist_mat <- host_tree_dist_mat[rownames(host_tree_dist_mat) == dead_host,]
      
      # save the distance between each infected species and the dead host, for EVERY host individual that was infected in the database
      distances <- c()
      for (k in infected_hosts) { 
        my_dist <- host_tree_dist_mat[k]
        distances <- c(distances, my_dist)
      }
      my_df$EvoIso[i] <- mean(distances) # get mean value of all distances
    }
  }
  return(my_df)
}

# run getEvoIso separately on dead and normal hosts
paraXdeadHost_EvoIso_added <- get_EvoIso(paraXdeadHost, mal_x_host_table, all_host_tree)
paraXnormalHost_EvoIso_added <- get_EvoIso(paraXnormalHost, mal_x_host_table, all_host_tree)

# combine into one dataframe
all_evoIso <- rbind(paraXdeadHost_EvoIso_added, paraXnormalHost_EvoIso_added)
all_evoIso$dead_or_not <- c(rep("yes", nrow(paraXdeadHost)), rep("no", nrow(paraXnormalHost)))
all_evoIso$dead_or_not <- as.factor(all_evoIso$dead_or_not)

#######
#
# phylogenetic ANOVA
#
#######

# make new factor column, with two categories:
# not_dead, dead
# full dataset
all_evoIso_for_ANOVA <- all_evoIso
for_ANOVA <- case_when(all_evoIso_for_ANOVA$dead_or_not == "yes" ~ "dead",
                         all_evoIso_for_ANOVA$dead_or_not == "no" ~ "not_dead")
for_ANOVA <- as.factor(for_ANOVA)
names(for_ANOVA) <- all_evoIso_for_ANOVA$Host

# check that all hosts are in the phylogeny
setdiff(all_evoIso_for_ANOVA$Host, all_host_tree$tip.label)

# name a vector of the EvoIso values, and make a new tree that only has the host species in the dataframe
my_evoiso <- all_evoIso_for_ANOVA$EvoIso
names(my_evoiso) <- all_evoIso_for_ANOVA$Host
all_host_tree_ANOVA <- drop.tip(all_host_tree, setdiff(all_host_tree$tip.label, all_evoIso_for_ANOVA$Host))

## ANOVA
library(phytools)
my_phyANOVA <- phylANOVA(all_host_tree_ANOVA, for_ANOVA, my_evoiso)
my_phyANOVA

# plot
ggplot(all_evoIso_for_ANOVA, aes(x=for_ANOVA, y=EvoIso, fill=for_ANOVA)) + geom_boxplot() + theme_classic() + scale_fill_manual(values=c("darkgreen", "gold")) 

##
# test of relationship between host specificity (SESMPD) and naive/not naive
###

# average SES values of parasites that infected naive and not naive hosts
mean(filter(paraXdeadHost_EvoIso_added, naive == "Naive")$SESMPD)
mean(filter(paraXdeadHost_EvoIso_added, naive == "Not_Naive")$SESMPD)

t.test(filter(paraXdeadHost_EvoIso_added, naive == "Naive")$SESMPD, filter(paraXdeadHost_EvoIso_added, naive == "Not_Naive")$SESMPD) # where y and x are numeric

ggplot(paraXdeadHost_EvoIso_added, aes(x=naive, y=SESMPD, fill=naive)) + geom_boxplot() + theme_classic() + scale_fill_manual(values=c("darkgreen", "yellow"))


###
# test of the relationship between host specificity (SESMPD) and evolutionary isolation of the hosts the parasites infected
###
library(caper)

# parasite tree
parasite_tree <- read.nexus("MalAvi_Death_Lineages.tre")
parasite_tree$tip.label <- c("PADOM09", "BUBSCA02", "EUMIN01", "EUMIN02", "MELSTR01", "PYERY01", "SIAMEX01", "STVAR01", "TASCH01", "TUPHI01", "TURDUS2", "BUVIR06", "COCOR09", "MEGANT01", "STOCC16", "TUMER01", "AFTRU5", "ANTPAR01", "CATUST05", "COREG01", "DENPET03", "DENVID02", "FANTAIL01", "GASAN01", "GRW04", "GRW06", "GRW11", "LAIRI01", "LINN1", "NIDAV01", "PADOM10", "PADOM11", "PHPAT01", "SEIAUR01", "SGS1", "SPMAG01", "SPMAG02", "SPMAG04", "SPMAG06", "SPMAG07", "SPMAG08", "SPMAG09", "SPMAG10", "SPMAG11", "SYAT05", "TUMER12", "TUMER16", "TUPHI09", "TURALB01", "WW3")
parasite_tree <- drop.tip(parasite_tree, c("ANTPAR01", "BUBSCA02", "COREG01",  "EUMIN02",  "MEGANT01", "NIDAV01",  "PADOM10",  "SPMAG01","SPMAG02",  "SPMAG07", "SPMAG09",  "SPMAG10",  "SPMAG11",  "TUMER12",  "TUMER16", "TUPHI09"))

simplified_data <- paraXdeadHost_EvoIso_added %>% group_by(Lineage) %>% summarize(evo_mean = mean(EvoIso), ses_mean = mean(SESMPD)) %>% as.data.frame()
comp.data <- comparative.data(parasite_tree, simplified_data, names.col = Lineage, vcv.dim=2)
model2 <- pgls(evo_mean~ses_mean, data=comp.data, lambda = "ML")
summary(model2)

# plot
ggplot(simplified_data, aes(x=ses_mean, y=evo_mean)) + geom_point() + xlab("SES-MPD of Parasite") + ylab("Weighted Evolutionary Isolation of Host That Died") + theme_classic()  + 
  geom_smooth(method='lm')




