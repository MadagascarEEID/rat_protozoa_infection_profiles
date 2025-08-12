# script for filtering NC data and converting ASVs into OTUs
# Workflow:
# A. pre-filtering
#   1. filtering only for small mammals (all species)
#   2. for duplicated samples (two buffers) - excluding the NAP buffer replicate
#   3. changing samples with total reads < 500  - no nematodes
# B. converting ASVs into OTUs by 97% similarity
#   * the most common ASV of each OTU cluster is marked

# output:
# NC_wide_otu.csv - matrix of individual rats (rows) and protozoa OTUs (cols)
# NC_long_otu.csv - final host-OTUs data (including all infected/uninfected hosts)
# NC_taxonomy.csv - taxonomy of OTUs and the matching ASV


# Load necessary libraries
library(tidyverse)
library(DECIPHER)
library(Biostrings)
library(rstudioapi)
library(ape)

rm(list=ls())


# defining the species of interest
sm_species <- c("Rattus rattus", "Suncus murinus", "Mus musculus", "Microgale brevicaudata", "Mus musculus", "Eliurus", "Eliurus minor", "Nesogale talazaci", "Microgale prolixacaudata", "Nesomys audeberti", "Microgale longicaudata", "Microgale cf. talazaci", "Setifer setosus", "Eliurus webbi", "Suncus etruscus", "Tenrec ecaudatus", "Nesomys rufus", "Microgale parvula", "Eliurus cf. tanala")

# reading small mammals data
rats <- read_csv("data/data_raw/small_mammals/Terrestrial_Mammals.csv") %>% 
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id))) %>% 
  dplyr::select(host_ID, field_identification) %>%
  filter(field_identification == "Rattus rattus")

# organizing raw NC asv data
NC_wide_asv <- read_csv("data/data_raw/NC/ALL_NC_1pct_RRA_NIH_filtered_samples.csv") %>% 
  filter(Species %in% sm_species) %>%
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", Sample_Name))) %>% # readable name
  filter(note!="NAP Replicate" | is.na(note)) %>% # remove replicates
  select(c(host_ID, reads, starts_with("ASV"))) %>% 
  mutate(across(-1, ~ if_else(reads <= 500, 0, .))) %>% # reads under 500 are set to 0, no nematodes
  select(-which(colSums(., na.rm = T) == 0)) # remove columns where sum is 0 (ASVs with no links)

# sanity check - no ASVs with no reads
(colSums(NC_wide_asv[,-1]) <= 0) %>% sum(na.rm = T) 

# vector of ASVs names
nc_asv_to_cluster <- colnames(NC_wide_asv[-c(1:2)])


# Convert ASVs to OTUs by 97% Similarity 

# reading the full dna sequences (fasta file)
seq_fa <- read.FASTA(file="data/data_raw/NC/ASV_NC_all_village_1pct.fa")

# filtering the fasta file for relevant ASVs
seq_fa2 <- seq_fa[names(seq_fa) %in% nc_asv_to_cluster]
write.FASTA(seq_fa2, file = "data/data_raw/NC/NC_filtered.fa") # writing the filtered fasta file

seqs <- readDNAStringSet("data/data_raw/NC/NC_filtered.fa")
seqs <- OrientNucleotides(seqs)

# clustering for 97% similarity
set.seed(123)
clusters <- DECIPHER::Clusterize(
  seqs, 
  method = "overlap",
  penalizeGapLetterMatches = NA,
  cutoff = 0.03) # use `cutoff = 0.03` for a 97% OTU 

# saving the clustering
pathogen_cluster <- clusters %>% 
  rownames_to_column("asv_ID") %>% 
  dplyr::rename("otu_ID"="cluster")

length(unique(pathogen_cluster$otu_ID))

NC_OTUs_97 <- pathogen_cluster


# assigning taxonomy

NCASVs_taxonomy80 <- read_tsv("data/data_raw/NC/NCASVs_taxonomy80.tsv") %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species) %>%
  dplyr::rename(asv_ID = ASV)

# creating a table of the most common asv and its taxa (per OTU) 
NC_taxonomy <- NC_wide_asv %>%
  pivot_longer(starts_with("ASV"), names_to="asv_ID", values_to="relative_reads") %>%
  filter(relative_reads > 0) %>%
  count(asv_ID) %>% 
  left_join(NC_OTUs_97, by="asv_ID") %>%
  group_by(otu_ID) %>%
  mutate(asv_ID_top = asv_ID[which.max(n)]) %>%
  select(otu_ID, asv_ID, asv_ID_top) %>%
  unique() %>%
  left_join(NCASVs_taxonomy80, by=c("asv_ID_top" = "asv_ID")) %>%
  arrange(otu_ID) %>%
  ungroup() %>%
  mutate(otu_ID = as.character(otu_ID))


# sanity check - that the OTUs are not crazy taxonomically
NC_OTUs_97 %>%
  left_join(NCASVs_taxonomy80, by="asv_ID") %>%
  group_by(otu_ID) %>%
  summarise(dist_orders = n_distinct(Order, na.rm = TRUE),
            dist_families = n_distinct(Family, na.rm = TRUE),
            dist_genues = n_distinct(Genus, na.rm = TRUE),
            dist_species = n_distinct(Species, na.rm = TRUE)) %>%
  filter(dist_orders > 1 | dist_families > 1 | dist_genues > 1 | dist_species > 1)


# making final OTU table

NC_long_otu <- NC_wide_asv %>%
  filter(host_ID %in% rats$host_ID) %>%
  pivot_longer(starts_with("ASV"), names_to="asv_ID", values_to="relative_reads") %>%
  left_join(NC_taxonomy, by="asv_ID") %>%
  group_by(host_ID, otu_ID) %>%
  summarise(relative_reads = sum(relative_reads)) %>%
  ungroup() %>%
  left_join(NC_taxonomy %>% select(-asv_ID) %>% unique, by="otu_ID") %>%
  mutate(link = if_else(relative_reads > 0,1,0))

NC_wide_otu <- NC_long_otu %>%
  select(host_ID, otu_ID, link) %>%
  pivot_wider(names_from="otu_ID", values_from="link") %>%
  select(-which(colSums(., na.rm = T) == 0)) # remove columns where sum is 0 (OTUs with no links to rats, but other SMs...)

NC_long_otu <- NC_long_otu %>%
  filter(otu_ID %in% colnames(NC_wide_otu)[-1]) # remove OTUs with no rat links

length(unique(NC_long_otu$host_ID)) # 879 hosts
length(unique(NC_long_otu$otu_ID)) # 94 OTUs

# writing final tables
write_csv(NC_wide_otu, "data/data_processed/NC_wide_otu.csv")
write_csv(NC_long_otu, "data/data_processed/NC_long_otu.csv")
write_csv(NC_taxonomy, "data/data_processed/NC_taxonomy.csv")
