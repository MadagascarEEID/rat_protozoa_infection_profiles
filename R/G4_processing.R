# script for filtering G4 protozoa data and converting ASVs into OTUs
# Workflow:
# A. pre-filtering
#   1. filtering only for small mammals (all species)
#   2. for duplicated samples (two buffers) - excluding the NAP buffer replicate
#   3. removing samples with total reads < 500
# B. converting ASVs into OTUs by 97% similarity
#   * the most common ASV of each OTU cluster is marked
# C. filtering parasitic OTUs
#   1. filtering only OTUs that their central ASV is parasitic 

# output:
# G4_wide_otu.csv - matrix of individual rats (rows) and protozoa OTUs (cols)
# G4_long_otu.csv - final host-OTUs data (including all infected/uninfected hosts)
# G4_taxonomy.csv - taxonomy of OTUs and the matching ASV


# Load necessary libraries
library(tidyverse)
library(DECIPHER)
library(Biostrings)
library(ape)

rm(list=ls())


# defining the species of interest
sm_species <- c("Rattus rattus", "Suncus murinus", "Mus musculus", "Microgale brevicaudata", "Mus musculus", "Eliurus", "Eliurus minor", "Nesogale talazaci", "Microgale prolixacaudata", "Nesomys audeberti", "Microgale longicaudata", "Microgale cf. talazaci", "Setifer setosus", "Eliurus webbi", "Suncus etruscus", "Tenrec ecaudatus", "Nesomys rufus", "Microgale parvula", "Eliurus cf. tanala")

# reading small mammals data
rats <- read_csv("data/data_raw/small_mammals/Terrestrial_Mammals.csv") %>% 
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", animal_id))) %>% 
  dplyr::select(host_ID, field_identification) %>%
  filter(field_identification == "Rattus rattus")


# filtering raw G4 asv data
G4_wide_asv <- read_csv("data/data_raw/G4/ALL_G4_1pct_RRA_NIH_filtered_samples.csv") %>% 
  filter(Species %in% sm_species) %>% # only small mammals
  mutate(host_ID = as.numeric(gsub(".*?([0-9]+).*", "\\1", Sample_Name))) %>% # readable name
  filter(note!="NAP Replicate" | is.na(note)) %>% 
  select(c(host_ID, reads, starts_with("ASV"))) %>%
  filter(reads > 500) %>% # reads <=500
  select(-which(colSums(., na.rm = T) == 0)) # remove columns where sum is 0 (ASVs with no links)

# sanity check - no ASVs with no reads
(colSums(G4_wide_asv[,-1]) <= 0) %>% sum() 

# vector of ASVs names
g4_asv_to_cluster <- colnames(G4_wide_asv[-c(1:2)])



# Convert ASVs to OTUs by 97% Similarity 

# reading the full dna sequences (fasta file)
seq_fa <- read.FASTA(file="data/data_raw/G4/ASV_G4_all_village_1pct_2.fa")

# filtering the fasta file for relevant ASVs
seq_fa2 <- seq_fa[names(seq_fa) %in% g4_asv_to_cluster]
write.FASTA(seq_fa2, file = "data/data_raw/G4/G4_filtered.fa") # writing the filtered fasta file

seqs <- readDNAStringSet("data/data_raw/G4/G4_filtered.fa")
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

G4_OTUs_97 <- pathogen_cluster


# assigning taxonomy

# reading G4 ASV data with parasite information
G4_asv_taxa <- read_csv("data/data_raw/G4/G4_reconciled_parasite_taxa_nonhuman.csv")


# creating a table of the most common asv and its taxa (per OTU) 
# Parasite == "Y" - known potential mammalian parasites 
G4_taxonomy <- G4_wide_asv %>%
  pivot_longer(starts_with("ASV"), names_to="asv_ID", values_to="relative_reads") %>%
  filter(relative_reads > 0) %>%
  count(asv_ID) %>% 
  left_join(G4_OTUs_97, by="asv_ID") %>%
  group_by(otu_ID) %>%
  mutate(asv_ID_top = asv_ID[which.max(n)]) %>%
  select(otu_ID, asv_ID, asv_ID_top) %>%
  unique() %>%
  left_join(G4_asv_taxa, by=c("asv_ID_top" = "ASV")) %>%
  arrange(otu_ID) %>%
  ungroup() %>%
  mutate(otu_ID = as.character(otu_ID)) %>%
  filter(Parasite == "Y")


# use this taxa list to review the literature and decide which are *protozoa* parasites
G4_taxonomy %>%
  select(final_id) %>%
  unique() %>%
  pull()

# only the parasitic protozoa (excluding non-protozoa ASVs)
g4_sp <- c("Hypotrichomonas","Entamoeba","Hexamastix","Tritrichomonas","Simplicimonas_similis","Tetratrichomonas","Pentatrichomonas_hominis","Eimeria","Blastocystis","Cryptosporidium")


G4_long_otu <- G4_wide_asv %>%
  filter(host_ID %in% rats$host_ID) %>% 
  pivot_longer(starts_with("ASV"), names_to="asv_ID", values_to="relative_reads") %>%
  left_join(G4_taxonomy, by="asv_ID") %>%
  filter(final_id %in% g4_sp) %>%
  group_by(host_ID, otu_ID,final_id) %>%
  summarise(relative_reads = sum(relative_reads)) %>%
  ungroup() %>%
  mutate(link = if_else(relative_reads > 0,1,0)) %>%
  select(host_ID, otu_ID, link, final_id)


# final tables
G4_wide_otu <- G4_long_otu %>%
  select(-final_id) %>%
  pivot_wider(names_from="otu_ID", values_from="link") %>%
  select(-which(colSums(., na.rm = T) == 0)) # remove columns where sum is 0 (OTUs with no links to rats, but other SMs...)

G4_long_otu <- G4_long_otu %>%
  filter(otu_ID %in% colnames(G4_wide_otu)[-1]) # remove OTUs with no rat links

length(unique(G4_long_otu$host_ID)) # 841 hosts
length(unique(G4_long_otu$otu_ID)) # 41 OTUs
# connectance
sum(G4_long_otu$link) / nrow(G4_long_otu) # 0.045

# writing final tables
write_csv(G4_wide_otu, "data/data_processed/G4_wide_otu.csv")
write_csv(G4_long_otu, "data/data_processed/G4_long_otu.csv")
write_csv(G4_taxonomy, "data/data_processed/G4_taxonomy.csv")



