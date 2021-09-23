library(readr)
library(ggplot2)
library(Rmisc)
library(dplyr)
library(biomartr)
library(reutils)
library(tidyverse)
library(stringr)

db <- read_delim("CRISPR_crisprdb.csv", 
                              "\t", escape_double = FALSE, trim_ws = TRUE)

CC_loci <- read_delim("CRISPRcity_loci.tsv", 
                              "\t", escape_double = FALSE, col_names = c('x', 'n', 'name', 'seq_id', 'type'), 
                              trim_ws = TRUE)[,3:5]

CC_arrays <- read_delim("CRISPRcity_CRISPR_arrays.tsv", 
                                "\t", escape_double = FALSE, trim_ws = TRUE)
CC_arrays$crispr_id <- rownames(CRISPRcity_arrays)


colnames_mapping <- c('SequenceName',	'SequenceRole',	'AssignedMolecule',	'AssignedMoleculeLocationType',	'GenBankAccn',	'Relationship',	'RefSeqAccn',	'AssemblyUnit',	'SequenceLength', 'UCSC_style_name')
genbank_refseq_mapping <- read_delim("genbank_refseq_mapping.txt", 
                                     "\t", escape_double = FALSE, trim_ws = TRUE, col_names = colnames_mapping)

CC_loci_mapped <- left_join(CC_loci, genbank_refseq_mapping.f, by = c('seq_id' = 'gbk_accession'))

# SPACEROM_spacers_1 <- read_delim("spacerom_1.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)
# SPACEROM_spacers_2 <- read_delim("spacerom_2.tsv", "\t", escape_double = FALSE, trim_ws = TRUE) %>% mutate(spacer_start = NA) %>% 
#   select(genbank_id, c_start, c_end, spacer_n, spacer_start, length, sequence)

SPACEROM_spacers <- rbind(SPACEROM_spacers_1, SPACEROM_spacers_2) %>% mutate(gbk_id = gsub("\\..*", "", genbank_id))
SPACEROM_spacers$spacer_id <- rownames(SPACEROM_spacers)

SPACEROM_annotated <- left_join(SPACEROM_spacers, CC_ranges, by = c('gbk_id' = 'genbank_id', 'c_start' = 'start', 'c_end' = 'end'))
write.table(SPACEROM_annotated, 'SPACEROM_mapped_to_CRISPRicity.tab', col.names = T, row.names = F, quote = F)

SPACEROM_annotated %>% nrow()
SPACEROM_annotated %>% filter(!is.na(crispr_id), length < 40, length > 15) %>% write.table('SPACEROM_mapped_to_CRISPRicity.clean.tab', col.names = T, row.names = F, quote = F)

spy_arrays <- CC_arrays %>% dplyr::filter(grepl('Streptococcus_pyogenes', name), str_detect(type, 'CAS-II-'))
SPACEROM_annotated %>% filter(!is.na(crispr_id), length < 40, length > 15, crispr_id %in% spy_arrays$crispr_id) %>% write.table('SPACEROM_mapped_to_CRISPRicity.clean.spy.tab', col.names = T, row.names = F, quote = F)


# 23844 spacers not mapped

# ------------ WE HAVE MAPPED CRISPR arrays to spacers -------------------------------------
# ------ Start selecting arrays by type

arrays = list()
arrays[['CAS-II']] <- CC_arrays %>% filter(str_detect(type, '^CAS-II-'))
arrays[['CAS-I']] <- CC_arrays %>% filter(str_detect(type, '^CAS-I-'))
arrays[['CAS-III']] <- CC_arrays %>% filter(str_detect(type, '^CAS-III-'))
arrays[['CAS-VI']] <- CC_arrays %>% filter(str_detect(type, '^CAS-VI-'))

CRISPRcity_arrays %>% group_by(type) %>% summarize(n = n())

# ---------------- ISLANDS
# 1501564490	40023	40727	+	Streptococcus_pyogenes_Lacen4_GCA_000685545.1	JHTO01000050	11	cd09752	cas5,G5,CASCADE	234	Seed	CLUSTER_5709	CAS-I-C
CC_islands <- read_delim("CRISPRcity_islands.flat.tsv", 
                            "\t", escape_double = FALSE, col_names = colnames(CC_arrays)[1:13], 
                            trim_ws = TRUE)

CC_islands <- left_join(CC_islands, CC_arrays)

CC_islands %>% group_by(name, genbank_id, type) %>% group_indices(name, genbank_id, type)

# --------------- ALL elements ----------------------------------------
gbk_assembly_taxid_mapping <- read_delim("assembly_summary_genbank.txt", 
                                         "\t", escape_double = FALSE, trim_ws = TRUE) %>% select(assembly_accession, taxid, species_taxid)

CX_islands <- read_delim("CRISPRcity_islands.flat.id.tsv", 
                         "\t", escape_double = FALSE, col_names = c(colnames(CC_arrays)[1:13], 'crispr_id'), 
                         trim_ws = TRUE)
CX_islands.arrays <- CX_islands %>% filter(profile == 'CRISPR')
CX_islands.arrays <- CX_islands.arrays %>% dplyr::mutate(assembly_id = str_extract(name, "GCA.+"))
CX_islands.arrays <- left_join(CX_islands.arrays, gbk_assembly_taxid_mapping, by = c('assembly_id' = 'assembly_accession'))
CX_ranges <- CX_islands.arrays %>% select(start, end, genbank_id, crispr_id, taxid, species_taxid)

CX_spacers <- left_join(SPACEROM_spacers, CX_ranges, by = c('gbk_id' = 'genbank_id', 'c_start' = 'start', 'c_end' = 'end'))
CX_spacers$spacer_id <- rownames(CX_spacers)
CX_spacers %>% filter(!is.na(crispr_id), length < 40, length > 15) %>% write.table('CX_spacers.clean.tab', col.names = T, row.names = F, quote = F)

# ----------------------------
CX_spacers <- read_delim("CX/CX_spacers.clean.tab", 
                         "\t", escape_double = FALSE, 
                         trim_ws = TRUE)

CX_spacers.taxid <- left_join(CX_spacers, CX_ranges, by = c('gbk_id' = 'genbank_id', 'c_start' = 'start', 'c_end' = 'end')) %>% select(-one_of('crispr_id.y'))
colnames(CX_spacers.taxid)[10] <- c('crispr_id')

# --------------------------- FIND ARRAYS WITH CAS4 -----
CX_islands.cas4 <- CX_islands %>% group_by(crispr_id) %>% summarise(has_cas4 = 'cd09637' %in% profile) %>% filter(has_cas4 == TRUE)
CX_islands.cas4 <- CX_islands %>% filter(crispr_id %in% CX_islands.cas4$crispr_id)
CX_spacers.cas4 <- CX_spacers.taxid %>% filter(crispr_id %in% CX_islands.cas4$crispr_id)

CX_spacers.cas4 %>% write.table('CX/CX_spacers.cas4.tab', col.names = T, row.names = F, quote = F, sep = '\t')
CX_islands.cas4 %>% write.table('CX/CX_islands.cas4.tab', col.names = T, row.names = F, quote = F, sep = '\t')

CX_islands.arrays %>% write.table('CX/CX_islands.all.tab', col.names = T, row.names = F, quote = F, sep = '\t')
CX_spacers.taxid %>% write.table('CX/CX_spacers.all.tab', col.names = T, row.names = F, quote = F, sep = '\t')


CX_islands %>% filter(grepl('Streptococcus_pyogenes', name)) %>% write.table('test/CX_islands.spy.tab', col.names = T, row.names = F, quote = F, sep = '\t')
