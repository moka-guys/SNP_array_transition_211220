library(tidyverse)
library(dplyr)

################################
# Using the tsv output from moka_array_legacy_data_export.py
# This script has been created to:
#  1) Tidy and edit legacy data from Moka 
#  2) Export data for liftover to hg38
#  3) Merge liftover coordinates with patient information 
#  4) Export data to be loaded into UCSC
##################################


##### Load data  ####### ====================================================================

setwd("/home/erin/Documents/Work/SNP_array_liftover/moka_export")
dataframe <- read.delim("moka_array_export_220215.tsv", sep = '\t' , header = TRUE)

### Tidy data ######### ====================================================================

array_df <- dataframe %>% 
  select(-c(X)) %>%  # Drop Python export unique ID column
  group_by(PatientID) %>%
  mutate(phenotype_string = toString(Phenotype)) %>% # Group phenotypes into one row  per patient
  ungroup() %>% # ungroup to be able to remove column
  select(-c(Phenotype)) %>% # remove phenotype column 
  distinct() %>% # Return only one of each unique row 
  # Create a new column for sex, merge all the different sex/gender columns into one, label all the pts with no gender
  mutate(sex = case_when(Gender == "Male" | Sexed == "M" | Sexed == "M\n\n\nM" |  
                           Sexed == "M M" |  Sexed == "m"  |
                           BookinSex == "M" ~ "Male",
                        Gender == "Female" | Sexed == "F" | Sexed == "f" |  Sexed == "F F" |  Sexed == "Female" | BookinSex == "F" ~ "Female",
         TRUE ~ "blank"), #  highlight unknowns/blanks to remove
         # Create a new column grouping tissue types.
         tissue_groups = case_when(Referral == "General referral" | Referral == "Repeat general referral" ~ "General referral", 
                                   Referral == "Prenatal" | Referral == "Repeat prenatal" ~ "Prenatal", 
                                   TRUE ~ "Tissues"), #POC/Tissues, repeat POC/Tissues & pseudorush all go in this group
         # Calculate patients age when test was requested in days
         age_days = as.double(difftime(lubridate::ymd_hms(RequestedDate), 
                                  lubridate::ymd_hms(DoB), 
                                  units = "days")),
         # Only ages for gen refs are relevant, other tissues are the age of the mother
         age_years_days = case_when(tissue_groups == "General referral" & age_days > 364 ~ str_c(as.character(as.integer(age_days/365)), " years"), 
                               tissue_groups == "General referral" & between(age_days, 30, 364) ~ str_c(as.character(as.integer(age_days/30.437)), " months"), # Avg days in a month 
                               tissue_groups == "General referral" & age_days <= 30 ~ "< 1 month",
                                TRUE ~ "Age not revelant for tissue type"), 
         #  Create a new column, turn numbers into strings they represent
         pathogenic_string = case_when(Pathogenic == 1202218788 ~ "Pathogenic", 
                            Pathogenic ==  1202218783 ~ "Likely pathogenic", 
                            Pathogenic == 1202218781 ~ "Uncertain significance"),
         # Group changes into cnv groupings 
         cnvtype = case_when(Copies == 1190384936 | Copies == 1190384935 ~ "Copy number loss",  # change = 'x1 OR x0'
                             Copies == 1190384938 | Copies == 1190384940 | Copies == 1190384941 |
                             Copies == 1190384942 ~ "Copy number gain", # change = 'x2 OR  x3 OR x4 OR other. Sally Walsh confirmed all 'others' to go in with the gains 
                             Copies == 1190384937 | Copies == 1190384943 ~ "Mosaic loss", # change = 'x0~1 OR x1~2'
                             Copies == 1190384939 | Copies == 1190384949 ~ "Mosaic gain"), # change = 'x2~3 OR x2~4'
         # Create colours for UCSC depending on CNV type and pathogenic
         itemRgb = case_when((cnvtype == "Copy number loss" & pathogenic_string == "Pathogenic") ~ "139,0,0", # Dark red
                             (cnvtype == "Copy number loss" & pathogenic_string == "Likely pathogenic") ~ "255,0,0", # Red
                             (cnvtype == "Copy number loss" & pathogenic_string == "Uncertain significance") ~ "240,128,128", # Light red
                             (cnvtype == "Copy number gain" & pathogenic_string == "Pathogenic") ~ "0,0,128", # Dark blue
                             (cnvtype == "Copy number gain" & pathogenic_string == "Likely pathogenic") ~ "0,191,255", # Blue
                             (cnvtype == "Copy number gain" & pathogenic_string == "Uncertain significance") ~ "173,216,230", # Light blue
                             (cnvtype == "Mosaic loss" & pathogenic_string == "Pathogenic") ~ "50,23,77", # Dark purple 
                             (cnvtype == "Mosaic loss" & pathogenic_string == "Likely pathogenic") ~ "143,0,255", # Purple 
                             (cnvtype == "Mosaic loss" & pathogenic_string == "Uncertain significance") ~ "230,230,250", # Light purple
                             (cnvtype == "Mosaic gain" & pathogenic_string == "Pathogenic") ~ "18,53,36", # Dark green  
                             (cnvtype == "Mosaic gain" & pathogenic_string == "Likely pathogenic") ~ "141,182,0", # Green 
                             (cnvtype == "Mosaic gain" & pathogenic_string == "Uncertain significance") ~ "208,240,192", # Light green 
                             ),  
        # Create the detailed string of patient information to go into UCSC
        details_string = str_c("<br /><br /><strong>", PatientResult, "</strong><br /><br />", # Add line breaks
                                  DESCtwo, "<br /><br />", "Tissue group: ", tissue_groups ,"<br /><br />", 
                                  "Inheritance: " , PatientIDInheritance,  "<br /><br />", 
                                  "Sex: ", sex , "<br /><br />", "Age: ", age_years_days, "<br /><br />",
                                  "CNV group: "  , cnvtype, "<br /><br />", "Pathogenicity: ", pathogenic_string,
                                  "<br /><br />" , "Phenotypes: ", phenotype_string,"<br /><br />"),
        # Create columns for UCSC requirements 
         start192 = Start19, # Some columns need to be repeated (and used for ID when re joining tables)
         stop192 = Stop19,
         patientid_2 = PatientID,
         strand = ".", # Create strand column for bed file requirement
         score = 0, # blank score column 
         ID = row_number()) %>% # create a unique ID column to join this data frame to data after lift over 
  filter(sex != "blank")  # Remove patients without a sex (64 results over 44 patients are removed at this step)


# Make a new df ready to export for liftover via command line 
bed_track <- array_df %>% 
  select(Chromo, Start19, Stop19, PatientID, start192, stop192, ID)  # Select needed columns

### Save data for liftover to hg38 ######### ====================================================================
# Liftover software doesn't like the variable size of the phenotype/detailed string columns so minimal information exported for liftover 
# To be re merged with patient information after liftover 
# Write df to file 
write.table(bed_track, "array_data_to_liftover_20220215.bed", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

### Run liftover  ######### ====================================================================

# Run liftover via command line, newFile is the lifted over results
# over.chain date at download hg19ToHg38.over.chain.gz 2013-12-31 23:08 
# liftOver version downloaded from UCSC website 2022-01-27, no version listed 
# bedPlus=3 File is bed N+ format (i.e. first N fields conform to bed format)
# Command run 
# /home/erin/Documents/Software/UCSC_liftover/liftOver array_data_to_liftover_20220215.bed /home/erin/Documents/Software/UCSC_liftover/hg19ToHg18.over.chain 
# newFile_20220215 unMapped_20220215 bedPlus=3

### Load in liftover  data ######### ====================================================================

# Import lifted over data 
array_liftover_done <- read.delim("newFile_20220215", sep = '\t'  , header = FALSE)

### Merge lift over data with patient information  ######### ====================================================================

merged_lifted_over <- array_liftover_done %>% 
  rename("hg38_chr" = "V1", # Change headers for lifted over hg38 coordinates
         "hg38_str" = "V2",
         "hg38_stp" = "V3",
         "PatientID" = "V4",
         "start192" = "V5", # change headers to old hg19 coordinates for ID whilsts joining 
         "stop192" = "V6", 
         "ID" = "V7") %>% 
  inner_join(array_df,  by = c("PatientID", "start192", "stop192", "ID")) %>% # add phenotype / patient details back on, match by start(hg19) & stop(hg19) coordinates, PRU and ID
  mutate(hg38_str_2 = hg38_str,
         hg38_stp_2 = hg38_stp)  # Create duplicate str and stp for UCSC 


### Sub select dfs for exporting  ######### ====================================================================

# Prepare final df for upload to UCSC for whole chromosome CNVs (to include dups and dels)
ucsc_bedfile_whole_chromo <- merged_lifted_over %>% 
  ungroup() %>% 
  filter(WholeChromosome == 1202218774) %>% # whole chromosome = yes 
  select(hg38_chr, hg38_str, hg38_stp, PatientID, score, strand, hg38_str_2, hg38_stp_2, itemRgb, patientid_2 ,details_string) # Columns needed for UCSC

# Prepare two dfs for upload to UCSC for non whole chromosome CNVs
ucsc_bedfile_cnv_gains_only <- merged_lifted_over %>% 
  ungroup() %>% 
  filter(WholeChromosome == -1128799521 | is.na(WholeChromosome)) %>% # whole chromosome = no OR is blank 
  filter(cnvtype == "Copy number gain" | cnvtype == "Mosaic gain" ) %>% 
  select(hg38_chr, hg38_str, hg38_stp, PatientID, score, strand, hg38_str_2, hg38_stp_2, itemRgb, patientid_2 ,details_string) # Columns needed for UCSC

ucsc_bedfile_cnv_losses_only <- merged_lifted_over %>% 
  ungroup() %>% 
  filter(WholeChromosome == -1128799521 | is.na(WholeChromosome)) %>% # whole chromosome = no or is blank 
  filter(cnvtype == "Copy number loss" | cnvtype == "Mosaic loss" ) %>% 
  select(hg38_chr, hg38_str, hg38_stp, PatientID, score, strand, hg38_str_2, hg38_stp_2, itemRgb, patientid_2 ,details_string) # Columns needed for UCS

### Export data to be uploaded to UCSC  ######### ====================================================================

write.table(ucsc_bedfile_whole_chromo, "ucsc_array_cnv_moka_legacy_whole_chromo_20220215.bed", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE )

write.table(ucsc_bedfile_cnv_gains_only, "ucsc_array_cnv_moka_legacy_cnv_gains_only_20220215.bed", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE )

write.table(ucsc_bedfile_cnv_losses_only, "ucsc_array_cnv_moka_legacy_cnv_losses_only_20220215.bed", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE )




