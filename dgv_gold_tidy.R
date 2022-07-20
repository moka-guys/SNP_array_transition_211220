library(tidyverse) # V1.3.1

########
# This script has been created to tidy DGV Gold Standard Variants 
# downloaded from http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3 release date 2016-05-15
# to be put into custom UCSC tracks to be used by the array team 
########

##### Load data  ####### ====================================================================

setwd("/home/erin/Documents/Work/SNP_array_liftover/dgv_gold_track_38") 
dgv_gold_df <- read.delim("DGV.GS.hg38.gff3", sep = '\t' , header = FALSE) 

### Tidy data ######### ====================================================================

# Column names for separate below
col_names <- c("ID", "Name", "variant_type", "variant_sub_type", "outer_start", "inner_start", "inner_end", "outer_end", "inner_rank", "num_variants",
               "variants", "num_studies", "Studies", "num_platforms", "Platforms", "number_of_algorithms", "algorithms","num_samples", "samples", 
               "Frequency","PopulationSummary","Number_of_unique_samples_tested")

# Column names for pivout longer, with "ID" & "variant_sub_type" removed, so each row is uniquely identified 
col_names_two_less <- c( "Name", "variant_type",  "outer_start", "inner_start", "inner_end", "outer_end", "inner_rank", "num_variants",
               "variants", "num_studies", "Studies", "num_platforms", "Platforms", "number_of_algorithms", "algorithms","num_samples", "samples", 
               "Frequency","PopulationSummary","Number_of_unique_samples_tested")


tidy_dgv_gold_df <- dgv_gold_df %>% 
  separate(V9, col_names, sep = ";") %>% # Separate the final column into multiple columns
  pivot_longer(all_of(col_names_two_less), names_to = "name", values_to = "value") %>% # Flip df to easily manipulate data
  mutate(value_alone = str_remove(value, ".*=")) %>% # Remove string up to "="
  select(-c(value)) %>% # Drop unnecessary column 
  pivot_wider(names_from = name, values_from = value_alone) %>% 
  mutate(variant_sub_type = str_remove(variant_sub_type, ".*=")) %>% 
  rename(chrom = V1, # rename 
         chromStart = outer_start, # Array team have requested to have both 
         chromEnd = outer_end, 
         thickStart = inner_start,
         thickEnd = inner_end) %>% 
  mutate(calc = str_c(num_samples, Number_of_unique_samples_tested, sep = "/"), # Make new columns for UCSC, some are required in bed file format
         frequency_round = round(((as.integer(num_samples)/as.integer(Number_of_unique_samples_tested))*100), digits = 1),
         Frequency = str_c(frequency_round, "%"),
         ID = str_c(Frequency, "(",calc, ")"),
         ID_2 = ID, # Variable duplicated for UCSC requirements for bedDetail format
         strand = ".",
         score = 0,  
details_string = str_c("<br /><h3>  Frequency ", Frequency, " (", num_samples,  " out of ", Number_of_unique_samples_tested,
                       " unique samples tested) </h3><a href=http://dgv.tcag.ca/gb2/gbrowse_details/dgv2_hg38?name=", Name ,
                       ">DGV Gold call ID" , Name ," </a> (click to go to DGV)<br /><br />DGV Gold Standard Variants Release Date 2016-05-15.<br />"),
itemRgb = case_when(variant_sub_type == "Gain" ~ "47,63,235", # Blue for gain
                    variant_sub_type == "Loss" ~ "255,0,0")) %>% # Reds for loss
  select(chrom, chromStart, chromEnd, ID, score, strand, thickStart, thickEnd, itemRgb, ID_2, details_string) %>% # Select required columns
  distinct() # Some variants have inner and outer coordinates from the raw data, leading to duplicated entries & this distinct removes the duplicates 


# One track needed per variant type 
tidy_dgv_gold_gain <- tidy_dgv_gold_df %>% 
  filter(itemRgb == "47,63,235") # df of gains only 

tidy_dgv_gold_loss <- tidy_dgv_gold_df %>% 
  filter(itemRgb == "255,0,0") # df of losses only 

### Save data to be uploaded to UCSC  ######### ====================================================================

write.table(tidy_dgv_gold_gain, "ucsc_hg38_DGVgold_gain_2016-05-15.bed", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(tidy_dgv_gold_loss, "ucsc_hg38_DGVgold_loss_2016-05-15.bed", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)



