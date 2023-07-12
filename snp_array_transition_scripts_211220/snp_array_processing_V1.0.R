library(tidyverse) # V1.3.1
########################
# This script has been created to find the range, mean and other stats from the 
# output of the snp_array_v2.0R script 
# 1) Takes the input of a csv file  
# 2) Pivots the data  
# 3) Creates stats 
# 4) Saves the output
##########################

setwd("/home/erin/Documents/Work/SNP_array_liftover")
means_df <- read.delim("df_means_211218.csv", sep = '\t' , header = TRUE)   # Import in data to be processed 

means_df_longer <- means_df %>% 
  pivot_longer(cols = c("chr1",  "chr2",  "chr3"  ,"chr4"  ,"chr5" , "chr6" , "chr7",  "chr8",
                        "chr9" , "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16",
                        "chr17" ,"chr18" , "chr19", "chr20", "chr21" ,"chr22", "chrX" , "chrY") , 
               names_to = "chromosome", values_to = "num_probes") %>% # Pivot data frame longer for easier calculations 
  group_by(chromosome) %>% # Calculate per chromosome 
  mutate(lower_range = min(num_probes, na.rm = TRUE), # Get stats, ignoring NA values
         upper_range = max(num_probes, na.rm = TRUE),
         mean = signif(mean(num_probes, na.rm = TRUE),3),
         median = median(num_probes, na.rm = TRUE),
         chromosome = factor(chromosome, levels = ordered_chrs, ordered = TRUE),
         less_than_three = sum(num_probes <= 3,  na.rm = TRUE),
         row_count = sum(!is.na(num_probes)),
         perc_less_than_three = signif((less_than_three/row_count)*100,3)) %>% 
  select(-c(num_probes)) %>% # Remove unneeded columns 
  distinct() %>%  # Select one row per column 
  arrange(chromosome)


write.table(means_df_longer, "means_range_211218.txt", sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE) # Save as text file 

