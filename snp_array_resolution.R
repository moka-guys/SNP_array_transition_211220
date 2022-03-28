library(tidyverse) # V1.3.1
source("config.R") # This file contains the probe_window variable which can be changed
####################
# This script has been built to find the genomic resolution of the new HT-CMA array bed file 
# over a number of different probe windows
# 1) Import data 
# 2) Create a slice of the bed file which is a rolling N (users choice) probe window
  # value taken from the config.R file 
# 3) Count the number of base pairs in each N probe window
# 4) Create an average base pair range, across the genome
# 5) Save to df to a txt file
##################

## Import data ##
setwd("/")
snp_array_bed <- read.delim("HT-CMA-hg38.bed", sep = '\t' , header = FALSE)  
print("data loaded")

names(snp_array_bed)[1:4] <-c("chr", "start", "end", "probe") # rename columns

chromosomes_list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                      "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                      "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                      "chr19", "chr20", "chr21", "chr22","chrX", "chrY") # ordered chromosomes for a factor

# Ensure dataframe is in chromosome and then ascending order by start
# Make chromosome column into a factor so it sorts in the correct order 
sorted_snp_array_bed <- snp_array_bed %>% 
  mutate(chr = factor(chr, levels = chromosomes_list, ordered = T)) %>% # make sure choromosomes are in the correct order 
  group_by(chr) %>% 
  arrange(chr, start)  # order by chromosome, then by start 

df_list = list() # Empty list to sort dfs in 
df_list_genome  = list()
## Processing loop ## ====================================================================

print("starting loop!")
for (chromosome in chromosomes_list) { # filter to one chromosome at a time
  
  df_list = list() # New list for each chromosome
  start_row <- 1 # Create the first rolling window, but start fresh for every new chromosome
  end_row <- probe_window # end row is equal to the probe_window being investigated 

  print(paste0("Starting ", chromosome))
  
  filtered_df <-  sorted_snp_array_bed %>%
    filter(chr == chromosome) %>% # Filter the old df to the current chromosome
    arrange(start) # order
  
  count_rows_stop <- as.integer(nrow(filtered_df)) # get a count of the rows to stop the loop

  while (start_row<(count_rows_stop-probe_window)){  # Stop before the last probe window 
    
    df <- filtered_df %>% 
      slice(start_row:end_row) %>%  # Slice the df by the rows in the current window
      mutate(start_min =min(start),
             stop_max = max(end),
             bp_size = stop_max - start_min) %>% 
      select(chr, bp_size) %>% 
      distinct()
     
     start_row <- start_row+1 # Progress on to the next window
     end_row <- end_row+1
     
     df_list[[start_row]] <- df # Add this mean bp for this window to the list
  } 
  
  print(paste0("Chromsome being added to df: ", chromosome))
  
  df_bound_chr <- dplyr::bind_rows(df_list) %>%  # Bind all the columns for this chromosome together
    mutate(average_bp = mean(bp_size)) %>% # Find the mean bp_size per probe window
    select(chr, average_bp) %>% 
    distinct() # Return the average bp and chromosome on one row
   
  df_list_genome[[chromosome]] <- df_bound_chr # add the result to the df 
}

## Binding df and saving output ## ========================================
print("Loop finished, making dataframe")

df_bound_genome <- dplyr::bind_rows(df_list_genome) %>%  # Bind all the dfs in the list into one df 
  ungroup()  %>%
  mutate(overall_mean = mean(average_bp))  

mean(df_bound_genome$average_bp) # Print out to command line 

setwd("/output")

write.table(df_bound_genome , (paste0("new_array_genome_resolution_",probe_window,"_probe_window.txt")), sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE) # Save df 
print("df saved")

