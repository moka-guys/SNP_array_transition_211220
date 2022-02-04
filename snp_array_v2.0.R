library(tidyverse) # V1.3.1
library(bedtoolsr) #  V2.30.0.1 bedtools R package 
options(bedtools.path = "/") # Set new path to bedtoools 
# Bedtools Linux version 2.26.0


####################
# This script has been built to find out how many probes in the new HT-CMA array bed file 
# overlap with a rolling three probe window from the old CGH array bed file 
# 1) Import data 
# 2) Filter the old bed file to the chromosome being looped over, order it by start 
# 3) Create a slice of the old bed file which is a rolling 3 base pair window
# 5) Run bedtools intersect to find the probes from the new array which overlap with this window 
# 6) Count the number of rows in the intersected df, 
#    which is the number of new probes which can be found in the current old 3 probe window 
# 7) When then end of the loop is reached add this set of numbers to a list 
# 8) Bind all the row counts for each chromosome into one data frame 
##################



### Importing data etc ####
setwd("/")
old_array <- read.delim("CGH_4.3_085030_D_BED_20170809.bed_hg38_liftover.bed", sep = '\t' , header = FALSE)  # Load data 
new_array <- read.delim("HT-CMA-hg38.bed", sep = '\t' , header = FALSE)  
print("data loaded")

names(old_array)[1:4] <-c("chr", "start", "end", "probe") # rename columns

names(new_array)[1:4] <-c("chr", "start", "end", "probe") 

chromosomes_list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                      "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                      "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                      "chr20", "chr21", "chr22","chrX", "chrY") # list of chromosomes to be looped over

df_list = list() # Empty list to sort dfs in 

setwd("/output")


## Processing loop ##
print("starting loop!")
for (chromosome in chromosomes_list) { # filter to one chromosome at a time to speed up intersect
  
  print(paste0("Starting ", chromosome))
  
  row_counts <- integer() # Make a new list for each chromosome loop
  
  for (slice_row in 1:4953) { # Run the loop for the maximum number of rows (found in chr1) 
    
    sliced_df <-  old_array %>% 
      filter(chr == chromosome) %>% # Filter the old df to the current chromosome 
      arrange(start)  %>% # Arrange in size order 
      slice(slice_row:(slice_row+2)) %>%  # Slice the row which the counter is on, then add three to get window
      mutate(start_new = min(start), # Take the smallest number from the start column
             end_new = max(end)) %>% # Largest from the end, this becomes out window 
      select(chr, start_new, end_new) %>% # Make this into our temp bed files 
      distinct()
    
    new_array_filtered <- new_array %>% 
      filter(chr == chromosome) %>% # filter new array to the current chromosome 
      arrange(start) # Arrange in size order 
    
    
    intersect <- bedtoolsr::bt.intersect(new_array_filtered, sliced_df) # Intersect the new array by the slice of the old array
    count_rows <- as.integer(count(intersect)) # count the number of rows in the new array bed file which overlap with our 3 probe old array slice
    row_counts[[length(row_counts) + 1]] <- count_rows  # add this count of rows into a list
  }
  print(paste0("Chromsome being added to df: ", chromosome))
  
  df <- data.frame(row_counts) # Make the list of row counts into a df
  names(df)[1] <- chromosome # Make the header of the column the chromosome name
  df_list[[chromosome]] <- df # add df to list outside of the loop 
  print(paste0("Chromsome finished: ", chromosome))
  
  
}

## Binding df and saving output ## 
print("Loop finished, making dataframe")

df_bound <- dplyr::bind_cols(df_list) # Bind all the dfs in the list into one df 
write.table(df_bound , "df_means_211218.txt", sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE) # Save df 

print("Saved means DF")
