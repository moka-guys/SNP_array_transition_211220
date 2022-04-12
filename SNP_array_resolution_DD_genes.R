library(tidyverse) # V1.3.1
library(dplyr) # V1.0.8
library(biomaRt) # V2.46.3
library(bedtoolsr) # V2/30.0-1
source("config.R") # This file contains the probe_window variable which can be changed

####################
# This script has been built to find the genomic resolution of developmental delay (DD) genes 
# in the new HT-CMA array bed file 
# Genes list taken from Decipher downloaded as DDG2P (DDG2P_16_3_2022.csv)
# Containing 2539 genes 

## Script actions ##
# 1) Import bed file and DD csv 
# 2) Use BSgenome.Hsapiens.UCSC.hg38 to get genomic coordinates for DD genes  
# 3) Intersect the DD genes coordinates with the probes bed file
  # remove the genes which don't meet the window for probes specified 
# 4) Create a slice of the bed file which matches each genes coordinates 
# 5) Intersect a X probe window with the gene coordinates
 # Record the resolution of the probes
# 6) Move one alone to the next probe window, repeat step 5
# 7) When all genes have been processed, create a mean and median resolution per X windows 
##################

## Import data & make lists =============================================================

#setwd("/") # for use with docker image
setwd("/home/erin/Documents/Work/SNP_array_liftover/bed_comparison")
snp_array_bed <- read.delim("HT-CMA-hg38.bed", sep = '\t' , header = FALSE) # load in new array bed file 

names(snp_array_bed)[1:4] <-c("chr", "start", "end", "probe") # rename 

setwd("/home/erin/Documents/Work/SNP_array_liftover")
dd_uscs_csv <- read.delim("DDG2P_16_3_2022.csv", sep = ',' , header = TRUE)  # load in DD genes csv
print("data loaded")

chromosomes_list <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6",
                      "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                      "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                      "chr19", "chr20", "chr21", "chr22","chrX", "chrY")  # list for loop later 


### Getting coordinates of DD genes ==================================================

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # hsapiens_gene_ensembl = Human genes (GRCh38.p13) GRCh38.p13

dd_uscs_csv_hgnc_id_symbol <- dd_uscs_csv %>% # Take the dd csv df and filter to only HGNC ID and gene symbol 
  dplyr::select(hgnc.id, gene.symbol) %>% # there's a select in biomart too
  distinct() # filter to only return one entry for each gene

dd_uscs_csv_hgnc_id_symbol$hgnc.id <- paste("HGNC:", dd_uscs_csv_hgnc_id_symbol$hgnc.id, sep="") # Add HGNC: to search Ensembl database

# Query Ensembl database with DD genes list 
dd_ensembl_coordinates <- getBM(attributes = c('hgnc_id','hgnc_symbol', 'chromosome_name',
                                      'start_position', 'end_position'), # Return these values
                       filters = 'hgnc_id',
                       values = dd_uscs_csv_hgnc_id_symbol$hgnc.id, # Search by this column
                       mart = ensembl) # Use this database

## Tidying data ## ==============================

join_uscs_ensembl_dd <- dd_ensembl_coordinates %>% 
  inner_join(dd_uscs_csv_hgnc_id_symbol, by=c('hgnc_id'='hgnc.id')) # Join coordinates back to csv


  # some chromosomes returned are "CHR_HSCHR6_MHC_DBB_CTG1"
  # this section of code pulls out the number from the chromosome 
  # 30 genes are lost in this code as they have no chr returned
  greater_three_chr <- join_uscs_ensembl_dd %>% 
    filter(nchar(chromosome_name)>2) %>% # Bring back rows where chr has > 2 characters (e.g more than 'chr')
    group_by(hgnc_id) %>% 
    mutate(chr = str_split(chromosome_name, "_")[[1]][[2]]) %>% # Split the string at '_' 
    filter(str_detect(chr, regex("CHR"))) %>% # Find all the columns which have CHR
    mutate(chromosome_name = str_split(chr, "CHR")[[1]][[2]]) %>% # Split this again to get the chromosome number 
    dplyr::select(-c(chr)) 
 
  
ensemble_dd_genes_df <- join_uscs_ensembl_dd %>% # join our tidy df back 
  bind_rows(greater_three_chr) %>% 
  filter(nchar(chromosome_name)<=3) %>% 
    mutate(bp_size = end_position - start_position) %>%  
    group_by(hgnc_id) %>%  
    top_n(1, bp_size) %>% # Take the largest transcript in the case of duplicates
    ungroup() %>%  
    distinct(hgnc_symbol, bp_size, .keep_all=TRUE) %>%  # Remove duplicates where the coordinates were different but bp size was the same 
  dplyr::select(chromosome_name, start_position, end_position, hgnc_symbol) # Re order 

ensemble_dd_genes_df$chromosome_name <- paste("chr", ensemble_dd_genes_df$chromosome_name, sep="") # Add chr to the chromosome col

names(ensemble_dd_genes_df)[1:3] <-c("chr", "start", "end") # rename cols

#### Finding resolution of DD genes ============================================================

## Intersect to get probe count ### 

intersect <- bedtoolsr::bt.intersect(snp_array_bed, ensemble_dd_genes_df, wa = TRUE, wb = TRUE) # wa keeps all the cols from snp_bed, wb from DD genes

names(intersect)[1:8] <-c("chr_probe", "start_probe", "end_probe", "name_probe",
                          "chr", "start", "end", "name_gene") # rename

intersect_counts <- intersect %>% 
  group_by(name_gene) %>% 
  mutate(probe_count = n()) %>% # Get a count of probes in each gene 
  dplyr::select(chr, start, end,  name_gene, probe_count) %>% 
  distinct()
  
# Make a df of genes which have < X probes listed in config 
genes_less_than_probes <- intersect_counts %>% 
  filter(probe_count < probe_window) 

# Save these for later inspection 
write.table(genes_less_than_probes , (paste0("DD_genes_less_than_",probe_window,"_probe_window.txt")), sep ='\t', 
            col.names = TRUE, row.names = FALSE, quote = FALSE) 
print("less than probes genes df saved")


## X window probe resolution ## ===========================

gene_only_intersect_counts_filtered <- intersect_counts %>% # make a df of only the genes, with a no of probes > window 
  filter(probe_count >= probe_window) %>% 
  dplyr::select(chr, start, end, name_gene) %>% 
  distinct()
  

## Loop! ## ===================================================

df_list_genome = list()

count <- 1 

## Processing loop ##
print("Starting loop!")
for (chromosome in chromosomes_list) { # filter to one chromosome at a time to speed up intersect
  
  print(paste0("Starting ", chromosome))
  
  df_list = list() # Empty list to sort dfs in 
  
  dd_gene_df <-  gene_only_intersect_counts_filtered %>%
    filter(chr == chromosome) %>% # Filter the gene df to one chromosome 
    arrange(start)
  
  start_window <- 1 # Create the first rolling window
  end_window <- probe_window # end row is equal to the probe_window being investigated 
  
  for (i in 1:nrow(dd_gene_df)) { # For each gene in the filtered_df
    
    row <- dd_gene_df[i,] # Get the gene as a single row 
    
    print(paste0("Processing gene ", row$name_gene))
    slice_df_start <- row$start # Start coordinate of the gene 
    slice_df_end <- row$end # End coordinate of the gene 
    
    snp_array_bed_filtered <-  snp_array_bed %>% # filter the array bed file to match the limits of the gene being investigated 
      filter(chr == chromosome,
             start >= slice_df_start,
             end <= slice_df_end) %>% 
      arrange(start)
    
    count_rows_stop <- as.integer(nrow(snp_array_bed_filtered)) # get a count of the rows to stop the loop
    
    while (end_window <=(count_rows_stop)) { # Stop the loop when the row is the size of the final window 
      
      snp_array_bed_sliced <-  snp_array_bed_filtered %>% # slice the array bed file  to the window 
        slice(start_window:end_window)
      
      # Intersect the array bed file with the dd gene row, wb keeps information from DD genes
      intersect_dd_genes <- bedtoolsr::bt.intersect(snp_array_bed_sliced, dd_gene_df, wb =TRUE) 
      
      names(intersect_dd_genes)[1:8] <-c("chr_probe", "start_probe", "end_probe", "name_probe",
                                         "chr_gene", "start_gene", "end_gene", "name_gene") # rename 
      
      
      snp_array_res_per_gene <- intersect_dd_genes %>% 
        mutate(start_min =min(start_probe), 
               stop_max = max(end_probe),
               bp_size = stop_max - start_min) %>% # Get the bp range for the probes 
        dplyr::select(chr_probe, name_gene,  bp_size) %>% 
        distinct()
      
      count <- count+1 # Make each addition to the list unique to stop it being overwritten 
      
      df_list[[count]] <- snp_array_res_per_gene # add into the list 
      
      
      start_window <- start_window+1 # Move window along 
      end_window <- end_window+1
    }
  
  df_bound_chr <- dplyr::bind_rows(df_list) # Make the chromosome into a df 

  df_list_genome[[chromosome]] <- df_bound_chr # Put the chromosome df into another list 
  }
  
  print(paste0("Chromsome being added to df: ", chromosome))
}

## Binding df and saving output ## =====================================
print("Loop finished, making dataframe")

df_bound_genome <- dplyr::bind_rows(df_list_genome) %>%  # Bind all the dfs in the list into one df 
  ungroup()  %>%
  mutate(overall_mean = mean(bp_size),
         overall_median = median(bp_size)) %>% 
  select(overall_mean, overall_median)


#setwd("/") # for use with docker image
write.table(df_bound_genome , (paste0("new_array_dd_genes_resolution_",probe_window,"_probe_window.txt")), sep ='\t', col.names = TRUE, row.names = FALSE, quote = FALSE) # Save df 
print("df saved")

