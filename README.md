# SNP_array_transition_211220
Scripts facilitating SNP array transition to a new platform including pulling legacy from Moka and liftover of that data to GRCh38. 

In late 2021 the array service at Viapath started a transition from a CGH4.3 platform on GRCh37 to CytoScan™ HT-CMA on GRCh38. A number of scripts were required to facilitate this transition and they are held in this repo. 

# Scripts for comparing old probes and new probes bed files

These scripts were run on the backup workstation using a docker version of R with tidyverse:4.1.1
 
Run snp_array_v2.0.R
Run snp_array_processing_V1.0.R

# Scripts for pulling legacy data from MOKA

Run using python2.7

moka_array_cnv_export_V1.0.py

# Scripts for lifting over legacy data to GRCh38 

To be added 
