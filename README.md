# SNP_array_transition_211220
Scripts facilitating SNP array transition to a new platform including pulling legacy from Moka and liftover of that data to GRCh38. 

In late 2021 the array service at Viapath started a transition from a CGH4.3 platform on GRCh37 to CytoScan™ HT-CMA on GRCh38. A number of scripts were required to facilitate this transition and they are held in this repo. 

## Compare coverage of old probes and new probes 

Scripts for comparing old probes and new probes bed files. Run using Docker on the backup workstation.

rocker/tidyverse:4.1.2: snp_array_v2.0.R

rocker/tidyverse:4.1.2: snp_array_processing_V1.0.R

## Legacy data export from Moka 

#### 1 ) Filtering requirments for data 
- Each patient has to meet one of the criteria in each row to be included in the legacy data set(Query run 11/02/2022)

|Column being filtered  | Filters 
|-----------------------| ---------------------------------------------------------------|
| Referral types IS     |  POC/Tissues **OR** General referral **OR** Prenatal **OR** Repeat POC **OR** Repeat general referral **OR** Repeat prenatal **OR** Pseudorush 
| Pathogenicity IS      | Uncertain clinical significance (class 3) **OR** Likely to be pathogenic (class 4) **OR** Pathogenic / abnormal result (class 5)
| CNV/ Confirmation (from Array  Results page in Moka)  IS  | Reported **OR** Not reported (insufficient evidence) **OR** Not reported (unrelated to referral)
| Array Test Overall Result IS | Completed 
| Overall Results Code IS NOT * | Failed / testing not possible **OR** Preliminary result
| Array Date received (when the batch of slides used was received) IS | Between '2016-01-01 00:00:00' and '2021-12-31 00:00:00'
  
*The overall results statuses of the results included after this filtering are: Abnormality detected / Abnormality detected, inherited / Abnormality detected, de novo / Abnormality detected, inferred inherited / Abnormality detected, inferred de novo / Normal 

#### 2) Script for pulling legacy data from MOKA

python(v2.7): moka_array_legacy_data_export.py

#### 3) Script for tidying legacy data to be uploaded to UCSC & facilitating liftover to hg38

R(v4.0.5): moka_array_legacy_data_export_tidy.R
