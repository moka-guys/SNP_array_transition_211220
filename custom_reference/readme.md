# Custom Reference

This folder contains two scripts to create custom reference files, for the ChAS SNP array platform
Custom reference files (which are a bit like a panel of normals) are used to normalise the data and remove noise. Using a custom reference made from our own data can help improve the calling.

## file\_mover\_postnatal.py

This script takes a spreadsheet where one file contains a list of specimen numbers and copies CEL files from some hardcoded folders into the specified output folder.

The script can be run from either the S Drive or the Arrays VM (access added to the cloud by Synnovis IT). A pop up box will appear where the user can select
where they are running the script from. This will choose the correct hard-coded CEL file folder locations specified within the script.

The script will then look recursively through hardcoded folders for CEL files containing the specimen number. It will copy those to sex specific folders in the directory given to --output_folder. The script takes into account any duplicate files (i.e. if a file exists in multiple folder locations, it will only copy a single copy of each file). Additionally, if any specimen numbers are seen in multiple CEL files, all CEL files for those specimen numbers will be excluded from copying.

### Spec Numbers CSV File

The CSV file given as an argument to --spec_number_file should contain specimen numbers for all the files that need to be moved. This file should have a header row, with one column named "Specimen ID" (case sensitive) and the script will extract all spec numbers from that column.

Additionally, if there is also a column named "Result Type" then the script will attempt to decode the sex from this to seperate CEL files into male and female subfolders within the output location. If this column is empty, or the sex can't be determined from the content of this column the file will be saved into an "undetermined" subfolder.

* Specimen numbers must be in the form 22-12345 (with "-" rather than "/")

### CSV Spec Numbers to Exclude CSV File

If a CSV containing spec numbers to exclude was given as an input (--exclude_spec_numbers) then this file will also be parsed, looking for a column named "Specimen ID" (case sensitive) and extract this list. This list would be used to remove any matching specimen numbers from the --spec_numbers input.

* Specimen numbers must be in the form 22-12345 (with "-" rather than "/")

### Usage

The script can be run on either the S Drive or the arrays VM.

```bash
usage: file_mover_postnatal.py [-h] --output_dir OUTPUT_DIR --spec_number_file SPEC_NUMBER_FILE
                               [--exclude_spec_numbers_file EXCLUDE_SPEC_NUMBERS_FILE]

This script is used to collate CEL files for a list of spec numbers to create a custom reference file

options:
  -h, --help            show this help message and exit
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output folder to copy CEL files to
  --spec_number_file SPEC_NUMBER_FILE, -s SPEC_NUMBER_FILE
                        Path to CSV file where one column contains spec numbers
  --exclude_spec_numbers_file EXCLUDE_SPEC_NUMBERS_FILE, -e EXCLUDE_SPEC_NUMBERS_FILE
                        Path to CSV file where one column contains spec numbers to exclude
```

**N.B. the output location in the below commands should be changed each time the script is run**

#### S Drive

```
S:\Genetics_Data2\Array\Software\Python-3.6.5\python.exe "S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud\Software\SNP_array_transition_211220\custom_reference\file_mover_postnatal.py" -o "S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud\Custom reference project\custom_reference_2024" -s "S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud\Custom reference project\custom_reference_2024\GW list_all unfiltered_1 Jan to 30 Sep 2023_NEW_normal postnatal female and male only_SPECIMEN LIST.csv" -e "S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud\Custom reference project\custom_reference_2024\DMD list for custom reference exported 07022024.csv"
```

#### VM

```
\\GRPVCHASDB01\Genetics\In\Software\Python-3.6.5\python.exe "\\GRPVCHASDB01\Genetics\In\Software\SNP_array_transition_211220\custom_reference\file_mover_postnatal.py" -o "\\GRPVCHASDB01\Genetics\In\Custom reference project\custom_reference_2024" -s "\\GRPVCHASDB01\Genetics\In\Custom reference project\custom_reference_2024\GW list_all unfiltered_1 Jan to 30 Sep 2023_NEW_normal postnatal female and male only_SPECIMEN LIST.csv" -e "\\GRPVCHASDB01\Genetics\In\Custom reference project\custom_reference_2024\DMD list for custom reference exported 07022024.csv"
```

## file\_mover\_prenatal.py
This script is a bit more complex as prenatal samples are only analysed in certain regions. It is run in two parts:

1) take a list of specimen numbers, find the rhchp files and copy them into a subfolder. These can be used to create a msv to be used in step 2:
2) for the list of specimens, if they do not have a call within the a known syndrome region copy them into a subfolder


### Usage

The script should be run on the arrays VM to ensure that all folder locations are accessible to the script.


### step 1 - Collect rhchp files:
need to provide 3 arguments:

--output_folder','-o', output folder to copy files to
--input_folder','-i', folder to search for rhchp files
--spec_numbers','-s', csv file where one column contains spec numbers

This will loop through the csv file (comma seperated) containing specimen numbers. This should have a header row, with one column named "Specimen ID" (case sensitive) and will extract all spec numbers from within the rest of the file,
It will then look recursively through the --input\_folder path and copy to the --output_folder directory

### step 2 - identify files that can be used to make a custom reference
repeats step 1 (these steps include checks that files are present so shouldn't take much longer, or result in duplicate files as long as the same argument are provided).
Requires 3 additional arguments (optional - means if not provided only step 1 will be run):

--syndrome_regions,'-r', optional, BED file containing syndromic regions
--syndrome\_free\_files, '-f', optional , folder to contain anonymised rhchp files that do not have calls overlapping with syndromic regions
--multi\_sample\_viewer_output','-m',optional, output of multisample viewer containing calls in multiple samples

With the additional 3 arguments the script will then assess all the calls from the multisample viewer file. 
For each sample, if there are no calls which overlap the syndromes provided in the syndrome BED file (tab seperated in format chr, start, stop , type (Gain/Loss)) then it can be used to make the custom reference. 
To make a custom reference .CEL files are required. The script looks in (hardcoded) directories for the file, it will anonymise the file and copy into the directory provided to the --syndrome\_free\_files argument.
