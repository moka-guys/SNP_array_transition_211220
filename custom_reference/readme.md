# Custom Reference
This folder contains two scripts to create custom reference files, for the ChAS SNP array platform
Custom reference files (which are a bit like a panel of normals) are used to normalise the data and remove noise. Using a custom reference made from our own data can help improve the calling.

## file\_mover\_postnatal.py
This script takes a spreadsheet where one file contains a list of specimen numbers and copies CEL files from two hardcoded folders into the specified output folder.
need to provide 2 arguments:

--output_folder','-o', output folder to copy files to
--spec_numbers','-s', csv file where one column contains spec numbers

This will loop through the csv file (comma seperated) containing specimen numbers. This should have a header row, with one column named "Specimen ID" (case sensitive) and will extract all spec numbers from within the rest of the file,
It will then look recursively through two hardcoded folders for CEL files containing the specimen number and copy these to the directory given to --output_folder

If there are 0 or > 1 CEL file this is reported to stdout so the user can select the relevant file (or exclude any with duplicates).


## file\_mover\_prenatal.py
This script is a bit more complex as prenatal samples are only analysed in certain regions. It is run in two parts:

1) take a list of specimen numbers, find the rhchp files and copy them into a subfolder. These can be used to create a msv to be used in step 2:
2) for the list of specimens, if they do not have a call within the a known syndrome region copy them into a subfolder

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
