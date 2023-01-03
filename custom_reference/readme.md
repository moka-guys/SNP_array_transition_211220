# filemover.py
this script can be used to create a custom reference file (which is used a bit like a panel of normals) to normalise the data and remove noise

There are 2 use cases
1) take a list of specimen numbers, find the rhchp files and copy them into a subfolder. These can be used to create a msv to be used in step 2:
2) for the list of specimens, if they do not have a call within the a known syndrome region copy them into a subfolder

# step 1 - Collect rhchp files:
need to provide 3 arguments:

--output_folder','-o', output folder to copy files to
--input_folder','-i', folder to search for rhchp files
--spec_numbers','-s', csv file where one column contains spec numbers

This will loop through the csv file (comma seperated) containing specimen numbers. This should have a header row, with one column named "Specimen ID" (case sensitive) and will extract all spec numbers from within the rest of the file,
It will then look recursively through the --input_folder path and copy to the --output_folder directory

# step 2 - identify files that can be used to make a custom reference
repeats step 1 but requires 3 additional arguments:

--syndrome_regions,'-r', optional, BED file containing syndromic regions
--syndrome_free_files, '-f', optional , folder to contain anonymised rhchp files that do not have calls overlapping with syndromic regions
--multi_sample_viewer_output','-m',optional, output of multisample viewer containing calls in multiple samples

With the additional 3 arguments the script will then assess all the calls from the multisample viewer file. 
For each sample, if there are no calls which overlap the syndromes provided in the syndrome BED file (tab seperated in format chr, start, stop , type (Gain/Loss)) then it can be used to make the custom reference. 
To make a custom reference .CEL files are required. The script looks in (hardcoded) directories for the file, it will anonymise the file and copy into the directory provided to the --syndrome_free_files argument.
