#
# This script is used to collate CEL files for a list of spec numbers to create a custom reference file
# Unlike prenatals this does not require anonymisation or comparison with a BED file.
#
import argparse
import os
import sys
import shutil
import re


def get_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder','-o',help='output folder to copy files to')
    parser.add_argument('--spec_numbers','-s',help='csv file where one column contains spec numbers')
    return parser.parse_args(args)

def check_for_output_dir(args):
    """
    Checks if output folders exist and if not create them
    """
    for folder in [args.output_folder]:
        if folder and not os.path.exists(folder):
            os.mkdir(folder)
    

def create_list_of_spec_numbers(args):
    """
    Reads the spec number file provided as an argument
    From the header determine which column contains spec numbers (titled "Specimen ID") 
    Returns a list of specimen numbers
    """
    spec_number_list=[]
    with open(r'%s' % args.spec_numbers) as input_file:
        file_list=input_file.readlines()
        # get the column number containing the specimen number
    for count,header in enumerate(file_list[0].split(",")):
        if header=="Specimen ID":
            spec_column=count
    # extract the spec number from rest of the lines in file
    for line in file_list[1:]:
        spec_number_list.append(line.split(",")[spec_column])
    #summarise number of specimens
    print("%s in spec number list" % len(spec_number_list)) 
    return spec_number_list

def find_files(parsed_args,spec_number_list):
    """
    Searches (recursively) through two hardcoded folders for a CEL file for each speciment number in the list.
    If found it will copy the file to the given output folder
    If there are 0 or more than 1 CEL file found it will report the count for that spec number
    """
    for spec_number in spec_number_list:
        count=0
        folders = [r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\Archive\CEL and ARR files do not delete", r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud"]
        for folder in folders:
            for root,dirs,files in os.walk(r'%s' % folder):
                for file in files:
                    # file ends with .CEL to exclude some other file types
                    if re.match(r'(%s).*(.CEL$)' % (spec_number), file):
                        if not os.path.isfile(os.path.join(parsed_args.output_folder,file)):
                            # copy the file into the provided subfolder
                            shutil.copyfile(os.path.join(root,file),os.path.join(parsed_args.output_folder,file))
                        count+=1
            if count != 1:
                print("Warning - %s CEL files found for spec number %s" %(count, spec_number))

        
    
def main(args):
    parsed_args=get_args(args)
    check_for_output_dir(parsed_args)
    spec_number_list=create_list_of_spec_numbers(parsed_args)
    find_files(parsed_args,spec_number_list)


if __name__ == "__main__":
    main(sys.argv[1:])