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
    parser.add_argument('--output_folder','-o',help='output folder to copy files to', required=True)
    parser.add_argument('--spec_numbers','-s',help='csv file where one column contains spec numbers', required=True)
    parser.add_argument('--exclude_spec_numbers','-e',help='csv file where one column contains spec numbers to exclude', required=False)
    return parser.parse_args(args)

def check_for_output_dir(args): 
    """
    Checks if output folders exist and if not create them
    """
    for folder in [args.output_folder]:
        male_folder=os.path.join(folder,"male")
        female_folder=os.path.join(folder,"female")
        if male_folder and not os.path.exists(male_folder):
            os.mkdir(male_folder)
        if female_folder and not os.path.exists(female_folder):
            os.mkdir(female_folder)

def create_list_of_spec_numbers(args):
    """
    Reads the spec number file provided as an argument
    From the header determine which column contains spec numbers (titled "Specimen ID") and which contains the sex (Result Type)
    If sex column can't be determined return string "None" in that field
    Returns a list of tuples (specimen number, sex)
    """
    spec_number_list=[]
    sex_col=False
    with open(r'%s' % args.spec_numbers) as input_file:
        file_list=input_file.readlines()
        # get the column number containing the specimen number
    for count,header in enumerate(file_list[0].split(",")):
        if header.rstrip()=="Specimen ID":
            spec_column=count
        if header.rstrip()=="Result Type":
            sex_col=count

    # extract the spec number from rest of the lines in file
    for line in file_list[1:]:
        if sex_col:
            spec_number_list.append((line.split(",")[spec_column],line.split(",")[sex_col]))
        else:
            spec_number_list.append((line.split(",")[spec_column],"None"))
    #summarise number of specimens
    print("%s in spec number list" % len(spec_number_list)) 
    return spec_number_list


def filter_list_of_spec_numbers(args,spec_number_list):
    """
    Takes the list of specimen numbers read from the spec_numbers input (spec_number_list). This is a list of tuples (specimen number, sex)
    If a file of specimen IDs to exclude was given then remove any specimens that appear in this list from spec_number_list
    If exclude_spec_numbers was not given return spec_number_list without any filtering
    :param list spec_number_list: a list of tuples (specimen number, sex)
    """
    if not args.exclude_spec_numbers:
        #summarise number of specimens
        print("%s in spec number list post filtering" % len(spec_number_list)) 
        return spec_number_list
    else:
        to_exclude_list = []
        filtered_list = []
        excluded = []
        with open(r'%s' % args.exclude_spec_numbers) as input_file:
            file_list=input_file.readlines()
        # get the column number containing the specimen number
        for count,header in enumerate(file_list[0].split(",")):
            if header=="Specimen ID":
                spec_column=count

        # extract the spec number from rest of the lines in file
        for line in file_list[1:]:
            to_exclude_list.append((line.split(",")[spec_column]))
        
        for specimen in spec_number_list:
            spec_no, sex = specimen
            if spec_no not in to_exclude_list:
                filtered_list.append(specimen)
            else:
                excluded.append(specimen)
        assert len(excluded)+len(filtered_list)==len(spec_number_list)
        
        #summarise number of specimens
        print("%s in spec number list post filtering" % len(filtered_list)) 
        return filtered_list

def find_files(parsed_args,spec_number_list):
    """
    Searches (recursively) through two hardcoded folders for a CEL file for each speciment number in the list.
    If found it will copy the file to the given output folder
    If there are 0 or more than 1 CEL file found it will report the count for that spec number
    :param list spec_number_list: a list of tuples (specimen number, sex)
    """
    for spec_number in spec_number_list:
        spec_no, sex = spec_number
        sex_dir=False
        if "normal female" in sex:
            sex_dir = "female"
        elif "normal male" in sex:
            sex_dir = "male"
        else:
            sex_dir = "undetermined"
        count=0
        folders = [r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\Archive\CEL and ARR files do not delete", r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud",r"\\GRPVCHASDB01\Archive\CEL and ARR files do not delete"]
        output_dir = os.path.join(parsed_args.output_folder,sex_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        for folder in folders:
            if os.path.isdir(folder):
                for root,dirs,files in os.walk(r'%s' % folder):
                    for file in files:
                        # file ends with .CEL to exclude some other file types
                        if re.match(r'(%s).*(.CEL$)' % (spec_no), file):
                            if not os.path.isfile(os.path.join(parsed_args.output_folder,file)):
                                # copy the file into the provided subfolder
                                shutil.copyfile(os.path.join(root,file),os.path.join(parsed_args.output_folder,sex_dir,file))
                                count+=1
            else:
                print("%s is not a folder" % folder)
        if count != 1:
            print("Warning - %s CEL files found for spec number %s" %(count, spec_number))

        
    
def main(args):
    parsed_args=get_args(args)
    check_for_output_dir(parsed_args)
    spec_number_list=create_list_of_spec_numbers(parsed_args)
    filtered_spec_number_list = filter_list_of_spec_numbers (parsed_args,spec_number_list)
    find_files(parsed_args,filtered_spec_number_list)


if __name__ == "__main__":
    main(sys.argv[1:])
