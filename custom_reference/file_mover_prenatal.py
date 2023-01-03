import argparse
import os
import sys
import shutil
import random
import string
import re
from collections import defaultdict

#TODO make changes/new script to find files for postnatal
def get_args(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_folder','-o',help='output folder to copy files to')
    parser.add_argument('--input_folder','-i',help='folder to search for rhchp files')
    parser.add_argument('--spec_numbers','-s',help='csv file where one column contains spec numbers')
    parser.add_argument('--syndrome_regions','-r',required=False,help='BED file containing syndromic regions')
    parser.add_argument('--syndrome_free_files','-f',required=False,help='folder to contain anonymised cel files that do not have calls overlapping with syndromic regions')
    parser.add_argument('--multi_sample_viewer_output','-m',required=False,help='output of multisample viewer containing calls in multiple samples')
    return parser.parse_args(args)

def check_for_output_dir(args):
    """
    Checks if output folders exist and if not create them
    """
    for folder in [args.output_folder,args.syndrome_free_files]:
        if folder and not os.path.exists(folder):
            os.path.mkdir(folder)
    

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
    Searches (recursively) through the provided input folder for a rhchp file for each speciment number in the list.
    If found it will copy the file to the given output folder
    """
    for spec_number in spec_number_list:
        found=False
        for file in os.listdir(parsed_args.output_folder):
            # files end with .rhchp and can be anywhere within a folder tree
            if re.match(r'(%s).*(.rhchp$)' % (spec_number), file):
                found = True
        
        # if multiple files per specimen id only the first will be taken.
        while not found:
            for root,dirs,files in os.walk(r'%s' % parsed_args.input_folder):
                for file in files:
                    # file ends with .rhchp to exclude some other file types
                    if re.match(r'(%s).*(.rhchp$)' % (spec_number), file):
                        if not os.path.isfile(os.path.join(parsed_args.output_folder,file)):
                            # copy the file into the provided subfolder
                            shutil.copyfile(os.path.join(root,file),os.path.join(parsed_args.output_folder,file))
                        found=True
            if not found:
                print("Error when finding rhchp file - No file found for spec number %s" %(spec_number))
                break

def get_syndrome_regions(parsed_args):
    """
    Open BED file and create a list of tuples
    """
    syndrome_region=[]
    with open(parsed_args.syndrome_regions) as syndrome_BED:
        for region in syndrome_BED.readlines():
            if "start" in region or "stop" in region:
                pass
            else:    
                chr,start,stop,type=region.rstrip().split("\t")
                syndrome_region.append((str(chr.replace("chr","")),int(start),int(stop),type))
    return syndrome_region

def open_multisampleviewer_file(parsed_args):
    """
    multisample viewer file contains all the calls from all samples.
    This can be used to exclude any samples that have a call in a pathogenic region (to ensure these aren't used to normalise samples - would mask real calls)
    create a dictionary with the spec number as a key and list of tuples for value
    """
    #default dict handles creation of keys if already exists
    call_list=defaultdict(list)
    with open(r'%s' % parsed_args.multi_sample_viewer_output) as msvo:
        for line in msvo.readlines():
            if not line.startswith("File"):
                file,type,state,chr,start,stop = line.split("\t")[0:6]
                call_list[file].append((str(chr),int(start),int(stop),type[0:4]))
    return call_list
    
def copy_files_with_no_syndrome_overlaps(parsed_args,multi_sample_viewer_file,syndrome_region_list):
    """
    This function essentially creates a list of samples that can be used to create a new custom reference (AKA panel of normals)
    The files needed to create the reference are .CEL files
    For each sample in the multi sample viewer, check if there is a call which overlaps with a known syndrome.
    If there is such a call that sample should not be used as it's not 'normal'.
    If no calls find the CEL file and copy into folder, but anonymise so we don't know (by process of elimination) which prenatals have a call in a region that they shouldn't be looking in
    """
    skip=0
    not_skip=0
    for sample in multi_sample_viewer_file:
        if not os.path.exists(os.path.join(parsed_args.output_folder,sample)):
            print("rhchp file for sample %s not present in folder %s" % (sample,parsed_args.output_folder))
        else:
            sample_skip=False
            for call in multi_sample_viewer_file[sample]:
                if not sample_skip:
                    # check if call overlaps with any in the syndrome list. If so stop comparing and do not process any further
                    sample_chr,sample_start,sample_stop,sample_type = call
                    for syndrome in syndrome_region_list:
                        synd_chr,synd_start,synd_stop,synd_type = syndrome
                        if synd_chr==sample_chr and sample_start <synd_stop and sample_stop>synd_start and synd_type==sample_type:
                            sample_skip=True
                            skip+=1
            if not sample_skip:
                not_skip+=1
                #create random string to anonymise
                random_file_name = "".join(random.choices(string.ascii_letters + string.digits, k=8))
                cel_file_path = find_cel_file(sample)
                if cel_file_path:
                    shutil.copyfile(cel_file_path,os.path.join(parsed_args.syndrome_free_files,"%s.CEL" % (random_file_name)))
                else:
                    print("CEL file not found for %s" % sample)
    #print("skipped = %s, not skipped =%s" % (skip,not_skip))
    
def find_cel_file(sample_test_number):
    """
    Given a sample_test_number from the multisample viewer (eg 2128184_SNP_220302.1) find the associated cel file
    Note because the sample_test_number is used as opposed to spec number wouldn't expect duplicates.
    return the cel file if found, else return None.
    """
    # if multiple files per specimen id only the first will be taken.
    for root,dirs,files in os.walk(r'S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\Archive\CEL and ARR files do not delete'):
        for file in files:
            #print(files)
            # file ends with .CEL (case sensitive) to exclude some other file types
            if re.match(r'(%s).*(.CEL$)' % (sample_test_number.replace(".rhchp","")), file):
                return os.path.join(root,file)
    
    for root,dirs,files in os.walk(r'S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud'):
        for file in files:
            #print(files)
            # file ends with .CEL (case sensitive) to exclude some other file types
            if re.match(r'(%s).*(.CEL$)' % (sample_test_number.replace(".rhchp","")), file):
                return os.path.join(root,file)
    return None
        
    
def main(args):
    parsed_args=get_args(args)
    check_for_output_dir(parsed_args)
    spec_number_list=create_list_of_spec_numbers(parsed_args)
    find_files(parsed_args,spec_number_list)
    if parsed_args.syndrome_regions:
        syndrome_region_list = get_syndrome_regions(parsed_args)
        multi_sample_viewer_file = open_multisampleviewer_file(parsed_args)
        copy_files_with_no_syndrome_overlaps(parsed_args,multi_sample_viewer_file,syndrome_region_list)

if __name__ == "__main__":
    main(sys.argv[1:])
