import argparse
import os
import sys
import shutil
import re
import datetime
import pandas as pd
from pathlib import Path


def arg_parse():
    """
    Parse arguments supplied by the command line. Create argument parser,
    define command line arguments, then parse supplied command line arguments
    using the created argument parser
        :return (dict): Parsed command line attributes
    """
    parser = argparse.ArgumentParser(
        description=(
            "This script is used to collate CEL files for a list of spec numbers to create a custom reference file"
        )
    )
    parser.add_argument(
        "--output_dir",
        "-o",
        type=validate_path,
        help="Output folder to copy CEL files to",
        required=True,
    )
    parser.add_argument(
        "--spec_number_file",
        "-s",
        type=validate_path,
        help="Path to CSV file where one column contains spec numbers",
        required=True,
    )
    parser.add_argument(
        "--exclude_spec_numbers_file",
        "-e",
        type=validate_path,
        help="Path to CSV file where one column contains spec numbers to exclude",
        required=False,
        default=None,
    )
    # TODO add origin folders argument with list and file validation for each item in list
    #     self.cel_origin_folders = [
    #     r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\Archive\CEL and ARR files do not delete",
    #     r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud",
    # ]
    return vars(parser.parse_args())


def validate_path(path: str) -> str:
    """
    Check path exists
        :param path (str):  Input path to validate
    """
    if os.path.exists(path):
        return path
    else:
        raise argparse.ArgumentTypeError(f"{path} is not a valid path")


class CELMover:
    """
    Collates CEL files for a list of spec numbers to create a custom reference file
    Unlike prenatals this does not require anonymisation or comparison with a BED file

    Attributes
        output_dir (str):                           Directory to save CEL files to
        male_subdir (str):                          Output folder male subdirectory
        female_subdir (str):                        Output folder female subdirectory
        spec_number_file (str):                     Path to CSV file where one column contains spec numbers
        exclude_spec_numbers_file (str):            Path to CSV file where one column contains spec numbers to exclude
        spec_number_list
        (list(specimen_number (str), sex (str))):   List of tuples (specimen number, sex)
        cel_origin_folders (list):                  List of directories that contain CEL files to be searched for
                                                    copying
        filtered_spec_number_dict (dict):           Dictionary of specimen numbers with sex and subdirectories to save
                                                    the files in
        files_to_copy (dict):                       Dictionary of files to copy into the custom reference

    Methods
        create_output_subdirs()
            Checks if output subdirectories exist and if not create them
        create_list_of_spec_numbers()
            Reads the spec number file. From the header, determines which column contains the spec
            numbers ("Specimen ID") and which contains the sex ("Result Type"). If the sex column cannot
            be determined, return "None" in that field
        filter_list_of_spec_numbers()
            Takes the list of specimen numbers read from the spec_numbers input (spec_number_list).
            This is a list of tuples (specimen number, sex). If a file of specimen IDs to exclude was
            given then remove any specimens that appear in this list from spec_number_list. If
            exclude_spec_numbers_file was not given return spec_number_list without any filtering
        find_files()
            Searches (recursively) through two hardcoded folders for a CEL file for each specimen
            number in the list. If found it will copy the file to the given output folder
        remove_spec_nos()
            Where there are multiple files for a spec number, remove this file from the dictionary
            (requested by the Array team)
        copy_files()
            If CEL file does not already exist in the destination directory,
            copy the file from source to destination
    """

    def __init__(
        self,
        output_dir: str,
        spec_number_file: str,
        exclude_spec_numbers_file: str,
        cel_origin_folders: list,
    ):
        """
        Constructor for CELMover class
            :param output_dir (str):                    Directory to save CEL files to
            :param spec_number_file (str):              Path to CSV file where one column contains spec numbers
            :param exclude_spec_numbers_file (str):     Path to CSV file where one column contains
                                                        spec numbers to exclude
            :param cel_origin_folders (list):           List of directories that contain CEL files to be
                                                        searched for copying
        """
        self.output_dir = output_dir
        self.male_subdir = os.path.join(self.output_dir, "male")
        self.female_subdir = os.path.join(self.output_dir, "female")
        self.undetermined_subdir = os.path.join(self.output_dir, "undetermined")
        self.spec_number_file = spec_number_file
        self.exclude_spec_numbers_file = exclude_spec_numbers_file
        self.create_output_subdirs()
        self.spec_number_list = self.create_list_of_spec_numbers()
        self.cel_origin_folders = cel_origin_folders

        if self.exclude_spec_numbers_file:
            self.filtered_spec_number_dict = self.filter_list_of_spec_numbers()
        else:
            print(
                "Spec number filtering not specified (no file "
                "containing spec numbers to exclude was provided)"
            )
        self.files_to_copy = self.find_files()
        self.remove_spec_nos()
        self.copy_files()

    def create_output_subdirs(self) -> None:
        """
        Checks if output subdirectories exist and if not create them
        """
        for directory in [
            self.male_subdir,
            self.female_subdir,
            self.undetermined_subdir,
        ]:
            if not os.path.exists(directory):
                os.mkdir(directory)

    def create_list_of_spec_numbers(self) -> list:
        """
        Reads the spec number file. From the header, determines which column contains the spec
        numbers ("Specimen ID") and which contains the sex ("Result Type"). If the sex column cannot
        be determined, return "None" in that field
            :return spec_number_list (list(specimen_number (str), sex (str))):  List of tuples (specimen number, sex)
        """
        spec_number_list = []
        sex_col = False
        with open(self.spec_number_file, "r") as input_file:
            file_list = input_file.readlines()
            # get the column number containing the specimen number
        for count, header in enumerate(file_list[0].split(",")):
            if header.rstrip() == "Specimen ID":
                spec_column = count
            if header.rstrip() == "Result Type":
                sex_col = count

        for line in file_list[
            1:
        ]:  # Extract the spec number from rest of the lines in file
            if sex_col:
                spec_number_list.append(
                    (line.split(",")[spec_column], line.split(",")[sex_col])
                )
            else:
                spec_number_list.append((line.split(",")[spec_column], "None"))

        print(
            "%s spec numbers in spec number list" % len(spec_number_list)
        )  # Summarise number of specimens
        return spec_number_list

    def filter_list_of_spec_numbers(self) -> list:
        """
        Takes the list of specimen numbers read from the spec_numbers input (spec_number_list).
        This is a list of tuples (specimen number, sex). If a file of specimen IDs to exclude was
        given then remove any specimens that appear in this list from spec_number_list. If
        exclude_spec_numbers_file was not given return spec_number_list without any filtering
            :return filtered_spec_number_dict
            (list(specimen number (int), sex (str))):   List of spec numbers with those marked for exclusion removed
        """
        to_exclude_list, filtered_spec_number_dict, excluded = [], {}, []

        with open(self.exclude_spec_numbers_file, "r") as input_file:
            file_contents = input_file.readlines()
        # Get column number containing the specimen number
        for count, header in enumerate(file_contents[0].split(",")):
            if header == "Specimen ID":
                spec_column = count

        for line in file_contents[
            1:
        ]:  # Extract the specimen number from the rest of the lines in file
            to_exclude_list.append((line.split(",")[spec_column]))

        for specimen in self.spec_number_list:
            spec_no, sex = specimen
            if "normal female" in sex:
                sex_subdir = self.female_subdir
            elif "normal male" in sex:
                sex_subdir = self.male_subdir
            else:
                sex_subdir = "undetermined"
            if spec_no not in to_exclude_list:
                filtered_spec_number_dict[spec_no] = {
                    "sex": sex,
                    "sex_subdir": sex_subdir,
                }
            else:
                excluded.append(specimen)
        assert len(excluded) + len(filtered_spec_number_dict.keys()) == len(
            self.spec_number_list
        )
        # Summarise number of specimens
        print(
            "%s spec numbers in spec number list post filtering"
            % len(filtered_spec_number_dict)
        )
        return filtered_spec_number_dict

    def find_files(self) -> dict:
        """
        Searches (recursively) through hardcoded folders (self.cel_origin_folders) for a CEL file for
        each specimen number in the list. If found it will copy the file to the given output folder
            :return files_to_copy (dict):   Dictionary of files to copy into the custom reference
        """
        self.filtered_spec_number_dict
        files_to_copy = {}

        cel_files = []
        for folder in self.cel_origin_folders:
            cel_files.extend(
                list(Path(folder).rglob("*.CEL$"))
            )  # Extract paths of all CEL files in each folder

        for spec_no in self.filtered_spec_number_dict.keys():
            print("Spec number %s" % spec_no)
            files_to_copy[spec_no] = {}
            for cel_file in cel_files:
                if re.match(r"(%s).*(.CEL$)" % (spec_no), cel_file):
                    files_to_copy[spec_no][cel_file] = (
                        {  # Preliminary list of file names to copy for that spec no
                            "src": cel_file,
                            "dest": os.path.join(
                                self.filtered_spec_number_dict[spec_no]["sex_subdir"],
                                cel_file,
                            ),
                        }
                    )
                    cel_files.remove(
                        cel_file
                    )  # Remove from list to reduce the search burden
                    print(
                        "Spec no %s. Added file %s to files_to_copy" % (spec_no, file)
                    )
        print(files_to_copy)
        return files_to_copy

    def remove_spec_nos(self) -> None:
        """
        Where there are multiple files for a spec number, remove this file
        from the dictionary (requested by the Array team)
        """
        for spec_no in self.files_to_copy.keys():
            if self.files_to_copy[spec_no].keys() > 1:
                self.files_to_copy.pop(spec_no)
                print(
                    f"Multiple files identified for spec number {spec_no} ({self.files_to_copy[spec_no].keys()}). "
                    "These files have been excluded from files to copy"
                )

    def copy_files(self) -> None:
        """
        If CEL file does not already exist in the destination directory,
        copy the file from source to destination
        """
        for spec_no in self.files_to_copy.keys():
            if not os.path.isfile(self.files_to_copy[spec_no]["dest"]):
                print(
                    f"Copying CEL file. Src: {self.files_to_copy[spec_no]['src']}. "
                    f"Dest: {self.files_to_copy[spec_no]['dest']}"
                )
                shutil.copyfile(  # Copy the file into the provided subfolder
                    self.files_to_copy[spec_no]["src"],
                    self.files_to_copy[spec_no]["dest"],
                )


if __name__ == "__main__":

    parsed_args = arg_parse()

    CELMover(
        parsed_args["output_dir"],
        parsed_args["spec_number_file"],
        parsed_args["exclude_spec_numbers_file"],
    )
