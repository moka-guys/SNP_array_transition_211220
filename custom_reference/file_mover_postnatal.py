import argparse
import os
import sys
import shutil
import re
import datetime
import pandas as pd
from pathlib import Path
from logger import Logger
import tkinter as tk
from functools import partial
from tkinter import messagebox


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
        spec_number_dict (dict):					Dictionary of spec nos containing their sex and
                                                    sex subdir for the file to be saved in
        cel_origin_folders (list):                  List of directories that contain CEL files to be searched for copying
        filtered_spec_number_dict (dict):  			spec_number_dict with specimens marked for exclusion removed
        cel_files (list):                           List of all CEL files in supplied folder locations and their subdirectories
        files_to_copy (dict):                       Dictionary of files to copy into the custom reference

    Methods
        create_output_subdirs()
            Checks if output subdirectories exist and if not create them
        create_spec_number_dict()
            Reads the spec number file and creates dictionary of spec numbers containing sex, and sex subdir
            to save the CEL file to. These strings are obtained from the file using the header to extract
            data from the correct columns (spec number = "Specimen ID", sex = "Result Type". If the sex
            column cannot be determined, sex entry is defined as "None"
        filter_spec_number_dict()
            Takes the specimen number dict and removes any entries that are specified in the
            exclude_spec_numbers_file. If the exclude_spec_numbers_file was not supplied,
            the unchanged dictionary is returned
        find_files()
            Searches (recursively) through hardcoded folders (self.cel_origin_folders) to identify all CEL files and
            return these as a list
        identify_files_to_copy()
            Generate a dictionary of CEL files to be copied into the custom reference by identifying CEL files that
            match each specimen number in self.filtered_spec_number_dict, identifying which of those should be copied
            and adding those to the self.files_to_copy dict using self.add_to_files_to_copy_dict()
        add_to_files_to_copy_dict(spec_no, run_no, cel_file_str)
            Add files to self.files_to_copy dictionary
        remove_spec_nos()
            Where there are multiple files for a spec number, remove this spec number from the dictionary. Keeps files
            that are identical but found in multiple folder locations, but keeps only a single version of the file
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
        self.spec_number_dict = self.create_spec_number_dict()
        self.cel_origin_folders = cel_origin_folders

        if self.exclude_spec_numbers_file:
            logger.info(
                "A CSV file containing spec numbers for inclusion has been provided"
            )
            self.filtered_spec_number_dict = self.filter_spec_number_dict()
        else:
            logger.warning(
                "Spec number filtering not specified (no file "
                "containing spec numbers to exclude was provided)"
            )
        self.cel_files = self.find_files()
        self.files_to_copy = {}
        self.identify_files_to_copy()
        self.remove_spec_nos()
        self.copy_files()
        logger.info("Script complete")

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

    def create_spec_number_dict(self) -> list:
        """
        Reads the spec number file and creates dictionary of spec numbers containing sex, and sex subdir
                to save the CEL file to. These strings are obtained from the file using the header to extract
        data from the correct columns (spec number = "Specimen ID", sex = "Result Type". If the sex
                column cannot be determined, sex entry is defined as "None"
            :return spec_number_dict (dict):	Dictionary of spec nos containing their sex and
                                                                                                sex subdir for the file to be saved in
        """
        spec_number_dict = {}
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
            spec_no = line.split(",")[spec_column]
            if sex_col:
                spec_number_dict[spec_no] = {"sex": line.split(",")[sex_col]}
            else:
                spec_number_dict[spec_no] = {
                    "sex": "None",
                }
            if "normal female" in spec_number_dict[spec_no]["sex"]:
                sex_subdir = self.female_subdir
            elif "normal male" in spec_number_dict[spec_no]["sex"]:
                sex_subdir = self.male_subdir
            else:
                sex_subdir = "undetermined"

            spec_number_dict[spec_no]["sex_subdir"] = sex_subdir
        logger.info(
            "%s spec numbers in spec number dictionary" % len(spec_number_dict.keys())
        )  # Summarise number of specimens
        return spec_number_dict

    def filter_spec_number_dict(self) -> list:
        """
        Takes the specimen number dict and removes any entries that are specified in the
        exclude_spec_numbers_file. If the exclude_spec_numbers_file was not supplied,
        the unchanged dictionary is returned
            :return filtered_spec_number_dict (dict):   spec_number_dict with specimens marked for exclusion removed
        """
        to_exclude_list, excluded = [], []
        filtered_spec_number_dict = self.spec_number_dict.copy()

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

        list(set(to_exclude_list))

        for spec_no in self.spec_number_dict.keys():
            if spec_no in to_exclude_list:
                del filtered_spec_number_dict[spec_no]
                excluded.append(spec_no)

        logger.info(f"{len(excluded)} specimen numbers excluded")
        assert len(excluded) + len(filtered_spec_number_dict.keys()) == len(
            self.spec_number_dict
        )
        # Summarise number of specimens
        logger.info(
            "%s spec numbers in spec number dictionary post filtering"
            % len(filtered_spec_number_dict)
        )
        return filtered_spec_number_dict

    def find_files(self) -> dict:
        """
        Searches (recursively) through hardcoded folders (self.cel_origin_folders) to identify all CEL files and
        return these as a list
            :return cel_files (list):   List of all CEL files in supplied folder locations and their subdirectories
        """
        self.filtered_spec_number_dict

        cel_files = []
        for folder in self.cel_origin_folders:
            if os.path.isdir(folder):
                logger.info(f"The input directory exists: {folder}")
                cel_files.extend(
                    list(Path(folder).rglob("*.CEL"))
                )  # Extract paths of all CEL files in each folder
            else:
                logger.warning(
                    f"The input directory does not exist, skipping: {folder}"
                )
        logger.info("%s CEL files identified in input folders" % len(cel_files))
        return cel_files

    def identify_files_to_copy(self) -> None:
        """
        Generate a dictionary of CEL files to be copied into the custom reference by identifying CEL files that
        match each specimen number in self.filtered_spec_number_dict, identifying which of those should be copied
        and adding those to the self.files_to_copy dict using self.add_to_files_to_copy_dict()
        """
        file_count_to_copy = 0
        for spec_no in self.filtered_spec_number_dict.keys():
            self.files_to_copy[spec_no] = {}
            for cel_file in self.cel_files:
                cel_file_str = str(cel_file)
                if str(spec_no) in cel_file_str:
                    run_no = (cel_file_str.split("\\")[-1]).split("_")[1]
                    # If there is already an entry for this run number
                    if run_no in self.files_to_copy[spec_no].keys():
                        logger.info(
                            "Spec number with the same run number as an another CEL file has been identified. "
                            f"Spec no: {spec_no}. Run no: {run_no}"
                        )
                        # If the existing spec_no dictionary entry is a CEL file with an identical run name,
                        # skip over the CEL file as it is likely the same file stored in a different location
                        run_names = [
                            k.split("\\")[-1]
                            for k in self.files_to_copy[spec_no].keys()
                        ]  # Extract run names
                        if any(
                            cel_file_str.split("\\")[-1] in run_name
                            for run_name in run_names
                        ):
                            logger.info(
                                "CEL file already exists in the dictionary but is stored in multiple "
                                f"locations. Skipped. ({cel_file_str})"
                            )
                        else:  # Spec no. and run no. not in dict yet
                            self.add_to_files_to_copy_dict(
                                spec_no, run_no, cel_file_str
                            )
                            file_count_to_copy += 1
                    else:  # Spec no. and run no. not in dict yet
                        self.add_to_files_to_copy_dict(spec_no, run_no, cel_file_str)
                        file_count_to_copy += 1
                    self.cel_files.remove(
                        cel_file
                    )  # Remove from list to reduce the search burden
        logger.info(f"Identified {file_count_to_copy} files matching specimen numbers")

    def add_to_files_to_copy_dict(
        self, spec_no: str, run_no: str, cel_file_str: str
    ) -> None:
        """
        Add files to self.files_to_copy dictionary
            :param spec_no (str):       Specimen number
            :param run_no (str):        Run number
            :param cel_file_str (str):  CEL file path
        """
        self.files_to_copy[spec_no][cel_file_str] = (
            {  # Preliminary list of file names to copy for that spec no
                "run_no": run_no,
                "src": cel_file_str,
                "dest": os.path.join(
                    self.filtered_spec_number_dict[spec_no]["sex_subdir"],
                    cel_file_str.split("\\")[-1],
                ),
            }
        )

    def remove_spec_nos(self) -> None:
        """
        Where there are multiple files for a spec number, remove this spec number from the dictionary. Keeps files
        that are identical but found in multiple folder locations, but keeps only a single version of the file
        """
        final_dict = self.files_to_copy.copy()
        for spec_no in self.files_to_copy.keys():
            if len(self.files_to_copy[spec_no].keys()) > 1:
                files_to_exclude = ", ".join(self.files_to_copy[spec_no].keys())
                final_dict.pop(spec_no)
                logger.info(
                    f"Multiple files identified for spec number {spec_no} ({files_to_exclude}). "
                    "These files have been excluded from files to copy"
                )
            if len(self.files_to_copy[spec_no].keys()) == 0:
                logger.info(f"No files were identified for spec number {spec_no}.")
                final_dict.pop(spec_no)
        self.files_to_copy = final_dict

    def copy_files(self) -> None:
        """
        If CEL file does not already exist in the destination directory,
        copy the file from source to destination
        """
        for spec_no in self.files_to_copy.keys():
            if len(self.files_to_copy[spec_no].keys()) == 1:
                for file_name in self.files_to_copy[spec_no].keys():
                    if not os.path.isfile(
                        self.files_to_copy[spec_no][file_name]["dest"]
                    ):
                        logger.info(
                            f"Copying CEL file. Src: {self.files_to_copy[spec_no][file_name]['src']}. "
                            f"Dest: {self.files_to_copy[spec_no][file_name]['dest']}"
                        )
                    shutil.copyfile(  # Copy the file into the provided subfolder
                        self.files_to_copy[spec_no][file_name]["src"],
                        self.files_to_copy[spec_no][file_name]["dest"],
                    )


class MessageBox:
    """
    Class for creating a message box for user input

    Attributes
        root (tk.Tk):           Main tkinter window
        var (tk.StringVar):     Tkinter variable to store selected choice
        choice (str | None):    Variable to store selected choice

    Methods
        setup_window()

        submit()
    """

    def __init__(self):
        """
        Constructor for the MessageBox class
        """
        self.root = tk.Tk()  # Create the main tkinter window
        self.var = (
            tk.StringVar()
        )  # Create a tkinter variable to store the selected choice
        self.choice = None  # Variable to store selected choice

    def setup_window(self) -> None:
        """
        Configure the tkinter window
        """
        self.root.title(
            "Please select where you are running the script from"
        )  # Set the title of the window
        label = tk.Label(self.root, text="Select an option:")  # Create a label widget
        label.pack()  # Display the label in the window
        option1 = tk.Radiobutton(
            self.root, text="VM", variable=self.var, value="VM"
        )  # Create the first radio button
        option1.pack()  # Display the first radio button in the window
        option2 = tk.Radiobutton(
            self.root, text="S Drive", variable=self.var, value="S Drive"
        )  # Create the second radio button
        option2.pack()  # Display the second radio button in the window
        submit_button = tk.Button(
            self.root, text="Submit", command=self.submit
        )  # Create the Submit button
        submit_button.pack()  # Display the Submit button in the window
        self.root.mainloop()  # Start the tkinter event loop

    def submit(self) -> None:
        """
        Runs upon submission of message box submit button. Displays the selected choice
        in a message box, saves the choice, and closes the window
        """
        self.choice = self.var.get()
        # You can save the choice as a variable or perform any other action here
        messagebox.showinfo("Selection", f"You selected: {self.choice}")
        logger.info(f"Choice saved: {self.choice}")
        self.root.destroy()


if __name__ == "__main__":
    parsed_args = arg_parse()
    logfile_path = os.path.join(
        os.getcwd(),
        f"{datetime.datetime.now():%y%m%d_%H%M%S}_file_mover_postnatal.log",
    )
    logger = Logger(logfile_path).logger
    logger.info("Running file_mover_postnatal.py")
    message_box = MessageBox()
    message_box.setup_window()

    if message_box.choice == "S Drive":
        cel_origin_folders = [
            r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\Archive\CEL and ARR files do not delete",
            r"S:\Genetics_Data2\Array\Geneworks - Viapath Cloud sync folder\UploadToCloud",
        ]

    elif message_box.choice == "VM":
        cel_origin_folders = [
            r"\\GRPVCHASDB01\Archive\CEL and ARR files do not delete",
            r"\\GRPVCHASDB01\Genetics\In",
        ]

    logger.info(f"CEL origin folders set as: {', '.join(cel_origin_folders)}")
    CELMover(
        parsed_args["output_dir"],
        parsed_args["spec_number_file"],
        parsed_args["exclude_spec_numbers_file"],
        cel_origin_folders,
    )
