import argparse
import os
import sys
import shutil
import re
import datetime
import pandas as pd
from pathlib import Path
import logging
from logger import Logger
import tkinter as tk
from functools import partial
from tkinter import messagebox, font
from concurrent.futures import ProcessPoolExecutor


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
        spec_number_file (str):                     Path to CSV file where one column contains spec numbers
        exclude_spec_numbers_file (str):            Path to CSV file where one column contains spec numbers to exclude
        cel_origin_folders (list):                  List of directories that contain CEL files to be
                                                    searched for copying
        spec_col_name (str):                        Specimen number column header
        sex_col_name (str):                         Sex column header
        sex_subdirs (dict):                         Name of subdirectories to place CEL files in
        files_to_copy_dict (dict):                  Dictionary of files to copy into the custom reference
        specimen_number_dict (dict):			    Dictionary of spec nos containing their sex and sex subdir for the
                                                    file to be saved in, with specimen numbers in the
                                                    exclude_spec_numbers_file filtered out
        cel_files (list):                           List of all CEL files in supplied folder locations and their
                                                    subdirectories which contain desired specimen numbers

    Methods
        get_specimen_number_dict()
            Return filtered specimen number dictionary if exclude_spec_numbers_file has been supplied, else
            return unfiltered specimen number dictionary
        create_specimen_number_dict()
            Reads the spec number file and creates dictionary of spec numbers containing sex, and the sex subdir
            to save the CEL file to. These strings are obtained from the file using the header to extract
            data from the correct columns (spec number = "Specimen ID", sex = "Result Type". If the sex column
            cannot be determined, sex entry is defined as "None"
        get_csv_lines(file)
            Extract lines from CSV file and return as list of strings
        get_column_index(col_name)
            Return column index for input column name
        filter_specimen_numbers(specimen_number_dict)
            Takes the specimen number dict and removes spec numbers specified in the exclude_spec_numbers_file
        find_cel_files()
            Searches (recursively) through hardcoded folders (self.cel_origin_folders) to identify all CEL files and
            return these as a list
        add_cel_files_to_dict()
            Identify CEL files that match a spec number in the self.specimen_number_dict and call
            self.add_file_to_dict() to add the CEL file to the dictionary. If file is not a duplicate of a file that
            has already been added to sthe dictionary, add file to dictionary
        add_file_to_dict(spec_no, cel_file, file_name)
            Add file to dictionary of files to copy
        remove_spec_nos_in_multiple_runs()
            Remove CEL files from the dictionary that are for specimen numbers that appear in multiple runs
        create_output_subdirs()
            Checks if output subdirectories exist and if not create them
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
        self.spec_number_file = spec_number_file
        self.exclude_spec_numbers_file = exclude_spec_numbers_file
        self.cel_origin_folders = cel_origin_folders
        self.spec_col_name = "Specimen ID"
        self.sex_col_name = "Result Type"
        self.sex_subdirs = {
            "normal female": os.path.join(self.output_dir, "female"),
            "normal male": os.path.join(self.output_dir, "male"),
            "undetermined": os.path.join(self.output_dir, "undetermined"),
        }
        self.files_to_copy_dict = {}
        self.specimen_number_dict = self.get_specimen_number_dict()
        self.cel_files = self.find_cel_files()
        if self.cel_files:
            self.add_cel_files_to_dict()
            self.remove_spec_nos_in_multiple_runs()
            self.create_output_subdirs()
            self.copy_files()
        logger.info("Script complete")

    def get_specimen_number_dict(self) -> dict:
        """
        Return filtered specimen number dictionary if exclude_spec_numbers_file has been supplied, else
        return unfiltered specimen number dictionary
            :return specimen_number_dict (dict):   Dictionary of specimen numbers
        """
        specimen_number_dict = self.create_specimen_number_dict()
        if self.exclude_spec_numbers_file:
            logger.info(
                "A CSV file containing spec numbers for inclusion has been provided"
            )
            specimen_number_dict = self.filter_specimen_numbers(specimen_number_dict)
        else:
            logger.warning(
                "Spec number filtering not specified (no file "
                "containing spec numbers to exclude was provided)"
            )
        return specimen_number_dict

    def create_specimen_number_dict(self) -> dict:
        """
        Reads the spec number file and creates dictionary of spec numbers containing sex, and the sex subdir
        to save the CEL file to. These strings are obtained from the file using the header to extract
        data from the correct columns (spec number = "Specimen ID", sex = "Result Type". If the sex column
        cannot be determined, sex entry is defined as "None"
            :return specimen_number_dict (dict):	Dictionary of spec nos containing their sex and
                                                    sex subdir for the file to be saved in
        """
        specimen_number_dict = {}
        file_lines = self.get_csv_lines(self.spec_number_file)
        for line in file_lines[
            1:
        ]:  # Extract the spec number from rest of the lines in file
            spec_no = line.split(",")[self.get_column_index(self.spec_col_name)]
            sex = line.split(",")[self.get_column_index(self.sex_col_name)]
            for subdir in self.sex_subdirs.keys():
                if subdir in sex:
                    specimen_number_dict[spec_no] = {
                        "sex": sex,
                        "sex_subdir": self.sex_subdirs[subdir],
                    }
        logger.info(
            f"{len(specimen_number_dict.keys())} spec numbers in spec number dictionary"
        )  # Summarise number of specimens
        return specimen_number_dict

    def get_csv_lines(self, file: str) -> list:
        """
        Extract lines from CSV file and return as list of strings
            :param file (str):          Path to file
            :return file_list (list):   List of strings, one per line of the file
        """
        with open(file, "r") as input_file:
            file_list = input_file.readlines()
        return file_list

    def get_column_index(self, col_name: str):
        """
        Return column index for input column name
            :param col_name (str):      Column name string
            :return col_index (int):    Index of column
        """
        file_lines = self.get_csv_lines(self.spec_number_file)
        for count, header in enumerate(file_lines[0].split(",")):
            if header.rstrip() == col_name:
                col_index = count
        return col_index

    def filter_specimen_numbers(self, specimen_number_dict: dict):
        """
        Takes the specimen number dict and removes spec numbers specified in the exclude_spec_numbers_file
            :param specimen_number_dict (dict):	            Dictionary of spec nos containing their sex and
                                                            sex subdir for the file to be saved inmber_dict
            :return filtered_specimen_number_dict (dict):   specimen_number_dict with specimens marked for
                                                            exclusion removed
        """
        to_exclude_list = []
        excluded = 0
        filtered_specimen_number_dict = specimen_number_dict.copy()
        file_lines = self.get_csv_lines(self.exclude_spec_numbers_file)
        spec_col_index = self.get_column_index(self.spec_col_name)

        # Extract the specimen number from the rest of the lines in file
        for line in file_lines[1:]:
            to_exclude_list.append((line.split(",")[spec_col_index]))

        to_exclude_list = list(set(to_exclude_list))
        logger.info(f"{len(to_exclude_list)} specimen numbers in to exclude file")

        for spec_no in specimen_number_dict.keys():
            if spec_no in to_exclude_list:
                del filtered_specimen_number_dict[spec_no]
                excluded += 1

        logger.info(f"{excluded} specimen numbers excluded")
        assert excluded + len(filtered_specimen_number_dict.keys()) == len(
            specimen_number_dict
        )
        # Summarise number of specimens
        logger.info(
            f"{len(filtered_specimen_number_dict)} spec numbers in spec number dictionary post filtering"
        )
        return filtered_specimen_number_dict

    def find_cel_files(self) -> list:
        """
        Searches (recursively) through hardcoded folders (self.cel_origin_folders) to identify all CEL files and
        return these as a list
            :return cel_files (list):   List of all CEL files in supplied folder locations and their subdirectories
        """
        cel_files = []
        validated_folders = []
        for folder in self.cel_origin_folders:
            if os.path.exists(folder):
                logger.info(f"The input directory exists: {folder}")
                folder_cel_files = list(Path(folder).rglob("*.CEL"))
                logger.info(
                    f"Identified {len(folder_cel_files)} CEL files in input folder: {folder}"
                )
                cel_files.extend(
                    folder_cel_files
                )  # Extract paths of all CEL files in each folder
                validated_folders.append(folder)
            else:
                logger.warning(
                    f"The input directory does not exist, skipping: {folder}"
                )
        self.cel_origin_folders = validated_folders
        logger.info(f"{len(cel_files)} total CEL files identified in input folders")
        return cel_files

    def add_cel_files_to_dict(self) -> None:  # Docstring done
        """
        Identify CEL files that match a spec number in the self.specimen_number_dict and call self.add_file_to_dict()
        to add the CEL file to the dictionary. If file is not a duplicate of a file that has already been added to
        the dictionary, add file to dictionary
        """
        logger.info(
            f"Identifying CEL files that contain specimen numbers in the filtered specimen number list"
        )
        match_no_spec_nos = 0
        duplicate_file_count = 0
        spec_nos = list(self.specimen_number_dict.copy().keys())
        for cel_file in self.cel_files:
            cel_file_str = cel_file.__str__()
            if any(spec_no in cel_file_str for spec_no in spec_nos):
                for spec_no in spec_nos:
                    file_name = cel_file_str.split("\\")[-1]
                    if spec_no in cel_file_str:
                        if file_name in self.files_to_copy_dict:
                            duplicate_file_count += 1
                        else:  # File is not a duplicate, so add to the dict
                            self.add_file_to_dict(spec_no, cel_file, file_name)
            else:
                match_no_spec_nos += 1
        logger.info(
            f"{duplicate_file_count} files were discounted as they were determined to be likely duplicates"
        )
        logger.info(
            f"{match_no_spec_nos} files were discounted as they do not match input spec numbers"
        )
        logger.info(  # TODO update this with the correct way of counting this!
            f"{len(self.files_to_copy_dict)} CEL files identified for copying that match specimen numbers in the filtered specimen number list"
        )

    def add_file_to_dict(self, spec_no: str, cel_file: Path, file_name: str) -> None:
        """
        Add file to dictionary of files to copy
            :param spec_no (str):			Specimen number for the CEL file
            :param cel_file (Path):         Specimen number for the CEL file
            :param file_name (str):         Name of CEL file
        """
        duplicate_file = False
        cel_file_str = cel_file.__str__()
        self.files_to_copy_dict[file_name] = {
            "src": cel_file,
            "dest": os.path.join(
                self.specimen_number_dict[spec_no]["sex_subdir"],
                cel_file.__str__().split("\\")[-1],
            ),
            "spec_no": (cel_file_str.split("\\")[-1]).split("_")[0],
            "run_no": (cel_file_str.split("\\")[-1]).split("_")[1],
        }

    def remove_spec_nos_in_multiple_runs(self) -> dict:
        """
        Remove CEL files from the dictionary that are for specimen numbers that appear in multiple runs
        """
        spec_nos = list(self.specimen_number_dict.copy().keys())
        cel_files_removed = 0
        for spec_no in spec_nos:
            duplicate_spec_no_files = [
                file_name
                for file_name in self.files_to_copy_dict.keys()
                if self.files_to_copy_dict[file_name]["spec_no"] == spec_no
            ]
            if len(duplicate_spec_no_files) > 1:
                # Remove CEL files from dictionary for this spec number as they are from multiple runs
                for file_name in duplicate_spec_no_files:
                    self.files_to_copy_dict.pop(file_name)
                    cel_files_removed += 1
        logger.info(
            f"{cel_files_removed} CEL files removed as the specimen number was present in CEL files from multiple runs"
        )
        logger.info(f"{len(self.files_to_copy_dict)} files remain for copying")

    def create_output_subdirs(self) -> None:
        """
        Checks if output subdirectories exist and if not create them
        """
        for dirtype, dirpath in self.sex_subdirs.items():
            if not os.path.exists(dirpath):
                os.mkdir(dirpath)
                logger.info(f"Created directory: {dirpath}")

    def copy_files(self) -> None:
        """
        If CEL file does not already exist in the destination directory,
        copy the file from source to destination
        """
        logger.info(
            f"{len(self.files_to_copy_dict)} total CEL files to copy. Starting."
        )
        for file_name in self.files_to_copy_dict:
            if not os.path.isfile(self.files_to_copy_dict[file_name]["dest"]):
                shutil.copyfile(  # Copy the file into the provided subfolder
                    self.files_to_copy_dict[file_name]["src"],
                    self.files_to_copy_dict[file_name]["dest"],
                )
        logger.info(f"CEL file copying complete")


class MessageBox:
    """
    Class for creating a message box for user input

    Attributes
        root (tk.Tk):           Main tkinter window
        var (tk.StringVar):     Tkinter variable to store selected choice
        choice (str | None):    Variable to store selected choice

    Methods
        setup_window()
            Configure the tkinter window
        on_closing()
            Write to log and exit script upon closing the message box
        submit()
            Runs upon submission of message box submit button. Displays the selected choice
            in a message box, saves the choice, and closes the window
    """

    def __init__(self, logger):
        """
        Constructor for the MessageBox class
        """
        self.logger = logger
        self.root = tk.Tk()  # Create the main tkinter window
        self.var = (
            tk.StringVar()
        )  # Create a tkinter variable to store the selected choice
        self.choice = None  # Variable to store selected choice

    def setup_window(self) -> None:
        """
        Configure the tkinter window
        """
        helv = font.Font(family="Helvetica", size=20, weight=font.BOLD)
        self.root.geometry("520x300")
        self.root.title(
            "Please select where you are running the script from"
        )  # Set the title of the window
        label = tk.Label(self.root, text="Location", font=helv)  # Create a label widget
        label.pack()  # Display the label in the window
        option1 = tk.Radiobutton(
            self.root, text="VM", variable=self.var, value="VM", font=helv
        )  # Create the first radio button
        option1.pack()  # Display the first radio button in the window
        option2 = tk.Radiobutton(
            self.root, text="S Drive", variable=self.var, value="S Drive", font=helv
        )  # Create the second radio button
        option2.pack()  # Display the second radio button in the window
        submit_button = tk.Button(
            self.root, text="Submit", command=self.submit, font=helv
        )  # Create the Submit button
        submit_button.pack()  # Display the Submit button in the window
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        self.root.mainloop()  # Start the tkinter event loop

    def on_closing(self) -> None:
        """
        Write to log and exit script upon closing the message box
        """
        if not self.choice:
            self.logger.info("No input location provided, exiting script")
            sys.exit(1)

    def submit(self) -> None:
        """
        Runs upon submission of message box submit button. Displays the selected choice
        in a message box, saves the choice, and closes the window
        """
        self.choice = self.var.get()
        # You can save the choice as a variable or perform any other action here
        messagebox.showinfo("Selection", f"You selected: {self.choice}")
        self.root.destroy()
        if self.choice:
            logger.info(f"Choice saved: {self.choice}")
        else:
            logger.info(f"No input was selected. Exiting")
            sys.exit(1)


if __name__ == "__main__":
    parsed_args = arg_parse()
    logfile_path = os.path.join(
        os.getcwd(),
        f"{datetime.datetime.now():%y%m%d_%H%M%S}_file_mover_postnatal.log",
    )
    logger = Logger(logfile_path).logger
    logger.info("Running file_mover_postnatal.py")
    message_box = MessageBox(logger)
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
