"""
Entrypoint for overlap_map_bed_generator
"""
import os
import argparse
import overlap_map_bed_generator.overlap_map_bed_generator as om


def msg(name=None):
    return '''python3 -m overlap_map_bed_generator [-h] -N NUM_PROBES -O OUTDIR -G GENES_AED -P PROBES_BED'''


def is_valid_dir(parser: argparse.ArgumentParser, arg: str) -> str:
    """
    Check directory path is valid
        :param parser (argparse.ArgumentParser):    Holds necessary info to parse cmd
                                                    line into Python data types
        :param arg (str):                           Input argument
    """
    if not os.path.isdir(arg):
        parser.error(f"The directory {arg} does not exist!")
    else:
        return arg  # Return argument


def is_valid_file(parser: argparse.ArgumentParser, arg: str) -> str:
    """
    Check file path is valid
        :param parser (argparse.ArgumentParser):    Holds necessary info to parse cmd
                                                    line into Python data types
        :param arg (str):                           Input argument
    """
    if not os.path.exists(arg):
        parser.error(f"The file {arg} does not exist!")
    else:
        return arg  # Return argument


def arg_parse() -> dict:
    """
    Parse arguments supplied by the command line. Create argument parser, define command
    line arguments, then parse supplied command line arguments using the created
    argument parser
        :return (dict): Parsed command line attributes
    """
    parser = argparse.ArgumentParser(
        description=(
            "Generate the BED file for the overlap map array project which contains "
            "the positions of all probes to be masked within the CHAS array sofware"
            ),
        usage=msg()
    )
    parser.add_argument(
        "-N",
        "--num_probes",
        type=int,
        help=(
            "Number of probes away from the nearest coding region that we do not want "
            "to be masked. E.g. if a value of 50 is input, the BED output will contain "
            "all probes fom the 51st probe and higher away from the nearest coding "
            "region"
            ),
        required=True,
    )
    parser.add_argument(
        "-O",
        "--outdir",
        type=lambda x: is_valid_dir(parser, x),
        help="Directory to output created files to",
        required=True,
    )
    parser.add_argument(
        "-G",
        "--genes_aed",
        type=lambda x: is_valid_file(parser, x),
        help=(
            "Genes.aed file containing a curated list of genes and their coordinates "
            "which specifies whether the genes are coding or non-coding. This is used "
            "to identify which probes require masking"
        ),
        required=True,
    )
    parser.add_argument(
        "-P",
        "--probes_bed",
        type=lambda x: is_valid_file(parser, x),
        help="BED file containing all probe names and coordinates",
        required=True,
    )
    return vars(parser.parse_args())


args = arg_parse()

om.GenerateBed(
    args["num_probes"], args["probes_bed"], args["genes_aed"], args["outdir"]
)
