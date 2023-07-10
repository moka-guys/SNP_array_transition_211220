"""
This script is used to generate the BED file for the overlap map array project which
contains the positions of all probes to be masked within the CHAS array sofware.

The BED file is used with the 'overlap map' feature in CHAS which hides calls in regions
specified in a BED file. This is required because the array team needs a way to hide
'problematic' regions that have been called numerous times but contain no protein coding
genes - this reduces the number of calls needing interpretation so speeds up the
analysis and checking process.
"""
from __future__ import annotations
import os
import subprocess
from pathlib import Path
import pandas as pd
import duckdb


class GenerateBed():
    """
    This class contains the methods for creating the BED file containing coordinates of
    all probes to be masked

    Attributes
        num_probes (int):                   Number of probes away from the nearest
                                            coding region that shouldn't be masked
        probes_bed (str):                   Path to BED file containing all probe names
                                            and coordinates
        genes_aed (str):                    Path to Genes.aed file containing a curated
                                            list of genes and their coordinates which
                                            specifies whether the genes are coding or
                                            non-coding. This is used to identify which
                                            probes require masking
        outdir (str):                       Directory path to output files to
        grch38_fai_file (str):              Path to grch38 fasta index file
        grch38_genome_file (str):           Path to Genome file for bedtools. Generated
                                            from grch38_fai_file by
                                            self.get_genome_file()
        coding_regions_bed (str):           Path to BED file containing only coding
                                            regions
        coding_regions_bed_sorted (str):    Path to BED file containing single per-gene
                                            regions
        noncoding_regions_bed (str):        Path to BED file containing the non-coding
                                            regions outside protein-coding genes
        noncoding_probes (str):             Path to BED file containing only probes
                                            within the noncoding regions
        probes_to_mask_bed (str):           Path to BED file containing final probes
                                            that require masking in the CHAS software
        noncoding_regions_df (class):       Pandas dataframe containing the non-coding
                                            regions outside protein-coding genes
        noncoding_probes_df (class):        Pandas dataframe containing only probes
                                            within the noncoding regions
        probes_to_mask_df (class):          Pandas dataframe containing final probes
                                            that require masking in the CHAS software

    Methods
        get_genome_file()
            Run subprocess to generate .genome file from genome.fa.fai
        get_coding_regions_bed()
            Read the genes_aed file and manipulate to produce a BED file containing only
            coding regions
        filter_to_coding_regions(regions_df)
            Filter genes_aed dataframe so that it contains only coding regions, and drop
            category and strand columns
        condense_regions(coding_regions_df)
            Condense regions that overlap into single per-gene regions
        get_noncoding_regions_bed()
            Run bedtools complement to find the regions in the genome file that are not
            represented in the coding regions BED file (i.e. the non coding regions
            outside protein coding genes)
        get_noncoding_probes()
            Generate dataframe that contains all probes within the noncoding regions in
            the self.noncoding_regions_bed BED file
        remove_probes()
            Remove any probes that are within self.num_probes distance of a coding
            region
        match_probes_to_regions()
            Split all probes into 'per-region' groups for the noncoding regions in
            self.noncoding_regions_df
        cleanup_intermediate_files()
            Remove all intermediate BED files
        write_to_final_csv()
            Write to final csv, with header

    """
    def __init__(self, num_probes: int, probes_bed: str, genes_aed: str, outdir: str):
        """
        """
        self.num_probes = num_probes
        self.probes_bed = probes_bed
        self.genes_aed = genes_aed
        self.outdir = outdir
        self.grch38_fai_file = f"{Path(__file__).parent.resolve()}/data/genome.fa.fai"
        # Intermediate files
        self.grch38_genome_file = f"{outdir}/grch38.genome"
        self.coding_regions_bed = f"{outdir}/coding_regions_bed.bed"
        self.coding_regions_bed_sorted = f"{outdir}/coding_regions_bed_sorted.bed"
        self.noncoding_regions_bed = f"{outdir}/noncoding_regions_bed.bed"
        self.noncoding_probes = f"{outdir}/non_coding_probes.bed"
        self.probes_to_mask_bed = f"{outdir}/probes_to_mask.bed"
        # Call methods
        self.get_genome_file()
        self.get_coding_regions_bed()
        self.noncoding_regions_df = self.get_noncoding_regions_bed()
        self.noncoding_probes_df = self.get_noncoding_probes()
        self.probes_to_mask_df = self.remove_probes()
        self.cleanup_intermediate_files()
        self.write_to_final_csv()

    def get_genome_file(self) -> None:
        """
        Run subprocess to generate .genome file from genome.fa.fai. Extract column 1 and
        2 from the genome.fa.fai file and write to .genome file
            :return None:
        """
        genome_cmd = (
            "awk -v OFS='\t' {'print $1,$2'} " +
            f"{self.grch38_fai_file} > {self.grch38_genome_file}"
        )
        subprocess.run(genome_cmd, shell=True, check=True)

    def get_coding_regions_bed(self) -> None:
        """
        Read the genes_aed file and manipulate to produce a BED file containing only
        coding regions
            :return None:
        """
        headers = [
            "Chr",
            "Start",
            "Stop",
            "name",
            "value",
            "strand",
            "category",
            "spans",
            "cdsMin",
            "refseq:accessionNumber",
            "gnomad:percentHI",
            "refseq:geneName",
            "gnomad:pLI",
            "cdsMax",
        ]
        # Dataframe containing both coding and non-coding regions
        regions_df = pd.read_csv(
            self.genes_aed,
            sep="\t",
            skiprows=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],  # Don't read in lines 1-11
            names=headers,
            index_col="name",
        ).drop(  # Drop superfluous columns
            columns=[
                "spans",
                "cdsMin",
                "refseq:accessionNumber",
                "gnomad:percentHI",
                "refseq:geneName",
                "gnomad:pLI",
                "cdsMax",
                "value",
            ]
        )
        # If strand is negative reverse start and stop values
        regions_df.loc[regions_df["strand"] == "-", ["Start", "Stop"]] = regions_df.loc[
            regions_df["strand"] == "-", ["Stop", "Start"]
        ].values
        # Filter df to only coding regions
        coding_regions_df = self.filter_to_coding_regions(regions_df)
        # Condense coding regions to per-gene regions
        self.condense_regions(coding_regions_df)

    def filter_to_coding_regions(self, regions_df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter genes_aed dataframe so that it contains only coding regions, and drop
        category and strand columns
            :param regions_df (class):          Pandas dataframe containing all regions
                                                from the genes_aed file
            :return coding_regions_df (class):  Pandas dataframe containing only coding
                                                regions from the genes_aed file
        """
        coding_regions_df = regions_df[  # Keep only refseq/coding regions
            regions_df["category"].str.contains("refseq/coding")
        ].drop(
            columns=["category", "strand"]  # Drop unnecessary columns
            # Order and drop duplicate rows
            ).sort_values(by=['name']).drop_duplicates(inplace=False)
        return coding_regions_df

    def condense_regions(self, coding_regions_df: pd.DataFrame) -> None:
        """
        Condense regions that overlap into single per-gene regions and write to BED file
            :param coding_regions_df (class):   Pandas dataframe containing only coding
                                                regions
        """
        # Group by gene name and Chr, and aggregate into a single per-gene region, as we
        # want to find only regions at least x probes away from the nearest coding gene
        genes_start_stop = coding_regions_df.groupby(
            ["name", "Chr"], group_keys=False
            ).agg({"Start": min, "Stop": max})
        # Drop gene name column and re-index
        genes_start_stop.index = genes_start_stop.index.droplevel(0)

        # Write to BED file and sort
        genes_start_stop.to_csv(self.coding_regions_bed, sep="\t", header=False)
        subprocess.run(
            f"sort -k1,1V -k2,2n -k3,3n {self.coding_regions_bed} > "
            f"{self.coding_regions_bed_sorted}",
            shell=True,
            check=True,
        )

    def get_noncoding_regions_bed(self) -> None:
        """
        Run bedtools complement to find the regions in the genome file that are not
        represented in the coding regions BED file (i.e. the non coding regions) outside
        protein coding genes
            :return noncoding_regions_df (class):   Pandas dataframe containing the
                                                    non-coding regions outside protein-
                                                    coding genes
        """
        intersect_command = (
            f"bedtools complement -i {self.coding_regions_bed_sorted} -g "
            f"{self.grch38_genome_file} > {self.noncoding_regions_bed}"
        )
        subprocess.run(intersect_command, shell=True, check=True)
        # Read in noncoding regions dataframe
        noncoding_regions_df = pd.read_csv(
            self.noncoding_regions_bed, sep="\t", names=["Chr", "Start", "Stop"]
        )
        return noncoding_regions_df

    def get_noncoding_probes(self) -> pd.DataFrame:
        """
        Generate dataframe that contains all probes within the noncoding regions in the
        self.noncoding_regions_bed BED file
            :return noncoding_probes_df (class):    Pandas dataframe containing only
                                                    probes within the noncoding regions
        """
        intersect_cmd = (
            f"bedtools intersect -a {self.probes_bed} -b {self.noncoding_regions_bed} "
            f"> {self.noncoding_probes}"
        )
        subprocess.run(intersect_cmd, shell=True, check=True)

        # Read in noncoding probes dataframe generated by intersect command
        headers = ["Chr", "Pos", "Stop", "Probe"]
        noncoding_probes_df = pd.read_csv(
            self.noncoding_probes,
            sep="\t",
            names=headers,
        ).drop(columns=["Stop"])
        return noncoding_probes_df

    def remove_probes(self) -> pd.DataFrame:
        """
        Generate the final list of probes by removing any probes that are within
        self.num_probes_distance of a coding region. For each group, remove
        self.num_probes number of probes with the lowest and self.num_probes number of
        probes with the highest position, and write the remaining probes (probes to be
        masked) to BED file
            :return probes_to_mask_df (class):  Pandas dataframe containing final probes
                                                that require masking in the CHAS
                                                software
        """
        # Split probes into 'per-region' groups by noncoding region
        probes_within_regions = self.match_probes_to_regions()

        min_region_size = 2*self.num_probes
        # Remove any probes with size 2*self.num_probes (means that no probes are
        # included that are < self.num_probes away from the nearest coding region)
        probes_within_regions = probes_within_regions.groupby(["range", "Chr"]).filter(
            lambda x: len(x) >= min_region_size
            )
        # Use range (self.probes-1)-(self.probes-1) because index starts from 0. Meaning
        # only the slice from the self.probes probe away from the last coding region to
        # the self.probes probe away from the next coding region are retained
        probes_to_mask_df = probes_within_regions.groupby(
            ["range", "Chr"], group_keys=False).apply(
                lambda x: x.iloc[(self.num_probes-1):-(self.num_probes-1)]
        ).reset_index()

        probes_to_mask_df['Start'] = probes_to_mask_df['Pos']
        probes_to_mask_df['Stop'] = probes_to_mask_df['Pos']
        probes_to_mask_df.drop(["range", "Pos"], axis=1, inplace=True)

        return probes_to_mask_df

    def match_probes_to_regions(self) -> pd.DataFrame:
        """
        Split all probes into 'per-region' groups for the noncoding regions in
        self.noncoding_regions_df
            :return probes_within_regions (class):  Pandas dataframe containing each
                                                    probe mapped to the region within
                                                    which it falls
        """
        noncoding_probes_df = self.noncoding_probes_df
        noncoding_regions_df = self.noncoding_regions_df
        # Join the noncoding probes to the noncoding regions which the probe falls
        # within
        probes_within_regions = (
            duckdb.query(
                "SELECT * FROM noncoding_probes_df FULL OUTER JOIN noncoding_regions_df"
                " ON noncoding_probes_df.Chr = noncoding_regions_df.Chr WHERE (Pos >= "
                "Start AND Pos <= Stop)"
            )
            .df()
            .set_index("Probe")
        )  # Return a dataframe with Probe column as index
        # Add a range column
        probes_within_regions["range"] = (
            probes_within_regions["Start"].astype(str) + "," +
            probes_within_regions["Stop"].astype(str)
        )
        # Sort dataframe within each group by Pos column
        probes_within_regions = (
            probes_within_regions.groupby(["range", "Chr"])
            .apply(lambda x: x.sort_values(by=["range", "Chr", "Pos"], ascending=True))
            .drop(columns=["Start", "Stop", "Chr_2", "range", "Chr"])
        )
        return probes_within_regions

    def write_to_final_csv(self) -> None:
        """
        Write to final csv, with header
            :return None:
        """
        with open(self.probes_to_mask_bed, 'w', encoding="utf-8") as probes_bed:
            probes_bed.write('track db="hg38"\n')
        self.probes_to_mask_df.to_csv(
            self.probes_to_mask_bed, sep="\t", header=None, index=False, mode='a',
            columns=["Chr", "Start", "Stop", "Probe"],
            )

    def cleanup_intermediate_files(self) -> None:
        """
        Remove all intermediate BED files
            :return None:
        """
        for file in [
            self.coding_regions_bed, self.coding_regions_bed_sorted,
            self.noncoding_regions_bed, self.noncoding_probes, self.grch38_genome_file
        ]:
            os.remove(file)
