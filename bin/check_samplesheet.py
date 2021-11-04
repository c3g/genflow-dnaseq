#!/usr/bin/env python

# TODO nf-core: Update the script to check the samplesheet
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/dnaseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# TODO nf-core: Update the check_samplesheet function
def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,fastq_1,fastq_2
    SAMPLE_PE,SAMPLE_PE_RUN1_1.fastq.gz,SAMPLE_PE_RUN1_2.fastq.gz
    SAMPLE_PE,SAMPLE_PE_RUN2_1.fastq.gz,SAMPLE_PE_RUN2_2.fastq.gz
    SAMPLE_SE,SAMPLE_SE_RUN1_1.fastq.gz,

    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
    """

    id_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 2
        HEADER = ["Sample", "Readset", "Library", "RunType", "Run", "Lane", "Adapter1", "Adapter2", "QualityOffset", "BED", "FASTQ1", "FASTQ2", "BAM"]
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip('\n').split("\t")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {} but found {})!".format(len(HEADER), len(lspl)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, id, library, runtype, run, lane, adapter1, adapter2, qualityoffset, bed, fastq_1, fastq_2, bam = lspl[: len(HEADER)]
            id = id.replace(" ", "_")
            if not id:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, sample, library, runtype, run, lane, adapter1, adapter2, qualityoffset, bed, fastq_1, fastq_2]
            if id and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", sample, id, library, runtype, run, lane, adapter1, adapter2, qualityoffset, bed, fastq_1, fastq_2]
            elif id and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", sample, id, library, runtype, run, lane, adapter1, adapter2, qualityoffset, bed, fastq_1]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2 ] }
            if id not in id_mapping_dict:
                id_mapping_dict[id] = [sample_info]
            else:
                if sample_info in id_mapping_dict[id]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    id_mapping_dict[id].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(id_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["id", "single_end", "sample", "readset", "library", "runtype", "run", "lane", "adapter1", "adapter2", "qualityoffset", "bed", "fastq_1", "fastq_2"]) + "\n")
            for id in sorted(id_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == id_mapping_dict[id][0][0] for x in id_mapping_dict[id]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "ID: {}".format(id))

                for idx, val in enumerate(id_mapping_dict[id]):
                    fout.write(",".join(["{}_T{}".format(id, idx + 1)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
