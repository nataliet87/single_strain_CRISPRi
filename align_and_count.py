import os
from os.path import isfile, join

import pandas as pd
import pysam
import numpy
import argparse

from crispr_tools import read_fasta
from crispr_tools.counting import is_good_counts, count_bam_file, \
    good_enough_barcode_qual, good_enough_cigar, good_enough_qual, \
    good_enough_sgrna, write_counts, write_diagnostics, subread

# Running on python2.7 in pebble-dev env


def get_output_paths(sample):
    output_bam = os.path.join(sample["output_path"], "%s.bam" % (sample["Sample"]))
    output_sorted_bam = os.path.join(sample["output_path"], "%s.sorted.bam" % (sample["Sample"]))
    output_counts_file = os.path.join(sample["output_path"], "%s.counts" % (sample["Sample"]))

    output_synth_file = output_counts_file.replace(".counts", ".synth")
    output_diag_file = output_counts_file.replace(".counts", ".diagnostics")

    return output_bam, output_sorted_bam, output_counts_file, output_synth_file, output_diag_file


# # set up a csv here with a sample name that matches up with fastq file names. this will then
# # be used to generate output name paths etc
# def parse_metadata_file(metadata_path, fastq_path, reference_path, output_path):
#     df = pd.read_csv(metadata_path)
#
#     df['fastq_path'] = df['file'].apply(lambda x: os.path.join(fastq_path, x))
#     df["ref_path"] = reference_path
#     # df["ref_path"] = df["P7"].apply(lambda x: os.path.join(reference_path, x, x+"_barcodes.fa"))
#     df["output_path"] = output_path
#
#     if len(df[df['fastq_path'] == '']) > 0:
#         print("!!! The following expected fastqs are not found:")
#         print(df[df['fastq_path'] == '']['Sample'])
#
#     df = df[df['fastq_path'] != '']
#
#     return df


def align_fastqs(sample_df):
    """

    :param sample_df:
    :return:
    """
    for i in range(len(sample_df)):
        sample = sample_df.iloc[i]
        print(sample)

        fq = sample['fastq_path']
        ref_path = sample['ref_path']
        output_bam = os.path.join(sample['output_path'], str(sample['Sample'])+'.bam')

        # basename = os.path.split(sample['fastq_path'])[1].split('_R1_')[0]
        # output_bam = os.path.join(sample['output_path'], basename+'.bam')
        # fq = join(sample['fastq_path'], sample['fastq_name'])

        print("Aligning {} against {}".format(sample['Sample'], ref_path))
        subread.align(ref_path, fq, output_bam, type="dna", output_format="bam", maxMismatches=1)

        print("---> Fastq files aligned")

        return


def sort_and_index_bams(sample_df):
    """

    :param sample_df:
    :return:
    """

    for i in range(len(sample_df)):
        sample = sample_df.iloc[i]

        bam = join(sample["output_path"], str(sample["Sample"])+".bam")
        sorted_bam = join(sample["output_path"], str(sample["Sample"])+".sorted.bam")

        print(bam)
        pysam.sort("-o", sorted_bam, "-m", "2G", bam)
        pysam.index(sorted_bam)

    print("---> Bam files sorted and indexed")

    return


def _get_run_type_qc_functions(sample_df):
    """

    :param sample_df:
    :return:
    """

    if sample_df.run_type in ["barcode", "bc"]:
        matchQCFunc = None
        qualQCFunc = good_enough_barcode_qual
        sgrnaQCFunc = None

    else:
        matchQCFunc = good_enough_cigar
        qualQCFunc = good_enough_qual
        sgrnaQCFunc = good_enough_sgrna

    return matchQCFunc, qualQCFunc, sgrnaQCFunc


def count_and_save(sample, force=False, alt_names=False):
    """

    :param sample:
    :param force:
    :param alt_names:
    :return:
    """

    output_bam, output_sorted_bam, output_counts_file, output_synth_file, output_diag_file = get_output_paths(sample)

    try:

        total_reads = numpy.sum([eval('+'.join(l.rsplit('\t')[2:])) for l in pysam.idxstats(output_sorted_bam)])
        # total_reads = numpy.sum([eval('+'.join(l.rsplit('\t')[2:])) for l in
        #                         pysam.idxstats(output_sorted_bam).strip().split("\n")])

        # bam_idx = pysam.idxstats(output_sorted_bam)
        # total_reads = pd.DataFrame([i.split("\t") for i in bam_idx.strip("\n").split("\n")])[[2, 3]].astype(int).sum().sum()

        # Adjust total after trimming if not trimming
        discarded_trim = 0

        # Get the counts from the sorted output bam file
        if force or not is_good_counts([output_counts_file, output_synth_file, output_diag_file]):

            try:
                sgRNA_counts, synth_error_counts, diagnostics = count_bam_file(
                    output_sorted_bam, matchQCFunc=None, qualQCFunc=None, sgrnaQCFunc=None,
                    DEBUG=None, total_reads=total_reads, count_synth=False)

                # Save counts, synth files and diagnostics file
                if alt_names:
                    names = [rawname for rawname, seq in read_fasta(sample.ref_path)]
                    alt = []
                    for i in names:
                        alt.append(i.split(" ")[0])

                    write_counts(sgRNA_counts, output_counts_file, name_list=alt)

                else:
                    write_counts(sgRNA_counts, output_counts_file, SGRNA_LIBRARY_REV=sample.ref_path)

                del sgRNA_counts
                del synth_error_counts

                write_diagnostics(diagnostics, output_diag_file)
                del diagnostics

            except Exception as e:
                # raise Exception('Error counting file (%s):\n%s' % (output_sorted_bam, e))
                print('Error counting file (%s):\n%s' % (output_sorted_bam, e))
                pass

        else:
            print("Error or counts files already exist")

    except pysam.SamtoolsError:
        print("!!! WARNING: BAM file could not be found. SKIPPING {} !!!".format(output_sorted_bam))

    return output_counts_file


def main():

    parser = argparse.ArgumentParser(description="""Script to generate a merged counts file after aligning multiple 
    samples against their specific, anticipated, sgRNA construct to generate CRISPRi counts data. """)

    parser.add_argument('--metadata', dest='metadata', type=str, required=True,
                        help='Filepath to metadata .csv file with sequencing sample names matched to sgRNA targets')

    parser.add_argument("--sort", dest='sort', action='store_true',
                        help='Pass this flag if bam files still need to be sorted and indexed')

    parser.add_argument("--align", dest='align', action='store_true', help='Pass this flag if alignment should be run')

    args = parser.parse_args()
    sample_df = pd.read_csv(args.metadata)
    # sample_df = parse_metadata_file(args.metadata, args.raw_fastq, args.lib_dir, args.out_dir) #args.trim_dir,

    # process trimmed fastqs
    if args.align:
        align_fastqs(sample_df)

    if args.sort:
        sort_and_index_bams(sample_df)

    for i in range(len(sample_df)):
        sample = sample_df.iloc[i]
        # count_and_save(sample, force=True, alt_names=True)
        count_and_save(sample)


if __name__ == "__main__":
    main()
