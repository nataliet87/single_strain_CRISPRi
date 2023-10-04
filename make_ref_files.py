import os
import subprocess
import argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
## RUN THIS IN PYTHON 3.7 WGS CONDA ENV

def parse_metadata_file(filepath):
    """

    :param filepath: (str) path to metadata file for sequencing run for analysis
    :return: df with sample name, target, rev_complement, and rev_complement_construct
    """
    scaffold = "TGTAGCTTCTTTCGAGTACAAAAAC"
    promotor = "TCCCAGATTATATCTATCACTGATAGGGATCGCAAATCCCCAGGTCAGAGGGCTATTTTCCCTGGTCAGA"

    df = pd.read_excel(filepath, sheet_name='Sheet1')
    df = df[['sample', 'sgRNA Sequence']]

    target = []
    rev_complement = []
    rev_complement_construct = []

    for i in df.index:
        top_oligo = df.loc[i]["sgRNA Sequence"]
        if top_oligo[:4] == 'GGGA':
            guide = top_oligo[4:]
            rev = str(Seq(guide).reverse_complement())

            target.append(guide)
            rev_complement.append(rev)
            rev_complement_construct.append(scaffold + rev + promotor)

        else:
            target.append(top_oligo)
            rev_complement.append(top_oligo)
            rev_complement_construct.append(top_oligo)

    df['target'] = target
    df['rev_complement'] = rev_complement
    df['rev_complement_construct'] = rev_complement_construct

    return df


def write_fasta_file(record, path):
    with open(path, "w") as output_handle:
        SeqIO.write(record, output_handle, "fasta-2line")


def build_subread_index(basename, reference, gappedIndex=False, indexSplit=False, memory=8000,
                TH_subread=100, colorspace=False, stdout=None, stderr=None, force=False):
    """

    :param basename:
    :param reference:
    :param gappedIndex:
    :param indexSplit:
    :param memory:
    :param TH_subread:
    :param colorspace:
    :param stdout:
    :param stderr:
    :param force:
    :return:
    """
    # system:
    # ./subread-buildindex [options] -o <basename> {FASTA file1} [FASTA file2]

    # Check if index already exists
    index_ext_list = ["log", "files", "00.b.tab", "00.b.array", "reads"]
    if sum([os.path.exists("%s.%s" % (basename, ext)) for ext in index_ext_list]) == len(index_ext_list):
        if not force:
            print('{} index files already exist and not forcing. Skipping.'.format(os.path.split(basename)[1]))
            return

    indexsplit_str = " -B" if not indexSplit else ""
    gappedindex_str = " -F" if not gappedIndex else ""
    colorspace_str = " -c" if colorspace else ""
    command_list = ["subread-buildindex", "-M", "{mem}{color}{gapped}{split}", "-f", "{th}", "-o", "{base}", "{ref}"]

    cmd_kwargs = dict(
        # Optional arguments
        mem=memory, th=TH_subread,
        # Optional Flags
        color=colorspace_str, gapped=gappedindex_str, split=indexsplit_str,
        # Required arguments
        base=basename, ref=reference)

    command = [c.format(**cmd_kwargs) for c in command_list]
    print(command)
    print("` ".join(command))
    subprocess.run(command) #, stdout=stdout, stderr=stderr)


def write_library_reference_files(df, output_dir):
    # iterate through CRISPRi samples, write a directory for each
    # and write a .fasta, .rev.fasta, and associated index files into that directory
    for sample in df.index:
        guide = df.loc[sample]

        # make directory to hold ref files for this sample
        os.system("mkdir {}".format(os.path.join(output_dir, guide['sample'])))

        # write target and rev construct fasta files
        record = SeqRecord(Seq(guide["target"]), id=guide["sample"], description="")
        file_name = os.path.join(output_dir, guide["sample"], guide["sample"] + ".fasta")

        rev_record = SeqRecord(Seq(guide["rev_complement_construct"]), id=guide["sample"], description="")
        rev_file_name = os.path.join(output_dir, guide["sample"], guide["sample"]+".rev.fasta")

        write_fasta_file(record, file_name)
        write_fasta_file(rev_record, rev_file_name)

        build_subread_index(rev_file_name, rev_file_name)
        # subread.build_index(rev_fasta_path, rev_fasta_path, gappedIndex=False, indexSplit=False, memory=8000)
        print("Reverse fasta, and index files written for {}".format(guide["sample"]))


def main():
    parser = argparse.ArgumentParser(description="""Script to generate a merged counts file after aligning multiple 
    samples against their specific, anticipated, sgRNA construct to generate CRISPRi counts data. """)

    parser.add_argument('--metadata', dest='metadata', type=str, required=True,
                        help='Filepath to metadata .csv file with sequencing sample names matched to sgRNA targets')

    parser.add_argument('--library_dir', dest='lib_out', type=str, required=True,
                        help='Filepath to directory to write target libraries to')

    args = parser.parse_args()

    # read in and parse metadata sheet into sample: CRISPRi_target mapping
    library_df = parse_metadata_file(args.metadata)

    # write CRISPRi_targets data to .fa, .rev.fa, and associated subread index files
    write_library_reference_files(library_df, args.lib_out)


if __name__ == "__main__":
    main()
