import os
import pysam
import argparse
from crispr_tools.counting import subread


def align_fastq(fq, ref, out_dir):

    out_name = os.path.split(fq)[-1].split('_R1_00')[0]
    output_bam = os.path.join(out_dir, out_name+".bam")

    subread.align(ref, fq, output_bam, type="dna", output_format="bam", maxMismatches=1)
    pysam.sort("-o", output_bam+".sorted.bam", "-m", "2G", output_bam+"bam")
    pysam.index(output_bam+".sorted.bam")


def main():
    parser = argparse.ArgumentParser(description="""Script to align, sort and index a fastq against a reference using Subread. """)

    parser.add_argument('--in', dest='fq', type=str, required=True, help="Filepath to fastq")
    parser.add_argument('--ref', dest='ref', type=str, required=True, help="Filepath reference index base")
    parser.add_argument('--out_dir', dest='out_dir', type=str, required=True, help="Filepath to output directory")

    args = parser.parse_args()
    align_fastq(args.fq, args.ref, args.out_dir)


if __name__ == "__main__":
    main()
