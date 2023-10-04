from crispr_tools.counting import subread
from os import listdir
from os.path import isfile, join


def _get_trimmed_output_path(fastq_path):
    if fastq_path.endswith(".fastq"):
        return fastq_path.replace(".fastq", ".trimmed.fastq")
    elif fastq_path.endswith(".fastq.gz"):
        return fastq_path.replace(".fastq.gz", ".trimmed.fastq.gz")

def trim_and_align():
    fastq_dir = '../unmerged_fastq'
    output_dir = 'output_bams'

    fastqs = [f for f in listdir(fastq_dir) if isfile(join(fastq_dir, f))]
    
    for fq in fastqs:
        basename = fq.split('.')[0].rstrip("R1").rstrip("_")
        lib_path = join('references', basename+'.rev.fasta')
        output_bam = join(output_dir, basename+'.bam')


        trimmed_output = _get_trimmed_output_path(join(fastq_dir ,fq))


        upstream_adaptor = CRISPRI_UPSTREAM[:len(CRISPRI_UPSTREAM) - adaptor_size_left] ## i5 adaptor
        dowstream_adaptor = CRISPRI_DOWNSTREAM[adaptor_size_left:(adaptor_size_left + len(upstream_adaptor))] ## i7 adaptor
        _trim_sequences_from_fastq(fastq, trimmed_output,
                                   upstream_seq=upstream_adaptor,
                                   downstream_seq=dowstream_adaptor)


        subread.align(lib_path, join(fastq_dir, fq), output_bam, type="dna",
                       output_format="bam", maxMismatches=1)


def main():

    
    
if __name__ == "__main__":
    main()
