import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


SCAFFOLD = "TGTAGCTTCTTTCGAGTACAAAAAC"
PROMOTOR = "TCCCAGATTATATCTATCACTGATAGGGATCGCAAATCCCCAGGTCAGAGGGCTATTTTCCCTGGTCAGA"


def get_rev_constructs_from_top_oligo_df(df):
    target = []
    rev_complement = []
    rev_complement_construct = []

    for i in df.index:
        top_oligo = df.loc[i]["anticipated sgRNA"]
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


df = pd.read_excel("~/seds/projects/crispri_data_process/rodrigo_sgrna_validation/2101UNHX-0238.xlsx",
                   sheet_name='Stat')

df = get_rev_constructs_from_top_oligo_df(df)

rev_lib_file_name = "rogdrigo_sgrna_02_2021.rev.fasta"
sequences = []
for i in df.index:
    guide = df.loc[i]
    record = SeqRecord(Seq(guide["rev_complement_construct"]), id=guide["Sample ID"], description="")
    sequences.append(record)

with open(rev_lib_file_name, "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")
