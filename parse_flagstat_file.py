import argparse
import pandas as pd
import subprocess
import os


def parse_flagstat_file(file):
    file = open(file, mode='r', encoding='utf-8-sig')
    lines = file.readlines()
    file.close()

    data = []

    for i in range(0, len(lines), 14):
        data.append([lines[i], lines[i + 1], lines[i + 5]])

    df = pd.DataFrame(data, columns=['sample', 'total_reads', 'mapped_reads'])

    df['sample'] = df['sample'].apply(lambda x: os.path.split(x.rstrip())[1].split('.')[0])
    df['total_reads'] = df['total_reads'].apply(lambda x: int(x.split('+')[0].strip()))
    df['mapped_reads'] = df['mapped_reads'].apply(lambda x: int(x.split('+')[0].strip()))

    return df


def get_counts_combined_align(out_dir):
    """

    :param out_dir:
    :return:
    """
    counts_files = pd.Series(os.listdir(out_dir))[pd.Series(os.listdir(out_dir)).str.contains('counts')].tolist()

    df = pd.DataFrame(columns=['construct'])

    for file in counts_files:

        df_new = pd.read_csv(os.path.join(out_dir, file), delimiter='\t', names=['construct', file.rstrip('.counts')])
        df = df.merge(df_new, how='outer', on='construct')

    return df


def main():
    parser = argparse.ArgumentParser(description="""Script to write a raw flagstat report file (holding a list of bam
    filenames and the samtools flagstat output generated from those files), and then parse that flagstat report into 
    a readable .csv file.""")

    parser.add_argument('--bam_list', dest='bam_list', type=str, required=True,
                        help="File holding sorted bam filepaths")

    parser.add_argument('--flagstat_report_path', dest='flagstat_raw_path', default='./', type=str, required=False,
                        help="Directory to write raw flagstat_report.txt into if desired. (Writes flagstat_report.txt "
                             "to working directory by default)")

    parser.add_argument('--out', dest='output_name', type=str, required=True,
                        help="Filename to write reads csv file to")

    parser.add_argument('--counts', dest='counts', type=str, required=False,
                        help="Pass this flag with path to directory holding counts file, if want to generate a "
                             "combined counts.csv file merged with aligned read counts data; REVIEW THIS -- NOT SURE"
                             "IF IT'S WORKING")

    args = parser.parse_args()

    flagstat_file = os.path.join(args.flagstat_raw_path, "flagstat_report.txt")

    subprocess.run(["/Users/nataliethornton/PycharmProjects/single_strain_CRISPRi/generate_raw_flagstat_file.sh",
                    "-I", args.bam_list, "-O", flagstat_file])

    df = parse_flagstat_file(flagstat_file)

    df.to_csv(args.output_name, index=False)

    if args.counts:
        counts = get_counts_combined_align(args.counts)
        # counts.to_csv("{}_counts.csv".format(args.output_name.strip('csv').strip(".")), index=False)


if __name__ == "__main__":
    main()
