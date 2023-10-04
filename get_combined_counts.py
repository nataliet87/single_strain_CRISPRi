import pandas as pd
import os

def get_counts_combined_align(out_dir):
    counts_files = pd.Series(os.listdir(out_dir))[
        pd.Series(os.listdir(out_dir)).str.contains('counts')].tolist()

    df = pd.DataFrame(columns=['construct'])

    for file in counts_files:
        df_new = pd.read_csv(os.path.join(out_dir, file), delimiter='\t', 
                      names=['construct', file.rstrip('.counts')])

        df = df.merge(df_new, how='outer', on='construct')

    return df
