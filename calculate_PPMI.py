import pandas as pd
import numpy as np


def calculate_PPMI():

    co_occ = pd.read_csv("cofreqs_terms_perBin_1d.txt",
                         delimiter="\t", header=None)
    co_occ = co_occ.rename(
        columns={0: "term1", 1: "term2", 2: "co_occur_freq"})

    singlet = pd.read_csv("singlets_terms_perBin_1d.txt",
                          delimiter="\t", header=None)
    singlet = singlet.rename(columns={0: "term", 1: "singleton_freq"})

    intermed = co_occ.merge(singlet, left_on="term1", right_on="term")
    intermed = intermed.rename(
        columns={"singleton_freq": "term1_freq"}).drop(columns='term')

    out_df = intermed.merge(singlet, left_on="term2", right_on="term").rename(
        columns={"singleton_freq": "term2_freq"}).drop(columns='term')
    out_df['PMI'] = np.log((out_df['co_occur_freq'] * len(singlet)) /
                           (out_df['term1_freq'] * out_df['term2_freq']))
    out_df['PPMI'] = np.where(out_df['PMI'] < 0, 0, out_df['PMI'])

    return out_df
