import pandas as pd
import numpy as np
from random import sample
import pickle


def sample_terms(singlet_list, rate=0.01):

    num = round(len(singlet_list) * rate)

    return sample(singlet_list, num)


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


def create_pickled_data():

    singlet = pd.read_csv("singlets_terms_perBin_1d.txt",
                          delimiter="\t", header=None)

    singlet = singlet.rename(columns={0: "term", 1: "singleton_freq"})

    singlet_list = singlet['term'].to_list()

    sampled_singlet = sample_terms(singlet_list)
    print(len(sampled_singlet))

    out_df = calculate_PPMI()

    dataset = {}
    i = 0
    for singleton in sampled_singlet:
        if i % 50 == 0:
            print(i)

        out_list = out_df[out_df['term1'] == singleton][[
            'term2', 'PPMI', 'co_occur_freq']].values.tolist()
        out_tuples = [tuple(x) for x in out_list]

        dataset[singleton] = out_tuples
        i += 1

    with open('sub_neighbors_dict_ppmi_perBin_1.pkl', 'wb') as handle:
        pickle.dump(dataset, handle,
                    protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    create_pickled_data()
