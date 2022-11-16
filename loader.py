import os
import re
import sys
import pickle
import distance
import jellyfish
import numpy as np
import networkx as nx
from collections import Counter

import torch

import utils


def load_pretrain_vec(embed_filename, dim=None):
    word_dict = {}
    with open(embed_filename) as f:
        for idx, line in enumerate(f):
            L = line.strip().split()
            word, vec = L[0], L[1::]
            # word = L[0].lower()
            if dim is None and len(vec) > 1:
                dim = len(vec)
            elif len(vec) == 1:
                print('header? ', L)
                continue
            elif dim != len(vec):
                raise RuntimeError('Wrong dimension!')

            word_dict[word] = np.array(vec, dtype=np.float32)
            # assert(len(word_dict[word]) == input_dim)
    return word_dict


def load_pretrain_graph_embed(file_name):
    f = open(file_name).readlines()
    # print('Node embeddings: ', f[0].strip())
    node_num = int(f[0].strip().split()[0])
    node_embed = int(f[0].strip().split()[1])
    node_dict = {int(x.strip().split()[0]): np.array(x.strip().split()[1::], dtype=np.float32)
                 for x in f[1::]}

    return node_dict, node_embed, node_num


def w2v_mapping(list_phrases):
    phrases_list = [x.split() for x in list_phrases]
    words_list = [x.lower() for y in phrases_list for x in y]
    word_vocab = Counter(words_list)
    word_vocab = [x[0] for x in word_vocab.items() if x[1] > 10]
    word_to_id = {x: idx+1 for idx, x in enumerate(word_vocab)}
    word_to_id['<UNK>'] = len(word_to_id) + 1
    id_to_word = {v: k for k, v in word_to_id.items()}
    print('Find {0} words.'.format(len(word_to_id)))
    return word_vocab, word_to_id, id_to_word


def pre_train_word_mapping(word_to_id, args):

    def load_glove_vec(embed_filename, vocab, input_dim):
        word_dict = {}
        with open(embed_filename) as f:
            for idx, line in enumerate(f):
                L = line.split()
                word = L[0].lower()
                if word in vocab:
                    word_dict[word] = np.array(L[-input_dim::], dtype=np.float32)
                    # assert(len(word_dict[word]) == input_dim)
        return word_dict

    word_dict = load_glove_vec(args.embed_filename, word_to_id.keys(), args.word_embed_dim)

    initrange = 0.5 / args.word_embed_dim
    W = np.random.uniform(-initrange, initrange, (len(word_to_id) + 1, args.word_embed_dim))
    i = 0
    for cur_word in word_to_id.keys():
        if cur_word in word_dict:
            i += 1
            W[word_to_id[cur_word]] = word_dict[cur_word]

    print(i / float(len(word_to_id)))
    return W


def w2v_mapping_pretrain(list_phrases, args):
    pretrain_dict = load_pretrain_vec(args.embed_filename, args.word_embed_dim)
    pretrain_word_vocab = pretrain_dict.keys()

    phrases_list = [x.split() for x in list_phrases]
    words_list = [x.lower() for y in phrases_list for x in y]
    full_word_vocab = Counter(words_list)

    word_vocab = [x[0] for x in full_word_vocab.items() if x[0] in pretrain_word_vocab or x[1] > 20]
    # word_vocab = [x[0] for x in full_word_vocab.items() if x[1] > 10]

    word_to_id = {x: idx + 1 for idx, x in enumerate(word_vocab)}  # skip 0 id
    word_to_id['<UNK>'] = len(word_to_id) + 1
    id_to_word = {v: k for k, v in word_to_id.items()}

    initrange = 0.5 / args.word_embed_dim
    W = np.random.uniform(-initrange, initrange, (len(word_to_id) + 1, args.word_embed_dim))
    W[0] = np.zeros(args.word_embed_dim)
    i = 0
    for cur_word in word_to_id.keys():
        if cur_word in pretrain_dict:
            i += 1
            W[word_to_id[cur_word]] = pretrain_dict[cur_word]

    print('Find {0} words with pretrain ratio: {1}'.format(len(word_to_id), i / float(len(word_to_id))))

    return word_vocab, word_to_id, id_to_word, W


def c2v_mapping(list_words):
    # [str1, str2, str3, ...]
    chars_list = list(' '.join(list_words))
    char_vocab = Counter(chars_list)
    char_to_id = {x[0]: idx+1 for idx, x in enumerate(char_vocab.items())}
    char_to_id['<UNK>'] = len(char_to_id) + 1
    char_to_id['<s>'] = len(char_to_id) + 1
    char_to_id['</s>'] = len(char_to_id) + 1
    # counting frequency of characters
    id_to_char = {v: k for k, v in char_to_id.items()}
    assert((' ' in char_to_id) and (0 not in id_to_char))
    print('Find %i character.' % (len(char_to_id)))
    return char_vocab, char_to_id, id_to_char


def g2v_mapping_pretrain(list_words, args):
    gram_dict = load_pretrain_vec(args.ngram_embed_path, args.ngram_embed_dim)
    pretrain_gram_vocab = gram_dict.keys()  # 874474

    ngram_list = []
    list_words = list(set(list_words))
    for w in list_words:
        ngram_list += utils.get_single_ngrams(w, args.n_grams)

    # for n in args.n_grams:
    #     cur_list = [ngrams_pretrain(w, n) for w in list_words]
    #     ngram_list += [g for x in cur_list for g in x]

    full_ngrams_vocab = Counter(ngram_list)
    ngrams_vocab = [x[0] for x in full_ngrams_vocab.items() if x[0] in pretrain_gram_vocab or x[1] > 100]
    # ngrams_vocab = [x[0] for x in full_ngrams_vocab.items() if x[0] in pretrain_gram_vocab]

    ngrams_to_id = {x: idx + 1 for idx, x in enumerate(ngrams_vocab)}
    ngrams_to_id['<UNK>'] = len(ngrams_to_id) + 1
    id_to_ngrams = {v: k for k, v in ngrams_to_id.items()}

    initrange = 0.5 / args.ngram_embed_dim
    W = np.random.uniform(-initrange, initrange, (len(ngrams_to_id) + 1, args.ngram_embed_dim))
    W[0] = np.zeros(args.ngram_embed_dim)
    W[-1] = np.zeros(args.ngram_embed_dim)

    i = 0
    for cur_gram in ngrams_to_id.keys():
        if cur_gram in pretrain_gram_vocab:
            i += 1
            W[ngrams_to_id[cur_gram], :] = gram_dict[cur_gram]

    print('Find {0} grams with pretrain ratio: {1}'.format(len(ngrams_to_id), i / float(len(ngrams_to_id))))

    return ngrams_vocab, ngrams_to_id, id_to_ngrams, W


def node_mapping(args, node_to_id=None):
    f = open(args.node_embed_path).readlines()
    # print('Node embeddings: ', f[0].strip())
    node_dict = {int(x.strip().split()[0]): np.array(x.strip().split()[1::], dtype=np.float32)
                 for x in f[1::]}
    # print(len(node_dict))
    if not node_to_id:
        terms_mat = np.zeros((len(node_dict) + 1, int(f[0].strip().split()[1])))
        terms_to_idx = {-1: 0}
        idx_to_terms = {0: -1}
        for idx, term_id in enumerate(node_dict.keys()):
            terms_mat[idx + 1, :] = node_dict[term_id].reshape(1, -1)
            terms_to_idx[term_id] = idx + 1
            idx_to_terms[idx + 1] = term_id

        return terms_to_idx, idx_to_terms, terms_mat
    else:
        terms_mat = np.zeros((len(node_to_id), int(f[0].strip().split()[1])))
        for term_id, idx in node_to_id.items():
            terms_mat[idx, :] = node_dict[term_id].reshape(-1)

        return terms_mat
