#!/usr/bin/env bash

python -u main_testing.py \
--restore_model_path='./saved_models/saved_pretrained/snapshot_epoch_2000.pt' \
--rank_model_path='./saved_models/rank_model_perBin_1/best_epoch_18.pt' \
--restore_para_file='./saved_models/rank_model_perBin_1/Bin_1_pretrain_model_dict.pkl'\
--num_results=10 \
--cand_terms_path=''