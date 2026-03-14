import pandas as pd
import random
import os
import sys
import argparse
sys.path.insert(0, 'packages')
import order_pre_server
import order_pre_meta_server
from predictions_multiclass_kfold_server import all_func

## Variables
random_seed = 42
random.seed(random_seed)
top_path = "../../"
feat_table_path = top_path + "processing_sensitive_20072023/metaphlan/"
meta_path = top_path + "meta_after_r_266.csv"
top_save_path = top_path + "cst_predictions_output/"
prevalence_cutoff = 0.05
avg_abundance_cutoff = 0.05
unclassified_threshold = 80

## Parameters
parser = argparse.ArgumentParser(description="Run predictions using meta, metaphlan tables separated and merged for prediction of a pregnancy parameter")
parser.add_argument("remove_rare", help="rem if remove rare species function else not_rem (rem/not_rem)")
parser.add_argument("pred_col", help="column in metadata to be predicted, i.e., cst or lacto_bin")
args = parser.parse_args()

save_path = top_save_path + args.pred_col + "/rem_rare_" + str(args.remove_rare) + "_prevalence_cutoff_" + str(prevalence_cutoff) + "_avg_abundance_cutoff_" + str(avg_abundance_cutoff) + "_unclassified_thresh_" + str(unclassified_threshold)
if not os.path.exists(save_path):
    os.makedirs(save_path)

# Variables
selection = True
tuning = True
kmrmr = 20
k = 5
n_reps = 3
model_type = "RF"
run_clr = True
gini = False

# Save parameters
lines = ["Readme", "pred_col = " + args.pred_col, "selection = " + str(selection), "kmrmr = " + str(kmrmr), "k = " + str(k), "n_reps = " + str(n_reps),
        "model_type = " + str(model_type), "tuning = " + str(tuning), "clr = " + str(run_clr), "gini = " + str(gini),
        "remove_rare = " + str(args.remove_rare), "unclassified_threshold = " + str(unclassified_threshold), "prevalence_cutoff = " + str(prevalence_cutoff), "avg_abundance_cutoff = " + str(avg_abundance_cutoff),
	    "meta path = " + meta_path,
	    "feature table path = " + feat_table_path]
with open(save_path + '/readme.txt', 'w') as f:
    f.write('\n'.join(lines))

# Running options
META = True
DF_ALONE = True
DF_META = True

## Load data
meta = pd.read_csv(meta_path)

## Order data
pred_meta, pred_vag = order_pre_meta_server.order_meta(meta, args.pred_col, feat_table_path) # args.pred_col, feat_table_path)
print(f"Number of samples in meta: {len(pred_meta)}")

## Load metaphlan data
metaphlan_dict = {}
for level in ["sgb", "species", "genus", "phylum"]:
    organ_dict = {}
    for organ in ["rectal"]: # "vaginal", 
        path = feat_table_path + organ + "_metaphlan_" + level + ".tsv"
        df = pd.read_csv(path, sep="\t")
        ord_df = order_pre_server.order_metaphlan_df(df, feat_table_path, run_clr = run_clr, organ = organ, remove_rare=args.remove_rare, prevalence_cutoff=prevalence_cutoff, avg_abundance_cutoff=avg_abundance_cutoff, unclassified_threshold = unclassified_threshold)
        common = ord_df.index.intersection(pred_meta.index)
        ord_df  = ord_df.loc[common]
        print(f"In level {level} and organ {organ}, number of samples after matching with meta: {len(ord_df)} and number of features: {ord_df.shape[1]}")
        organ_dict[organ] = ord_df
    metaphlan_dict[level] = organ_dict

pred_meta = pred_meta.loc[common]
pred_vag = pred_vag.loc[common]
print(f"Number of samples in meta after matching with metaphlan: {len(pred_meta)}")

## Predictions
if META:
    res_key = "meta_no_selection"
    print(res_key)
    path = save_path + "/" + res_key + ".xlsx"
    meta_dict = all_func(pred_meta, pred_vag, path, kmrmr=kmrmr, k = k, n_reps = n_reps, select=False, tuning=tuning, gini = gini)
    print(res_key + ": " + str(meta_dict["all_df"].roc_auc.mean()))

if DF_ALONE:
    for level, lst in metaphlan_dict.items():
        df = lst["rectal"]
        print("Working on: " + level)
        res_key = "rectal_" + level
        path = save_path + "/" + res_key + ".xlsx"
        rectal_level_dict = all_func(df, pred_vag, path, kmrmr=kmrmr, k=k, n_reps=n_reps, select=selection, tuning=tuning, gini=gini)
        print(res_key + ": " + str(rectal_level_dict["all_df"].roc_auc.mean()))

if DF_META:
    for level, lst in metaphlan_dict.items():
        df = lst["rectal"]
        meta_df = pred_meta.merge(df, left_index=True, right_index=True)
        print("Number of features in meta + rectal " + level + ": " + str(meta_df.shape[1]))
        print("Working on: " + level)
        res_key = "meta_rectal_" + level
        path = save_path + "/" + res_key + ".xlsx"
        meta_rectal_level_dict = all_func(meta_df, pred_vag, path, kmrmr=kmrmr, k=k, n_reps=n_reps, select=selection, tuning=tuning, gini=gini)
        print(res_key + ": " + str(meta_rectal_level_dict["all_df"].roc_auc.mean()))