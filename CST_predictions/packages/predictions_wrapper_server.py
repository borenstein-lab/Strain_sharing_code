####################################################
# The following script contains various functions for running predictions of pregnancy parameters  
# Last updated: 13/11/2024 
# Implemented functions:  
# - create_specific_run_save_path
# - get_common_sample_id
# - cols_kmrmr_wrapper_func 
# - add_frame 
# - save_excel
# - save_pred_cols
####################################################

import os
import sys
import pandas as pd
sys.path.insert(0, 'my_packages')
import predictions_kfold_server


## Save path
def create_specific_run_save_path(out_dir, pred_col, kmrmr, k, n_reps, model_type, selection, tuning, run_clr, gini, feat_table_path, meta_path):
    save_path = out_dir + "/" + pred_col + "/kmrmr_" + str(kmrmr) + "_k" + str(k) + "_n_reps" + str(n_reps) + "_" + model_type
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    lines = ["Readme", "pred_col = " + pred_col, "selection = " + str(selection), "kmrmr = " + str(kmrmr), "k = " + str(k), "n_reps = " + str(n_reps), 
            "model_type = " + str(model_type), "tuning = " + str(tuning), "clr = " + str(run_clr), "gini = " + str(gini),
	    "meta path = " + meta_path, 
	    "feature table path = " + feat_table_path]
    with open(save_path + '/readme.txt', 'w') as f:
        f.write('\n'.join(lines))
    return save_path

## Common sample_id
def get_common_sample_id(features_df, target_df, pred_col):
    target_col = target_df.dropna()
    common = features_df.index.intersection(target_col.index)
    return features_df.loc[common], target_col.loc[common]

## All functions wrapper
def cols_kmrmr_wrapper_func(features_df, target_df, out_dir, res_key, feat_table_path, meta_path, kmrmr_lst, k_lst = [3, 5, 10], n_reps = 3, model_type = "RF", selection = True, tuning = False, run_clr = False, gini = True):
    pred_cols_res_dict = {}        
    pred_col = target_df.name # columns[0]
    all_pred_col_res_df = pd.DataFrame(columns = ["pred_col", "res_key", "kmrmr", "k", "n_reps", "mean_aucroc", "mean_auspr"])
    curr_features_df, curr_target_col = get_common_sample_id(features_df, target_df, pred_col)
    for kmrmr in kmrmr_lst:
        for k in k_lst:
            save_path = create_specific_run_save_path(out_dir, pred_col, kmrmr, k, n_reps, model_type, selection, tuning, run_clr, gini, feat_table_path, meta_path)
            file_save_path = save_path + "/" + res_key + ".xlsx"
            curr_dict = predictions_kfold_server.all_func(curr_features_df, curr_target_col, file_save_path, kmrmr, k, n_reps, model_type, selection, tuning, gini)
            curr_results_lst = [pred_col, res_key, kmrmr, k, n_reps, curr_dict["all_df"].auc.mean(), curr_dict["pr_df"].auc.mean()]
            all_pred_col_res_df.loc[len(all_pred_col_res_df.index)] = curr_results_lst
        pred_cols_res_dict[pred_col] = all_pred_col_res_df
    return pred_cols_res_dict

## Add dataframe values to each pred_col dataframe
def add_frame(add_dict, final_dict):
    for key, df in add_dict.items():
        final_df = pd.concat([final_dict[key], df])
        final_dict[key] = final_df
    return final_dict

## Save dictionary as excel
def save_excel(path, dict_res):
    writer = pd.ExcelWriter(path, engine='xlsxwriter')
    for key, value in dict_res.items():
        value.to_excel(writer, sheet_name=str(key))
    writer.save()
    writer.close()
    return

## Save pred_cols dictionary in each directory
def save_pred_cols_dict(final_dict, out_dir):
    for pred_col, df in final_dict.items():
        save_path = out_dir + "/" + pred_col 
        df.sort_values('mean_aucroc').to_csv(save_path + "/all_res_keys_" + pred_col + ".csv")
        df[["kmrmr", "mean_aucroc", "mean_auspr"]] = df[["kmrmr", "mean_aucroc", "mean_auspr"]].astype(float)
        print(df[["res_key", "kmrmr", "mean_aucroc", "mean_auspr"]].nlargest(3, 'mean_aucroc'))
    return