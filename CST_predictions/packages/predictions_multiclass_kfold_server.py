####################################################
# The following script contains various functions for running predictions of pregnancy parameters  
# Last updated: 30/12/2025
# Implemented functions:  
# - get_importance_array
# - multiclass_cros_val
# - add_mean_std_metrics
# - all_func
####################################################

import pandas as pd
import os
import numpy as np
from sklearn.metrics import roc_curve, auc, precision_recall_curve, average_precision_score
from sklearn.preprocessing import label_binarize
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedKFold, RepeatedStratifiedKFold, KFold, GridSearchCV
from sklearn.inspection import permutation_importance
from xgboost import XGBClassifier
from sklearn import metrics
from mrmr import mrmr_classif, mrmr_regression
from sklearn.multiclass import OneVsRestClassifier
import copy
from predictions_kfold_server import get_common_sample_id, feature_select, hyper_tuning, save_excel

# Create importance array
def get_importance_array(importance_df, X, X_train, X_test, y_test, classifier, label, fold, gini = True):
    importance_array = np.full(X.shape[1], np.nan)
    indices = [X.columns.get_loc(col) for col in X_train.columns]
    if gini:
        importance_values = list(classifier.feature_importances_)
    else:
        importance_bunch = permutation_importance(classifier, X_test, y_test, n_repeats=10, random_state=42, n_jobs=32)
        importance_values = importance_bunch["importances_mean"]
    importance_array[indices] = importance_values
    importance_dict = dict(zip(X.columns, importance_array))
    importance_row = pd.DataFrame([importance_dict], columns=importance_dict.keys())
    importance_row = pd.DataFrame({'label': [label], 'fold': [fold]}).join(importance_row)
    return pd.concat([importance_df, importance_row])

# Multiclass repeated Kfold cross validation
def multiclass_cros_val(X, y, k=5, n_reps=3, select=True, kmrmr=12, tuning = True, gini = False):
    cv = RepeatedStratifiedKFold(n_splits=k, n_repeats=n_reps, random_state=42)
    unique_classes = np.unique(y)
    y_bin = label_binarize(y, classes=unique_classes)
    label_to_bin_index = {label: idx for idx, label in enumerate(unique_classes)}
    model = RandomForestClassifier(random_state=42)
    y_reals, y_probas = [], []
    roc_data, pr_data = [], []
    importance_df = pd.DataFrame(columns = ["label", "fold"] + X.columns.to_list())
    for fold, (train, test) in enumerate(cv.split(X, y)):
        X_train = X.iloc[train]
        y_train = y.iloc[train]
        if select and X.shape[1] > kmrmr:
            X_train, y_train = feature_select(X_train, y_train, kmrmr)
        if tuning:
            model = hyper_tuning(model, "RF", X_train, y_train)
        model.fit(X_train, y_train)
        X_test = X.iloc[test].loc[:, X_train.columns]
        y_test = y.iloc[test]
        y_proba = model.predict_proba(X_test)
        for label, i in label_to_bin_index.items():
            fpr, tpr, _ = roc_curve(y_bin[test, i], y_proba[:, i])
            roc_auc = auc(fpr, tpr)
            precision, recall, _ = precision_recall_curve(y_bin[test, i], y_proba[:, i])
            pr_auc = average_precision_score(y_bin[test, i], y_proba[:, i])
            roc_data.extend(zip([label] * len(fpr), [fold] * len(fpr), fpr, tpr, [roc_auc] * len(fpr)))
            pr_data.extend(zip([label] * len(precision), [fold] * len(precision), precision, recall, [pr_auc] * len(precision)))
            importance_df = get_importance_array(importance_df, X, X_train, X_test, y_test, model, label, fold, gini=gini)
    # Append to all res lists
    roc_df = pd.DataFrame.from_records(roc_data, columns=["label", "fold", "fpr", "tpr", "roc_auc"])
    pr_df = pd.DataFrame.from_records(pr_data, columns=["label", "fold", "precision", "recall", "pr_auc"])
    y_reals.append(y_test)
    y_probas.append(y_proba)
    return roc_df, pr_df, importance_df, y_reals, y_probas

# Mean of fpr, tpr, and auc (or precision, recall, and pr auc)
def add_mean_std_metrics(df, metric_x, metric_y, auc_col):
    """Compute mean and standard deviation for metric_y over interpolated metric_x values.
    Works for both ROC (FPR/TPR) and Precision-Recall (Recall/Precision) curves."""
    mean_fold_name = df['fold'].max() + 1  # Assign a new fold index for mean values
    mean_std_data = []
    for label in df["label"].unique():
        mean_x = np.linspace(0, 1, 100)
        interp_y, aucs = [], []
        for fold in df["fold"].unique():
            subset = df[(df["label"] == label) & (df["fold"] == fold)]
            if subset.empty:
                continue
            subset = subset.sort_values(by=metric_x)
            fold_interp_y = np.interp(mean_x, subset[metric_x], subset[metric_y])
            interp_y.append(fold_interp_y)
            aucs.append(subset[auc_col].iloc[0])
        mean_y = np.mean(interp_y, axis=0)
        std_y = np.std(interp_y, axis=0)
        mean_auc = np.mean(aucs)
        if metric_x == "fpr" and metric_y == "tpr":
            mean_y[-1] = 1.0
        mean_std_data.extend(zip([label] * len(mean_x), [mean_fold_name] * len(mean_x), mean_x, mean_y, [mean_auc] * len(mean_x), std_y))
    mean_std_df = pd.DataFrame.from_records(mean_std_data, columns=["label", "fold", metric_x, metric_y, auc_col, "std"])
    return pd.concat([df, mean_std_df], ignore_index=True)

# All functions together
def all_func(features_df, target_col, path, kmrmr=12, k = 5, n_reps = 3, select=True, tuning = True, gini = True):
    y = target_col
    X = features_df
    X, y = get_common_sample_id(X, y)
    # K-fold cross validation
    roc_df, pr_df, importance_df, y_reals, y_probas = multiclass_cros_val(X, y, k, n_reps, select, kmrmr, tuning = tuning, gini = gini)
    # Add mean values
    roc_df = add_mean_std_metrics(roc_df, "fpr", "tpr", "roc_auc")
    pr_df = add_mean_std_metrics(pr_df, "recall", "precision", "pr_auc")
    # Create final dictionary
    save_dict = {'all_df': roc_df, 'pr_df': pr_df, "importance_df": importance_df}
    save_excel(path, save_dict)
    return save_dict

