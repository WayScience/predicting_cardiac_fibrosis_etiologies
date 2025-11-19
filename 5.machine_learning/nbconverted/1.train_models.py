#!/usr/bin/env python
# coding: utf-8

# # Train three models with subsequent shuffled baselines

# In[1]:


import pathlib
import sys
import warnings

import numpy as np
import pandas as pd
from joblib import dump
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.utils import parallel_backend

sys.path.append("../utils")
from training_utils import downsample_data, get_X_y_data


# In[2]:


# Set numpy seed for reproducibility
np.random.seed(0)

# Path to training/testing datasets
training_data_path = pathlib.Path("./data_splits")

# Find all training datasets
training_files = list(training_data_path.rglob("training_split.parquet"))

# Metadata column used for prediction class
label = "Metadata_cell_type"

# Directories for outputs
model_dir = pathlib.Path("./models")
model_dir.mkdir(exist_ok=True, parents=True)

encoder_dir = pathlib.Path("./encoder_results")
encoder_dir.mkdir(exist_ok=True, parents=True)

training_indices_dir = pathlib.Path("./training_indices")
training_indices_dir.mkdir(exist_ok=True, parents=True)

print(f"Found {len(training_files)} training datasets.")

# Dictionary to store loaded training datasets
training_data = {}

# Loop through and load each training dataset
for training_file in training_files:
    dataset_name = training_file.parent.name  # Use parent folder name as key
    print(f"Loading dataset: {dataset_name}")  # only print the model/folder name

    train_df = pd.read_parquet(training_file)
    training_data[dataset_name] = train_df


# In[3]:


# Loop through and downsample each loaded training dataset
for dataset_name, train_df in training_data.items():
    # Downsample to the smallest class
    downsample_df = downsample_data(data=train_df, label=label)

    # Replace or store the downsampled dataframe
    training_data[dataset_name] = downsample_df

    # Export sample indices used in training to CSV
    output_file = training_indices_dir / f"{dataset_name}_training_data_indices.csv"
    pd.DataFrame(downsample_df.index, columns=["Index"]).to_csv(
        output_file, index=False
    )

    print(f"CSV file created at {output_file} with {len(downsample_df.index)} entries.")
    print(downsample_df.shape)
    print(downsample_df[label].value_counts())


# In[4]:


# Collect all unique labels across all datasets
all_labels = set()
for dataset_name, train_df in training_data.items():
    all_labels.update(train_df[label].unique())

# Fit the LabelEncoder on the combined set of all labels
le = LabelEncoder()
le.fit(list(all_labels))

# Save the global label encoder for consistency
dump(le, encoder_dir / "label_encoder_global.joblib")

# Print the global class mapping
class_mapping = dict(zip(le.classes_, le.transform(le.classes_)))
print("Global Class Mapping:")
print(class_mapping)

# Process each dataset to get X and y
for dataset_name, train_df in training_data.items():
    # Non-shuffled data
    X_train, y_train = get_X_y_data(df=train_df, label=label, shuffle=False)
    y_train_encoded = le.transform(y_train)

    # Shuffled data
    X_shuffled_train, y_shuffled_train = get_X_y_data(
        df=train_df, label=label, shuffle=True
    )
    y_shuffled_train_encoded = le.transform(y_shuffled_train)

    # Store X and y in the dictionary
    training_data[dataset_name] = {
        "X_train": X_train,
        "y_train": y_train_encoded,
        "X_shuffled_train": X_shuffled_train,
        "y_shuffled_train": y_shuffled_train_encoded,
    }


# In[5]:


# Set folds for k-fold cross validation (default is 5, best for low sample size)
straified_k_folds = StratifiedKFold(n_splits=5, shuffle=False)

# Set Logistic Regression model parameters (use default for max_iter)
logreg_params = {
    "penalty": "elasticnet",
    "solver": "saga",
    "max_iter": 1000,
    "n_jobs": -1,
    "random_state": 0,
    "class_weight": "balanced",
}

# Define the hyperparameter search space for RandomizedSearchCV
param_dist = {
    "C": np.logspace(-2, 1, 7),  # values from 0.01 to 10, narrow for low sample size
    "l1_ratio": np.linspace(0, 1, 11),
}

# Set the random search hyperparameterization method parameters (used default for "cv" and "n_iter" parameter)
random_search_params = {
    "param_distributions": param_dist,
    "scoring": "f1_weighted",
    "random_state": 0,
    "n_jobs": -1,
    "cv": straified_k_folds,
}


# In[6]:


# Initialize Logistic Regression and RandomizedSearchCV
logreg = LogisticRegression(**logreg_params)
random_search = RandomizedSearchCV(logreg, **random_search_params)

# Loop through the training data dictionary
for dataset_name, data_dict in training_data.items():
    X_train = data_dict["X_train"]
    y_train = data_dict["y_train"]
    X_shuffled_train = data_dict["X_shuffled_train"]
    y_shuffled_train = data_dict["y_shuffled_train"]

    with parallel_backend("multiprocessing"):
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore", category=ConvergenceWarning, module="sklearn"
            )

            ########################################################
            # Train model on non-shuffled (final) training data
            ########################################################
            print(f"Training model for {dataset_name} (final)...")
            final_random_search = clone(random_search)
            final_random_search.fit(X_train, y_train)
            print(
                f"Optimal parameters for {dataset_name} (final):",
                final_random_search.best_params_,
            )

            # Save model
            final_model_filename = model_dir / f"{dataset_name}_final_downsample.joblib"
            dump(final_random_search.best_estimator_, final_model_filename)
            print(f"Model saved as: {final_model_filename}")

            ########################################################
            # Train model on shuffled training data
            ########################################################
            print(f"Training model for {dataset_name} (shuffled)...")
            shuffled_random_search = clone(random_search)
            shuffled_random_search.fit(X_shuffled_train, y_shuffled_train)
            print(
                f"Optimal parameters for {dataset_name} (shuffled):",
                shuffled_random_search.best_params_,
            )

            # Save model
            shuffled_final_model_filename = (
                model_dir / f"{dataset_name}_shuffled_downsample.joblib"
            )
            dump(shuffled_random_search.best_estimator_, shuffled_final_model_filename)
            print(f"Model saved as: {shuffled_final_model_filename}")

