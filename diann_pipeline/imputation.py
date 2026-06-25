"""Imputation, normalization and fold-change helpers (notebook cell 6).

Function bodies are copied verbatim from the notebook. The only changes are to
the signatures: ``group_columns`` and ``df_peptide`` are now explicit arguments
instead of globals, and the ``mode``/``df_peptide`` default values that the
notebook bound from globals are passed in by the caller. No numerical operation
is altered. Imputation remains stochastic (unseeded ``truncnorm.rvs`` /
``np.random.uniform``), exactly as before.
"""

import logging

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from scipy.stats import truncnorm


def _check_peptide_significance(protein, df_peptide, treated_group, control_group):
    """
    Performs a peptide-level t-test for a single protein to decide if imputation is warranted.
    This helper function contains the logic that is difficult to vectorize fully.

    Args:
        protein (str): The name of the protein to check.
        df_peptide (pd.DataFrame): DataFrame containing peptide-level data.
        treated_group (list): List of column names for the treated group.
        control_group (list): List of column names for the control group.

    Returns:
        bool: True if the protein's missing values should be imputed based on peptide data.
    """
    try:
        protein_peptides = df_peptide.loc[protein]
        # Ensure protein_peptides is a DataFrame
        if isinstance(protein_peptides, pd.Series):
            protein_peptides = protein_peptides.to_frame().T

        # Create a copy to avoid SettingWithCopyWarning and ensure data is numeric.
        protein_peptides = protein_peptides.copy()
        # Coerce data to numeric type, turning any non-numeric values into NaN.
        for col in treated_group + control_group:
            if col in protein_peptides.columns:
                protein_peptides[col] = pd.to_numeric(protein_peptides[col], errors='coerce')

        # Check if there is any data at all in the treated group for this protein's peptides
        if protein_peptides[treated_group].isna().all().all():
            logging.info(f"Protein {protein}: No peptide data in treated group. Imputing.")
            return True

        # Calculate p-values and fold changes for each peptide
        p_values = []
        fold_changes = []
        for i in range(len(protein_peptides)):
            peptide_row = protein_peptides.iloc[i]

            # Explicitly cast to float right before the test to prevent dtype errors.
            treated_vals = peptide_row[treated_group].dropna().astype(float)
            control_vals = peptide_row[control_group].dropna().astype(float)

            # Calculate fold change if possible
            if not treated_vals.empty and not control_vals.empty and control_vals.mean() != 0:
                fold_changes.append(treated_vals.mean() / control_vals.mean())
            else:
                fold_changes.append(np.nan)

            # Perform t-test if possible
            if len(treated_vals) > 1 and len(control_vals) > 1:
                p_values.append(ttest_ind(treated_vals, control_vals, equal_var=True, nan_policy='omit').pvalue)
            else:
                p_values.append(np.nan)

        p_values = np.array(p_values)
        fold_changes = np.array(fold_changes)

        # Condition 1: Median p-value is significant
        if not np.isnan(p_values).all():
            significant_p = np.nanmedian(p_values) < 0.05
        else:
            significant_p = False
        # Condition 2: All p-values are NaN, but at least one fold change could be calculated
        # This handles cases where t-tests fail (e.g., n=1) but there's still evidence of presence.
        failed_ttest_with_fc = np.isnan(p_values).all() and not np.isnan(fold_changes).all()

        if significant_p or failed_ttest_with_fc:
            if np.isnan(p_values).all():
                logging.info(f"Protein {protein}: Imputing based on peptide data with failed t-tests but fold changes present.")
            else:
                logging.info(f"Protein {protein}: Imputing based on significant peptide data (p-median: {np.nanmedian(p_values):.4f}).")

            return True
        else:
            #logging.info(f"Protein {protein}: Skipping imputation, peptide data not significant.")
            return False

    except KeyError:
        # This protein is not in the peptide dataframe
        logging.info(f"Protein {protein} not found in peptide data. Cannot perform peptide-level check.")
        return False


def imputation(df, treated_group_name, control_group_name, group_columns, df_peptide,
               peptide_count_cutoff=3, mode=0):
    """
    Performs data imputation on a proteomics dataset using vectorized operations.

    Args:
        df (pd.DataFrame): The main protein abundance DataFrame. Index should be protein names.
        treated_group_name (str): The key for the treated group in group_columns.
        control_group_name (str): The key for the control group in group_columns.
        group_columns (dict): Dictionary mapping group names to lists of column names.
        df_peptide (pd.DataFrame): DataFrame with peptide-level data. Index should be protein names.
        peptide_count_cutoff (int): The minimum number of peptides required for certain imputations.
        mode (int): The imputation mode. 0 for degradation (impute treated), 1 for enrichment (impute control).

    Returns:
        tuple: A tuple containing:
            - pd.DataFrame: The imputed DataFrame.
            - list: A list of proteins that were imputed under the high-missing-value condition.
    """
    # --- 1. SETUP AND INITIAL CALCULATIONS ---
    control_group = group_columns[control_group_name]
    treated_group = group_columns[treated_group_name]

    df_imputed = df.copy()
    imputation_list = []

    all_cols = control_group + treated_group
    for col in all_cols:
        if col in df_imputed.columns:
            df_imputed[col] = pd.to_numeric(df_imputed[col], errors='coerce')

    # Use the 1st percentile of the entire dataset for low-level imputation, as it's a more stable baseline
    low_1_percentile = np.nanpercentile(df_imputed[all_cols].values.flatten(), 1)
    peptide_counts = df_peptide.index.value_counts()

    logging.info(f"Global low 1st percentile: {low_1_percentile:.4f}")
    if mode == 0:
        logging.info("Running in Mode 0: Degradation Protocol")
    else:
        logging.info("Running in Mode 1: Enrichment Protocol")

    # --- 2. PRE-FILTERING ---
    filter_group = control_group if mode == 0 else treated_group
    filter_nan_frac = df_imputed[filter_group].isna().sum(axis=1) / len(filter_group)
    proteins_to_drop_mask = filter_nan_frac >= 0.5
    if proteins_to_drop_mask.any():
        logging.info(f"Dropping {proteins_to_drop_mask.sum()} proteins with >=50% missing values in the {filter_group} group.")
        df_imputed = df_imputed[~proteins_to_drop_mask]

    # --- 3. LOGICAL SEPARATION OF PROTEINS FOR IMPUTATION ---
    control_nan_frac = df_imputed[control_group].isna().sum(axis=1) / len(control_group)
    treated_nan_frac = df_imputed[treated_group].isna().sum(axis=1) / len(treated_group)

    target_nan_frac = treated_nan_frac if mode == 0 else control_nan_frac

    # Identify proteins for the special high-missing value protocol
    special_imputation_mask = target_nan_frac >= 0.5
    special_imputation_idx = df_imputed.index[special_imputation_mask]

    # Identify proteins for the standard protocol (all others)
    standard_imputation_mask = ~special_imputation_mask
    standard_imputation_idx = df_imputed.index[standard_imputation_mask]

    # --- 4. STANDARD IMPUTATION on proteins with <50% missing values ---
    if not standard_imputation_idx.empty:
        logging.info(f"Performing standard imputation on {len(standard_imputation_idx)} proteins.")
        # Impute control group NaNs
        means = df_imputed.loc[standard_imputation_idx, control_group].mean(axis=1).fillna(0)
        stds = df_imputed.loc[standard_imputation_idx, control_group].std(axis=1).clip(lower=1e-6).fillna(1e-6)
        for col in control_group:
            nan_mask = df_imputed.loc[standard_imputation_idx, col].isna()
            if nan_mask.any():
                df_imputed.loc[nan_mask.index[nan_mask], col] = truncnorm.rvs(
                    a=(0 - means[nan_mask]) / stds[nan_mask], b=np.inf,
                    loc=means[nan_mask], scale=stds[nan_mask], size=nan_mask.sum())

        # Impute treated group NaNs
        means = df_imputed.loc[standard_imputation_idx, treated_group].mean(axis=1).fillna(0)
        stds = df_imputed.loc[standard_imputation_idx, treated_group].std(axis=1).clip(lower=1e-6).fillna(1e-6)
        for col in treated_group:
            nan_mask = df_imputed.loc[standard_imputation_idx, col].isna()
            if nan_mask.any():
                df_imputed.loc[nan_mask.index[nan_mask], col] = truncnorm.rvs(
                    a=(0 - means[nan_mask]) / stds[nan_mask], b=np.inf,
                    loc=means[nan_mask], scale=stds[nan_mask], size=nan_mask.sum())

    # --- 5. SPECIAL IMPUTATION on proteins with >=50% missing values ---
    if not special_imputation_idx.empty:
        #logging.info(f"Performing special imputation on {len(special_imputation_idx)} proteins.")
        target_group = treated_group if mode == 0 else control_group

        protein_peptide_counts = peptide_counts.reindex(special_imputation_idx).fillna(0)
        has_enough_peptides = protein_peptide_counts > peptide_count_cutoff

        proteins_to_check_idx = special_imputation_idx[has_enough_peptides]
        if not proteins_to_check_idx.empty:
            impute_decision = proteins_to_check_idx.to_series().apply(
                _check_peptide_significance, args=(df_peptide, treated_group, control_group))
            proteins_to_impute_idx = impute_decision.index[impute_decision]

            if not proteins_to_impute_idx.empty:
                imputation_list.extend(proteins_to_impute_idx.tolist())
                nan_frac_subset = target_nan_frac.loc[proteins_to_impute_idx]

                # Case 1: All values are missing. Overwrite with low-level uniform data.
                all_missing_mask = nan_frac_subset == 1.0
                proteins_all_missing_idx = nan_frac_subset[all_missing_mask].index
                if not proteins_all_missing_idx.empty:
                    imputed_data = np.random.uniform(low=low_1_percentile * 0.5, high=low_1_percentile * 1.5,
                                                     size=(len(proteins_all_missing_idx), len(target_group)))
                    df_imputed.loc[proteins_all_missing_idx, target_group] = imputed_data

                # Case 2: Partially missing. Impute ONLY the NaNs based on existing data for that protein.
                partial_missing_mask = nan_frac_subset < 1.0
                proteins_partial_missing_idx = nan_frac_subset[partial_missing_mask].index
                if not proteins_partial_missing_idx.empty:
                    means = df_imputed.loc[proteins_partial_missing_idx, target_group].mean(axis=1).fillna(0)
                    stds = df_imputed.loc[proteins_partial_missing_idx, target_group].std(axis=1).clip(lower=1e-6).fillna(1e-6)

                    for col in target_group:
                        nan_mask = df_imputed.loc[proteins_partial_missing_idx, col].isna()
                        if nan_mask.any():
                            df_imputed.loc[nan_mask.index[nan_mask], col] = truncnorm.rvs(
                                a=(0 - means[nan_mask]) / stds[nan_mask], b=np.inf,
                                loc=means[nan_mask], scale=stds[nan_mask], size=nan_mask.sum())

    logging.info(f"Imputation completed. {len(imputation_list)} proteins imputed via special high-missing-value protocol.")
    return df_imputed, imputation_list


def normalize_by_specific_protein(df, uniprot_id, group_columns, scaling_factor=100000000):
    """
    Normalizes the intensities of all proteins in the DataFrame based on the
    intensity of a specific reference protein (e.g., a housekeeping protein).

    Formula: Normalized_Intensity = (Target_Protein_Intensity / Reference_Protein_Intensity) * scaling_factor

    Args:
        df (pd.DataFrame): Imputed DataFrame (rows=proteins, cols=samples).
        uniprot_id (str): The UniProt ID of the reference protein.
        group_columns (dict): Dictionary mapping group names to lists of column names.
        scaling_factor (float): Optional multiplier to keep values in a readable range
                                (e.g., 1e6 or the mean of the reference). Default is 1.0.

    Returns:
        pd.DataFrame: A new DataFrame with normalized intensities.
    """
    df_normalized = df.copy()

    # 1. Check if the reference protein exists
    if uniprot_id not in df.index:
        logging.warning(f"Normalization skipped: Reference protein '{uniprot_id}' not found in dataset index.")
        return df_normalized

    # 2. Extract reference protein intensity (only for experimental columns)
    all_experimental_cols = [col for group in group_columns.values() for col in group]
    ref_intensities_exp = df.loc[uniprot_id].reindex(all_experimental_cols)

    # 3. Identify valid and invalid samples for normalization
    valid_cols = ref_intensities_exp.dropna().index.tolist()
    invalid_cols = ref_intensities_exp[ref_intensities_exp.isna()].index.tolist()

    # 4. Perform Normalization on Valid Columns
    if valid_cols:
        # Extract the subset of the dataframe and the reference intensities for valid columns
        df_valid = df[valid_cols].copy()
        ref_valid = ref_intensities_exp.loc[valid_cols]

        # Debugging: Check for zero values
        zero_count = (ref_valid == 0).sum()
        if zero_count > 0:
             logging.warning(f"Normalization issue: Reference protein '{uniprot_id}' has {zero_count} zero values in valid samples (will result in Inf/NaN).")

        # Perform normalization: division by reference, multiplied by scaling factor
        df_valid_normalized = df_valid.div(ref_valid, axis=1) * scaling_factor

        # Overwrite the original columns in the copy with the normalized values
        df_normalized[valid_cols] = df_valid_normalized[valid_cols]

        logging.info(f"Normalization performed on {len(valid_cols)} samples using reference protein: {uniprot_id}.")

    # 5. Copy Original Values for Invalid Columns
    if invalid_cols:
        # Values in invalid_cols are copied as they were (unnormalized, potentially imputed/raw)
        df_normalized[invalid_cols] = df[invalid_cols]
        logging.info(f"Skipped normalization for {len(invalid_cols)} samples due to NaN in reference protein.")

    # 6. Final Debugging
    if uniprot_id in df_normalized.index:
        mean_normalized_intensity = float(df_normalized[valid_cols].mean().mean())
        logging.info(f"Post-normalization mean intensity on valid samples: {mean_normalized_intensity:.4f}")
    return df_normalized


# Calculate FC value and append to the right side
def calculate_average_FC_value(df, treated_group_name, control_group_name, group_columns):
    control_group = group_columns[control_group_name]
    treated_group = group_columns[treated_group_name]

    control_avg = df[control_group].mean(axis=1)
    treated_avg = df[treated_group].mean(axis=1)

    FC_values = treated_avg / control_avg
    FC_column_name = f'FC_{treated_group_name}_vs_{control_group_name}'
    df[FC_column_name] = FC_values

    log2FC_values = np.log2(FC_values.astype('float64'))

    log2FC_column_name = f'log2FC_{treated_group_name}_vs_{control_group_name}'
    df[log2FC_column_name] = log2FC_values

    return df
