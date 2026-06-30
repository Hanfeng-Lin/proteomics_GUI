"""Loading, decontamination, group assignment and control clean-up.

Ports notebook cells 2, 3 and 5. The numerical operations (contaminant
filtering, regex group matching, control-group drop rule) are unchanged; the
only difference is that inputs/outputs are passed explicitly instead of living
in the global namespace.
"""

import re
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def load_dataset(config, contaminants):
    """Read the DIA-NN protein/peptide matrices and remove contaminants.

    Equivalent to notebook cell 2. Returns ``(df, df_peptide)``.
    """
    pg_matrix = getattr(config, "pg_path", None) or (config.file + ".pg_matrix.tsv")
    pr_matrix = getattr(config, "pr_path", None) or (config.file + ".pr_matrix.tsv")

    df = pd.read_csv(pg_matrix, sep="\t", index_col=0)
    df_peptide = pd.read_csv(pr_matrix, sep="\t", index_col=0)
    logger.info("protein before decontamination: %s", df.shape)

    df = df[~df.index.isin(contaminants)]
    logger.info("protein after decontamination: %s", df.shape)

    logger.info("peptide before decontamination: %s", df_peptide.shape)
    df_peptide = df_peptide[~df_peptide.index.isin(contaminants)]
    logger.info("peptide after decontamination: %s", df_peptide.shape)

    return df, df_peptide


def assign_groups(df, group_names, group_patterns=None):
    """Map each group name to the matching sample columns.

    By default each group *name* is itself the regex (notebook cell 3 behaviour).
    If ``group_patterns`` is given, a group listed there is matched with its
    (usually anchored, exact) regex while the clean name stays the dict key -- so
    e.g. label ``NR`` can match only ``..._NR_...`` and not ``..._NR-TDxdR_...``.
    Returns the ``group_columns`` dict.
    """
    group_patterns = group_patterns or {}
    group_columns = {}
    for group in group_names:
        pattern = group_patterns.get(group, group)
        group_columns[group] = [x for x in df.columns if re.search(pattern, x)]
    # Only strip regex punctuation from the keys when the names themselves are the
    # patterns; explicit group_patterns means the names are already clean labels.
    if not group_patterns:
        group_columns = {re.sub(r'\\.|\|\^|\$', '', k): v for k, v in group_columns.items()}

    for key in group_columns:
        logger.info("%s: %d", key, len(group_columns[key]))
    return group_columns


def control_group_cleanup(df, config, group_columns):
    """Drop unstably-detected proteins from the control group (cell 5).

    Only active when ``mode == 0`` and ``comparison_matrix`` is empty, exactly as
    in the notebook. Mutates and returns ``df``.
    """
    if config.mode == 0 and not config.comparison_matrix:
        control_group = group_columns[config.reference_group]

        before_clean_up_number = str(df.shape)
        for protein in df.index:
            control_values = df.loc[protein, control_group]
            if control_values.isna().sum() >= len(control_group) * (1 - config.control_group_detection_threshold):
                # Discard this protein
                df.drop(index=protein, inplace=True)

        logger.info("protein before control group clean-up: %s", before_clean_up_number)
        logger.info("protein after control group clean-up: %s", df.shape)

    if config.mode == 1 and not config.comparison_matrix:
        # All the treatment groups will be reference group.
        logger.info("Control group clean-up will not be performed in mode 1.")

    return df
