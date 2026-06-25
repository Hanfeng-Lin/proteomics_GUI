"""Differential-statistics: limma (cell 9) and the empirical-FDR curve used by
the volcano plot (cell 10).

limma fits a per-protein linear model (lmFit) and shrinks variances with
empirical Bayes (eBayes), so its p-value is a *moderated* t-test. ``group_columns``
and ``output_adjpval`` are explicit arguments; R / rpy2 is initialised lazily.
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


# --------------------------------------------------------------------------- #
# limma (notebook cell 9) -- R is loaded lazily on first use
# --------------------------------------------------------------------------- #
_r = None
_robjects = None
_pandas2ri = None


def _ensure_r():
    """Import rpy2, activate pandas<->R conversion and load limma exactly once."""
    global _r, _robjects, _pandas2ri
    if _r is None:
        from rpy2.robjects import pandas2ri, r
        import rpy2.robjects as robjects
        # Activate automatic DataFrame conversion between pandas and R
        pandas2ri.activate()
        r('library(limma)')
        _r, _robjects, _pandas2ri = r, robjects, pandas2ri
    return _r, _robjects, _pandas2ri


def limma_differential_analysis(df, treated_group_name, control_group_name, group_columns,
                                min_valid_value=3):
    r, robjects, pandas2ri = _ensure_r()

    # Extract relevant sample columns
    treated_samples = group_columns[treated_group_name]
    control_samples = group_columns[control_group_name]
    selected_columns = treated_samples + control_samples

    # Filter out rows with insufficient valid values
    valid_rows = df.dropna(thresh=min_valid_value, subset=selected_columns)

    # Convert group names to valid R variable names
    treated_group_r = r('make.names')(treated_group_name)[0]
    control_group_r = r('make.names')(control_group_name)[0]

    # Construct the design matrix in R
    sample_info = [treated_group_r] * len(treated_samples) + [control_group_r] * len(control_samples)
    sample_info_r = robjects.StrVector(sample_info)
    r.assign("sample_info", sample_info_r)
    r("design_matrix <- model.matrix(~ 0 + factor(sample_info))")
    r("colnames(design_matrix) <- levels(factor(sample_info))")

    # Convert DataFrame to R object
    df_r = pandas2ri.py2rpy(valid_rows[selected_columns].astype(float))
    r.assign("expression_data", df_r)

    # Fit limma model
    r("fit <- lmFit(expression_data, design_matrix)")
    r("fit <- eBayes(fit)")

    # Fix: Ensure contrast formula is a valid R expression
    contrast_str = f'"{treated_group_r} - {control_group_r}"'
    r(f"contrast_matrix <- makeContrasts({contrast_str}, levels=design_matrix)")
    r("fit2 <- contrasts.fit(fit, contrast_matrix)")
    r("fit2 <- eBayes(fit2)")
    r("limma_results <- topTable(fit2, adjust='BH', number=Inf, sort.by='none')")

    # Convert results back to pandas. topTable returns both P.Value (raw,
    # moderated-t) and adj.P.Val (BH/FDR); we always keep both, as their own
    # columns -- Pvalue_<comp> (raw) and bh_FDR_<comp> (adjusted).
    pvalues_df = pandas2ri.rpy2py(r('limma_results'))
    raw_p_column_name = f'Pvalue_{treated_group_name}_vs_{control_group_name}'
    pvalue_column_name = f'bh_FDR_{treated_group_name}_vs_{control_group_name}'
    log10pvalue_column_name = f'-log_P_adj_{treated_group_name}_vs_{control_group_name}'

    pvalues_df[raw_p_column_name] = pvalues_df['P.Value']
    pvalues_df = pvalues_df.rename(columns={'adj.P.Val': pvalue_column_name})
    pvalues_df[log10pvalue_column_name] = -np.log10(pvalues_df[pvalue_column_name])
    return pvalues_df[[raw_p_column_name, pvalue_column_name, log10pvalue_column_name]]


# --------------------------------------------------------------------------- #
# Empirical-FDR curve (defined in notebook cell 10, used by the volcano plot)
# --------------------------------------------------------------------------- #
def empirical_fdr_curve(df, alpha=0.05, kappa=0.01, p_value_cutoff=1, mode="neg",
                        x_col='log2FC', p_col='pval',
                        x0_grid=None, c_grid=None):
    x = df[x_col].values
    mlogp = -np.log10(df[p_col].values)

    if x0_grid is None:
        x0_grid = np.round(np.linspace(0.0, 1.0, 26), 3)    # 0 ~ 1, step 0.04
    if c_grid is None:
        c_grid = np.round(np.linspace(0.5, 8.0, 76), 3)     # 0.5 ~ 8

    records = []
    for x0 in x0_grid:
        pos_denom = np.maximum(x - x0, 1e-12)      # Used only when x>=x0
        neg_denom = np.maximum(-x - x0, 1e-12)     # Used only when x<=-x0

        for c in c_grid:
            pos_mask = (x >= x0) & (mlogp >= c/pos_denom + p_value_cutoff)
            neg_mask = (x <= -x0) & (mlogp >= c/neg_denom + p_value_cutoff)
            npos = int(pos_mask.sum())
            nneg = int(neg_mask.sum())
            if mode == "neg":
                fdr_hat = (npos + kappa) / (nneg + kappa)
                feasible = (nneg > 0) and (fdr_hat <= alpha)
                records.append((x0, c, npos, nneg, fdr_hat, feasible))
            if mode == "pos":
                fdr_hat = (nneg + kappa) / (npos + kappa)
                feasible = (npos > 0) and (fdr_hat <= alpha)
                records.append((x0, c, npos, nneg, fdr_hat, feasible))

    res = pd.DataFrame(records, columns=['x0', 'c', 'npos', 'nneg', 'FDR', 'feasible'])
    feas = res[res['feasible']]
    if len(feas) == 0:
        best = res.sort_values(['FDR'], ascending=[True]).iloc[0].to_dict()
    else:
        if mode == "neg":
            best = feas.sort_values(['nneg', 'x0', 'c'], ascending=[False, True, True]).iloc[0].to_dict()
        if mode == "pos":
            best = feas.sort_values(['npos', 'x0', 'c'], ascending=[False, True, True]).iloc[0].to_dict()
    return res, best
