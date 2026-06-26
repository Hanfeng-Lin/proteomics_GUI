"""Pipeline orchestration (notebook cells 7 and 9).

These were the bare procedural cells that mutated globals (``imputed_dataframes``,
``imputation_dict``, ``df_final_results``). Here each becomes a function that
takes its inputs explicitly and returns its outputs. ``AnalysisResult`` bundles
the state the plotting functions need, replacing the global namespace.

``run_workflow`` (impute + normalize + fold change + Excel) then ``run_limma``
(differential statistics via R/limma -- a moderated t-test). The Student's
t-test fallback from the original notebook has been removed: statistics are
limma-only.
"""

import logging
from dataclasses import dataclass
from typing import Optional

import pandas as pd

from .config import AnalysisConfig
from .io import load_dataset, assign_groups, control_group_cleanup
from .imputation import (
    imputation,
    normalize_by_specific_protein,
    calculate_average_FC_value,
)
from .stats import limma_differential_analysis

logger = logging.getLogger(__name__)


@dataclass
class AnalysisResult:
    """All state produced by the core pipeline; consumed by the plotting funcs."""
    config: AnalysisConfig
    df_original: pd.DataFrame
    df_peptide: pd.DataFrame
    group_columns: dict
    imputed_dataframes: dict      # comparison_name -> imputed DataFrame
    imputation_dict: dict         # comparison_name -> list of imputed proteins
    summary: pd.DataFrame         # the notebook's df_final_results


# --------------------------------------------------------------------------- #
# cell 7: imputation + normalization + fold change + Excel export
# --------------------------------------------------------------------------- #
def _save_workflow_excel(config, df_original, group_columns, imputed_dataframes, df_final_results, imputation_dict):
    if config.imputation_option:
        if config.normalization_protein_id:
            output_excel_path = "final_analysis_results_imputed_normalized.xlsx"
        else:
            output_excel_path = "final_analysis_results_imputed.xlsx"
    else:
        if config.normalization_protein_id:
            output_excel_path = "final_analysis_results_normalized.xlsx"
        else:
            output_excel_path = "final_analysis_results.xlsx"

    # Identify metadata columns (columns that are not part of any experimental group).
    all_experimental_cols = [col for group in group_columns.values() for col in group]
    metadata_cols = [col for col in df_original.columns if col not in all_experimental_cols]

    logger.info(f"Saving results to {output_excel_path}...")
    with pd.ExcelWriter(output_excel_path, engine='openpyxl') as writer:
        # Sheet 1: the fold-change/stats summary. Surface the UniProt accession
        # (the index) as the first column.
        summary_out = df_final_results.copy()
        summary_out.index.name = "Protein.Group"
        summary_out.reset_index().to_excel(writer, sheet_name='Fold_Change_Summary', index=False)

        # Subsequent sheets: The raw imputed data for each comparison.
        for comparison_name, imputed_df in imputed_dataframes.items():
            treated_name, control_name = comparison_name.split('_vs_')
            fc_col = f'FC_{comparison_name}'
            log2fc_col = f'log2FC_{comparison_name}'
            control_cols = sorted(group_columns[control_name])
            treated_cols = sorted(group_columns[treated_name])

            # Build the subset from the imputed dataframe so imputed values are written.
            intensity_cols = metadata_cols + control_cols + treated_cols
            df_subset = imputed_df[intensity_cols].copy()
            df_subset[fc_col] = df_final_results[fc_col]
            df_subset[log2fc_col] = df_final_results[log2fc_col]

            # Include the limma statistics if present (this runs after run_limma).
            p_col = f'Pvalue_{comparison_name}'        # raw moderated-t p
            fdr_col = f'bh_FDR_{comparison_name}'      # adjusted P (BH/FDR)
            if p_col in df_final_results.columns:
                df_subset[p_col] = df_final_results[p_col]
            if fdr_col in df_final_results.columns:
                df_subset[fdr_col] = df_final_results[fdr_col]

            # Flag proteins imputed via the special high-missing-value protocol
            # (the same set marked orange on the volcano plot) as TRUE / FALSE.
            imputed_set = set(imputation_dict.get(comparison_name, []))
            df_subset["Imputed"] = [idx in imputed_set for idx in df_subset.index]

            # Truncate sheet name to Excel's 31-character limit.
            sheet_name = comparison_name[:31]

            # Save the subset, with the UniProt accession (index) as the first column.
            df_subset.index.name = "Protein.Group"
            df_subset.reset_index().to_excel(writer, sheet_name=sheet_name, index=False)

    print(f"All results saved to {output_excel_path}")
    return output_excel_path


def run_workflow(config, df_original, df_peptide, group_columns, save_excel=True):
    """Run imputation/normalization/FC for each comparison (cell 7).

    Returns ``(imputed_dataframes, imputation_dict, df_final_results)``.
    ``df_original`` should already be a clean copy (this function does not mutate it).
    """
    imputation_dict = {}
    imputed_dataframes = {}
    df_final_results = df_original.copy()  # stores the final summary with all FC columns.

    # Determine the list of comparisons to perform.
    if not config.comparison_matrix:
        comparisons = [(key, config.reference_group) for key in group_columns if key != config.reference_group]
    else:
        comparisons = config.comparison_matrix

    for pair in comparisons:
        treated_name, control_name = pair[0], pair[1]
        comparison_name = f"{treated_name}_vs_{control_name}"
        logger.info("Processing comparison: " + comparison_name)
        # Default to the original data if imputation is turned off.
        df_imputed_for_comparison = df_original.copy()

        if config.imputation_option:
            df_imputed_for_comparison, imputation_list = imputation(
                df_original,
                treated_name,
                control_name,
                group_columns,
                df_peptide,
                mode=config.mode,
            )
            imputation_dict[comparison_name] = imputation_list
            imputed_dataframes[comparison_name] = df_imputed_for_comparison
        else:
            logger.info("Imputation not enabled")
            imputed_dataframes[comparison_name] = df_original

        if config.normalization_protein_id:
            logger.info(f"Normalizing data for comparison '{comparison_name}' using protein '{config.normalization_protein_id}'...")
            df_imputed_for_comparison = normalize_by_specific_protein(
                df_imputed_for_comparison,
                config.normalization_protein_id,
                group_columns,
            )
            imputed_dataframes[comparison_name] = df_imputed_for_comparison
            df_imputed_for_comparison.to_csv(f"Normalized_imputed_data_{comparison_name}.tsv", sep="\t")
        else:
            logger.info("Normalization not enabled")

        # Calculate FC values on a COPY of the imputed data to avoid contamination.
        df_with_fc = calculate_average_FC_value(df_imputed_for_comparison.copy(), treated_name, control_name, group_columns)

        # Extract ONLY the new FC columns and add them to the final summary table.
        fc_col = f'FC_{comparison_name}'
        log2fc_col = f'log2FC_{comparison_name}'
        if fc_col in df_with_fc:
            df_final_results[fc_col] = df_with_fc[fc_col]
            df_final_results[log2fc_col] = df_with_fc[log2fc_col]

    if save_excel:
        _save_workflow_excel(config, df_original, group_columns, imputed_dataframes, df_final_results, imputation_dict)

    logger.info("Summary of imputed proteins per comparison:")
    for key, value in imputation_dict.items():
        logger.info(f"{key}: {len(value)} proteins imputed")

    return imputed_dataframes, imputation_dict, df_final_results


# --------------------------------------------------------------------------- #
# cell 9: limma (moderated t-test; the only statistics path)
# --------------------------------------------------------------------------- #
def run_limma(config, imputed_dataframes, df_final_results, group_columns, save_csv=True):
    # `save_csv` is kept for signature compatibility; the standalone summary CSV is
    # no longer written (the same table is the Excel's Fold_Change_Summary sheet).
    for comparison_name, imputed_df in imputed_dataframes.items():
        print(f"Running Limma for: {comparison_name}...")
        treated_name, control_name = comparison_name.split('_vs_')

        limma_results_df = limma_differential_analysis(
            imputed_df, treated_name, control_name, group_columns
        )

        if not limma_results_df.empty:
            df_final_results = df_final_results.merge(
                limma_results_df, left_index=True, right_index=True, how="left"
            )

    print("\nLimma analysis complete. Final summary has been updated.")
    return df_final_results


# --------------------------------------------------------------------------- #
# Export: significant up-/down-regulated entries across all comparisons
# --------------------------------------------------------------------------- #
def export_significant(result, direction="down", logfc_cutoff=1.0, fdr_cutoff=0.05, filename=None):
    """Write one Excel row per significant (comparison, protein) in one direction.

    direction "down": log2FC <= -logfc_cutoff;  "up": log2FC >= logfc_cutoff.
    Both also require adjusted P <= fdr_cutoff. Columns: comparison, gene, uniprot,
    log2FC, p, adjusted P, imputed (boolean). Returns (filename, DataFrame).
    """
    if direction not in ("down", "up"):
        raise ValueError("direction must be 'down' or 'up'")
    if filename is None:
        filename = "significant_downreg.xlsx" if direction == "down" else "significant_upreg.xlsx"

    summary = result.summary
    rows = []
    for comparison_name in result.imputed_dataframes.keys():
        log2fc_col = f'log2FC_{comparison_name}'
        fdr_col = f'bh_FDR_{comparison_name}'      # adjusted P (BH/FDR)
        p_col = f'Pvalue_{comparison_name}'        # raw moderated-t p
        if log2fc_col not in summary.columns or fdr_col not in summary.columns:
            continue
        imputed_set = set(result.imputation_dict.get(comparison_name, []))
        if direction == "down":
            fc_mask = summary[log2fc_col] <= -logfc_cutoff
        else:
            fc_mask = summary[log2fc_col] >= logfc_cutoff
        mask = fc_mask & (summary[fdr_col] <= fdr_cutoff)
        for idx, row in summary[mask].iterrows():
            rows.append({
                'comparison': comparison_name,
                'gene': row.get('Genes'),
                'uniprot': idx,
                'log2FC': row[log2fc_col],
                'p': row[p_col] if p_col in summary.columns else None,
                'adjusted P': row[fdr_col],
                'imputed': idx in imputed_set,
            })

    out = pd.DataFrame(rows, columns=['comparison', 'gene', 'uniprot', 'log2FC',
                                      'p', 'adjusted P', 'imputed'])
    out.to_excel(filename, index=False)
    logger.info("Wrote %d %s-regulated significant entries to %s.", len(out), direction, filename)
    print(f"{direction.capitalize()}-regulated significant entries ({len(out)}) saved to {filename}")
    return filename, out


def export_downregulated(result, logfc_cutoff=1.0, fdr_cutoff=0.05, filename=None):
    return export_significant(result, "down", logfc_cutoff, fdr_cutoff, filename)


def export_upregulated(result, logfc_cutoff=1.0, fdr_cutoff=0.05, filename=None):
    return export_significant(result, "up", logfc_cutoff, fdr_cutoff, filename)


# --------------------------------------------------------------------------- #
# Convenience: load -> group -> clean-up -> workflow -> stats
# --------------------------------------------------------------------------- #
def run_core(config: Optional[AnalysisConfig] = None, contaminants=None, save_outputs=True) -> AnalysisResult:
    """Run the full core pipeline (load -> group -> clean-up -> impute/FC -> limma)
    and return an AnalysisResult.
    """
    if config is None:
        config = AnalysisConfig()
    if contaminants is None:
        from .reference_data import contaminants as contaminants

    df, df_peptide = load_dataset(config, contaminants)
    group_columns = assign_groups(df, config.group_names)
    df = control_group_cleanup(df, config, group_columns)

    # Defer the Excel export until after limma so the per-comparison sheets can
    # include the p-value / adjusted-P columns (run_workflow runs before limma).
    imputed_dataframes, imputation_dict, summary = run_workflow(
        config, df, df_peptide, group_columns, save_excel=False
    )
    summary = run_limma(config, imputed_dataframes, summary, group_columns, save_csv=save_outputs)
    if save_outputs:
        _save_workflow_excel(config, df, group_columns, imputed_dataframes, summary, imputation_dict)

    result = AnalysisResult(
        config=config,
        df_original=df,
        df_peptide=df_peptide,
        group_columns=group_columns,
        imputed_dataframes=imputed_dataframes,
        imputation_dict=imputation_dict,
        summary=summary,
    )
    if save_outputs:
        export_downregulated(result)
        export_upregulated(result)
    return result
