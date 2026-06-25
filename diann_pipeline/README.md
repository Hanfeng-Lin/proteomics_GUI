# diann_pipeline

Importable-library version of `DataAnalysis_Jurkat.ipynb`. The notebook's single
flat global namespace is replaced by explicit objects and functions, but **the
analysis math is unchanged** — verified to reproduce the notebook's statistics to
~1e-15 (float round-off) on the Jurkat test data.

## Design

The notebook kept everything in module globals (`df`, `group_columns`,
`imputation_dict`, the loop variable `pair`, config flags, the reference lists).
Each function now takes what it needs as arguments, and two objects carry state:

- **`AnalysisConfig`** (`config.py`) — all settings (was cell 0). Defaults
  reproduce the original Jurkat run; construct a new one per dataset.
- **`AnalysisResult`** (`pipeline.py`) — dataset + `group_columns` +
  `imputed_dataframes` + `imputation_dict` + `summary` (the old `df_final_results`).

| Module | Notebook origin | Responsibility |
|---|---|---|
| `config.py` | cell 0 | `AnalysisConfig` dataclass |
| `reference_data.py` | cell 1 | contaminants, Pharos, kinase/ub, loop lists (verbatim) |
| `io.py` | cells 2,3,5 | load + decontaminate, group assignment, control clean-up |
| `imputation.py` | cell 6 | imputation, normalization, fold change |
| `stats.py` | cells 8,9,10 | t-test + BH, limma (lazy R), empirical-FDR curve |
| `pipeline.py` | cells 7,8,9 | `run_workflow` / `run_ttest` / `run_limma` / `run_core` |
| `plots/pca.py` | cell 4 | `generate_pca_plot` |
| `plots/volcano.py` | cell 10 | `volcano_plot` |
| `plots/bubble.py` | cell 12 | `bubble_dendro_plot` |
| `plots/gsea.py` | cell 16 | `gsea_analysis` |

What changed vs. the notebook (plumbing only, no math):
- Functions take `group_columns`, `df`, `config`, etc. as arguments instead of globals.
- The volcano imputation-highlight block, which read the loop var `pair`, was
  rewritten to its identical single form `imputation_dict[treatment+"_vs_"+control]`.
- R/rpy2 for limma is imported lazily, so importing the package doesn't need R.
- `print` → `logging` in the data-prep stages.
- **No random seed was added** — imputation is still stochastic, exactly as before.

## Usage

Run from the directory holding `<file>.pg_matrix.tsv` / `<file>.pr_matrix.tsv`,
with the parent of `diann_pipeline/` on `PYTHONPATH` (or run the example, which is
already there):

```python
from diann_pipeline import AnalysisConfig, run_core
from diann_pipeline.plots import generate_pca_plot, volcano_plot

cfg = AnalysisConfig()                    # or AnalysisConfig(file="myrun", group_names=[...], ...)
result = run_core(cfg)                    # load -> impute -> FC -> (t-test) -> limma

generate_pca_plot(result.df_original, result.group_columns, text=True)
for treated, control in cfg.comparison_matrix:
    volcano_plot(treated, control,
                 df=result.summary, group_columns=result.group_columns,
                 imputation_dict=result.imputation_dict, config=cfg,
                 FDR_cutoff=0.05, logFC_cutoff=1, xlim=[-10, 10])
```

`run_jurkat.py` (in the project root) is a complete runnable example.

### Analysing a different dataset

```python
cfg = AnalysisConfig(
    file="my_experiment",                 # my_experiment.pg_matrix.tsv / .pr_matrix.tsv
    group_names=["Ctrl", "DrugA", "DrugB"],
    reference_group="Ctrl",
    comparison_matrix=(["DrugA", "Ctrl"], ["DrugB", "Ctrl"]),
    mode=0, imputation_option=True, limma_option=True,
)
result = run_core(cfg)
```

### Optional stages

- **Bubble plot** (`plots.bubble_dendro_plot`) reads `<file>_analyzed.csv`. Write it
  first: `result.summary.to_csv(cfg.file.split(".")[0] + "_analyzed.csv")`, then pass
  an `SAR` dict whose entries match your treatment-group suffixes.
- **GSEA** (`plots.gsea_analysis`) makes network calls to Enrichr.

## Requirements

See the project `environment.yml` (conda). limma needs R + `R_HOME` + the Bioconductor
`limma` package; gseapy needs network access. Both are only imported when used.
