"""diann_pipeline -- importable library version of DataAnalysis_Jurkat.ipynb.

The notebook's flat global namespace is replaced by explicit objects:
    * AnalysisConfig    -- all settings (was cell 0)
    * AnalysisResult    -- dataset + per-comparison results (was the loose globals)
The analysis math is unchanged from the notebook.

Typical use::

    from diann_pipeline import AnalysisConfig, run_core
    from diann_pipeline.plots import generate_pca_plot, volcano_plot

    cfg = AnalysisConfig()                 # defaults reproduce the Jurkat run
    result = run_core(cfg)                 # load -> impute -> FC -> limma
    volcano_plot("IRAK1", "DMSO",
                 df=result.summary, group_columns=result.group_columns,
                 imputation_dict=result.imputation_dict, config=cfg,
                 FDR_cutoff=0.05, logFC_cutoff=1, xlim=[-10, 10])

See ``run_jurkat.py`` for an end-to-end example.
"""

from .config import AnalysisConfig
from .io import load_dataset, assign_groups, control_group_cleanup
from .imputation import (
    imputation,
    normalize_by_specific_protein,
    calculate_average_FC_value,
)
from .stats import (
    limma_differential_analysis,
    empirical_fdr_curve,
)
from .pipeline import (
    AnalysisResult,
    run_workflow,
    run_limma,
    run_core,
    export_downregulated,
)

__all__ = [
    "AnalysisConfig",
    "AnalysisResult",
    "load_dataset",
    "assign_groups",
    "control_group_cleanup",
    "imputation",
    "normalize_by_specific_protein",
    "calculate_average_FC_value",
    "limma_differential_analysis",
    "empirical_fdr_curve",
    "run_workflow",
    "run_limma",
    "run_core",
    "export_downregulated",
]
