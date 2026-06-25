"""Analysis configuration.

Replaces the notebook's module-level settings (cell 0) with an explicit,
passable object. Field values default to the original Jurkat notebook settings,
so ``AnalysisConfig()`` reproduces that run. Construct a new instance (or use
``dataclasses.replace``) to analyse a different dataset without touching code.

Field name -> original notebook global:
    mode                            -> MODE
    file                            -> FILE
    group_names                     -> GROUP_NAME
    reference_group                 -> reference_group
    comparison_matrix               -> comparison_matrix
    control_group_detection_threshold -> control_group_detection_threshold
    imputation_option               -> imputation_option
    normalization_protein_id        -> NORMALIZATION_PROTEIN_ID
    pharos_tcrd                     -> PharosTCRD
    limma_option                    -> limma_option
    output_adjpval                  -> output_adjPval
"""

from dataclasses import dataclass, field
from typing import List, Optional, Sequence


@dataclass
class AnalysisConfig:
    # 0 for "global proteomics", 1 for "positive enrichment pulldown, e.g. biotin"
    mode: int = 0
    # FILENAME stem of your <file>.pg_matrix.tsv / <file>.pr_matrix.tsv. No extension.
    file: str = "diann"
    group_names: List[str] = field(
        default_factory=lambda: ["DMSO", "Positive_Control", "IRAK1"]
    )
    # Reference group: should have the least missing values (most complete proteome).
    reference_group: str = "DMSO"
    # Each entry is [treatment_group, reference_group]. Empty => all groups vs reference.
    comparison_matrix: Sequence = field(
        default_factory=lambda: (["Positive_Control", "DMSO"], ["IRAK1", "DMSO"])
    )
    # 0.5 => keep proteins seen in >= half of control samples. 0 => discard none.
    control_group_detection_threshold: float = 0.5
    # Enable imputation on treatment and reference groups.
    imputation_option: bool = True
    # e.g. "Q13085" to normalise to one protein's abundance. "" => no normalisation.
    normalization_protein_id: str = ""
    # Colour the proteome by Pharos TCRD class (Tclin/Tchem/Tbio/Tdark).
    pharos_tcrd: bool = False
    # Use R-limma for statistics (needs R + R_HOME + limma). False => Student t-test + BH.
    limma_option: bool = True
    # Output adjusted p-value (FDR) instead of raw p-value.
    output_adjpval: bool = True
    # Optional explicit input paths. If set, they override the `file` stem when
    # loading -- lets a caller (e.g. the GUI file browser) load matrices from any
    # folder regardless of the stem. Leave None to use `<file>.pg_matrix.tsv` etc.
    pg_path: Optional[str] = None
    pr_path: Optional[str] = None
