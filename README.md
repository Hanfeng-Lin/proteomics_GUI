# DIA-NN Proteomics — Volcano Explorer (GUI)

Self-contained desktop app for the DIA-NN proteomics pipeline. Configure an
analysis, run it, and explore the volcano plot with **every** plot setting
exposed. The figure is shown live in the window and also saved to PNG.

## Download

First, get the code from GitHub:

1. Open **https://github.com/Hanfeng-Lin/proteomics_GUI**
2. Click the green **Code** button → **Download ZIP**
3. Unzip it somewhere convenient (e.g. your Desktop)

(Prefer git? `git clone https://github.com/Hanfeng-Lin/proteomics_GUI.git`)

Then follow **Requirements** and **Run** below. To **update later**, double-click
**`update.bat`** — it does a `git pull` if you cloned, otherwise it downloads and
applies the latest ZIP (your data files and `proteomics_GUI_output` are kept).

## Contents

```
proteomics_GUI/
├── start_gui.bat           ← double-click to launch (finds the 'proteomics' env)
├── update.bat              ← double-click to pull the latest version from GitHub
├── VERSION                 ← current version string (used for the update check)
├── gui.py                  ← the application
├── scripts/                ← analysis library the GUI calls
└── environment.yml         ← conda environment spec (recreates the 'proteomics' env)
```

DIA-NN matrices (`<stem>.pg_matrix.tsv` / `<stem>.pr_matrix.tsv`) are **not**
committed to the repo — supply your own via the **Browse...** button in tab 1.

## Requirements

**1. Python environment (conda).** Recreate the env the app was built and tested
in (named `proteomics`):

```bash
conda env create -f environment.yml
conda activate proteomics
```

**2. R + limma (required).** Statistics use R/limma (a moderated t-test), so you
must install R separately — it does **not** come with the conda env:

1. Download and install **R** from CRAN: https://cran.r-project.org/
2. In an R console, install the Bioconductor **limma** package:
   ```r
   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
   BiocManager::install("limma")
   ```
3. Set the **`R_HOME`** environment variable to your R install folder (e.g.
   `C:\Program Files\R\R-4.4.2`) so the app can find R.

Notes:

- Tkinter ships with Python — no install needed.
- `gseapy` (Enrichr/GSEA features) makes network calls when used.

Once the `proteomics` env exists, `start_gui.bat` finds and activates it for you.

## Run

Double-click **`start_gui.bat`**. It checks that Anaconda and the `proteomics`
env are present (and prints clear setup instructions if not), then compares the
local **`VERSION`** with the one on GitHub — if a newer version exists it asks
whether to update before launching. Or, from a shell with the right Python:

```bash
python gui.py
```

The left side is a 4-step notebook; the plot canvas (right) and Log (bottom) are
shared. **Hover over any field, checkbox, or button to see a hint** describing what
it does.

**Tab 1 — Analysis configuration.** Click **Browse...** next to *Protein file*
and pick your `*.pg_matrix.tsv`. The matching `*.pr_matrix.tsv` (same name/folder)
is filled in automatically; Browse for it separately if it differs. The
**working folder** is wherever those files live, shown under the path fields. Set
groups, reference group, and comparisons (`treated:control`, comma-separated —
**blank means every group vs the reference**), adjust the imputation / limma /
Pharos toggles, then click **Run analysis** (runs in the background; progress in
the Log).

Don't want to type the group names? Click **Auto-pick...** next to *Groups*: it
reads the sample-column headers and finds **all** the consensus regions they
share (a name like `…_Target_96Plate_03_…_Tech_…` has several), then offers each
variable field between them as a candidate with a preview of the group names it
would give — handling names that contain `_`, like `Positive_Control`. The most
likely field (usually the compound, first) is preselected; pick one and it fills
the Groups box. The **Group
assignments** box then shows how many samples matched each group — use **Preview
groups** to check this from the file's columns before running.

**Tab 2 — PCA.** Set a title / output filename, whether to label samples, and the
per-element **font sizes**, then click **Plot PCA**. PCA uses all samples with
per-protein mean imputation.

**Tab 3 — Volcano settings.** When the run finishes, a volcano for **every
comparison** is generated automatically, one per tab in the plot area (each saved
to PNG). Adjust any volcano parameter (thresholds, empirical-FDR curve, axes,
highlight sets, **label selection** — up / down / imputed — and **placement**
(adjustText on/off, repel force, arrows), **font sizes**) and click **Plot all
comparisons** to regenerate them, or **Plot selected** for just the
treatment/control pair in the dropdowns.

**Tab 4 — Bubbleplot settings.** Builds the clustered bubble/dendrogram plot over
significantly down-regulated proteins. Define **SAR groups** (one per line,
`label: treatmentA, treatmentB`; prefilled with your treatments after a run) and
the **suffix** appended to each (e.g. `_vs_DMSO`), tweak figure / colour /
highlight / **font-size** options, then click **Plot bubble**. (Needs at least two
proteins with log2FC < -1 and bh_FDR < 0.01 across the chosen treatments.)

**Tab 5 — Raw data lookup.** Enter a **UniProt accession** (e.g. `P51617`) or a
**gene symbol** (e.g. `IRAK1`) and pick a comparison (or *All samples*), then click
**Look up**. The right panel shows a table of the raw (pre-imputation) intensities
for that protein across the selected samples — the protein-group (**PG**) value
plus each **precursor (PR)** value, one column per precursor.

## Where outputs go

All outputs are written to a dedicated **`proteomics_GUI_output/`** folder created
**next to your data** (the working folder). Nothing is written into the
application folder. They include:

- the volcano / PCA PNGs and the main results Excel
  (**`final_analysis_results_imputed.xlsx`**): a `Fold_Change_Summary` sheet (FC,
  log2FC, raw `Pvalue_`, adjusted `bh_FDR_`) plus a per-comparison sheet each with
  the imputed intensities and an **`Imputed`** TRUE/FALSE column,
- **`significant_downreg.xlsx`** and **`significant_upreg.xlsx`** — one row per
  significant entry (log2FC ≤ −1 / ≥ 1, adjusted P ≤ 0.05) with columns
  *comparison, gene, uniprot, log2FC, p, adjusted P, imputed*,
- **`analysis_log.txt`** — the full run log (everything shown in the Log box,
  including which proteins were/weren't imputed and why).

## Using your own data

Just **Browse...** to your DIA-NN `*.pg_matrix.tsv` anywhere on disk — the app
works from that file's folder, so you do not need to copy data into this folder.
Keep the DIA-NN naming (`<stem>.pg_matrix.tsv` and `<stem>.pr_matrix.tsv` in the
same folder) so the precursor file auto-fills, and set **Groups**, **Reference
group**, and **Comparisons** to match your sample-column names.

## Notes

- The analysis math is identical to the original notebook
  (`DataAnalysis_Jurkat.ipynb`); the GUI only collects parameters and calls the
  library. See `scripts/README.md` for the library design.
