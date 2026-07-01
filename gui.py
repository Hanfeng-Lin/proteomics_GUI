"""Desktop GUI for the proteomics analysis (the `scripts` package).

Left side is a 3-step notebook:
  1. Analysis configuration  -> browse pg/pr files, set groups, run run_core()
     in a background thread (load -> impute -> fold change -> t-test/limma).
  2. Volcano settings        -> every volcano_plot() parameter; plots all
     comparisons (one tab each) or a single selected one.
  3. Bubbleplot settings     -> every bubble_dendro_plot() parameter.

The plot canvas (right) and log (bottom) are shared across the steps. The working
folder is wherever the chosen pg/pr files live, and all outputs are written to a
dedicated "<stem>_outputs" folder beside the data.

The GUI only collects parameters and calls the functions in `scripts`.
Run it with:  python gui.py
"""

import os
import re
import sys
import glob
import json
import queue
import logging
import threading
import traceback
import contextvars

import tkinter as tk
from tkinter import ttk, messagebox, filedialog, colorchooser

import numpy as np
import matplotlib
matplotlib.use("Agg")  # figures are embedded manually; no stray windows
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

try:
    from PIL import Image, ImageTk
    _HAS_PIL = True
except Exception:
    _HAS_PIL = False

# Make the package importable when running this file directly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scripts import AnalysisConfig, run_core
from scripts.plots.volcano import volcano_plot
from scripts.plots.pca import generate_pca_plot
from scripts.plots.bubble import bubble_dendro_plot


def _app_version():
    """The version string from the VERSION file next to this script, or ''."""
    try:
        here = os.path.dirname(os.path.abspath(__file__))
        with open(os.path.join(here, "VERSION"), encoding="utf-8") as f:
            return f.read().strip()
    except Exception:
        return ""


def _version_suffix():
    v = _app_version()
    return f"  (v{v})" if v else ""


# --------------------------------------------------------------------------- #
# Parsing helpers (pure functions -- testable without a display)
# --------------------------------------------------------------------------- #
def _opt_float(s):
    s = str(s).strip()
    return None if s == "" else float(s)


def _opt_int(s):
    s = str(s).strip()
    return None if s == "" else int(float(s))


def _req_float(s, default):
    s = str(s).strip()
    return default if s == "" else float(s)


def _req_int(s, default):
    s = str(s).strip()
    return default if s == "" else int(float(s))


def _as_float(s, name):
    s = str(s).strip()
    if s == "":
        raise ValueError(f"'{name}' is required")
    return float(s)


def _genes(s):
    return [x.strip() for x in str(s).split(",") if x.strip()]


def _num_pair(s, name):
    parts = [p.strip() for p in str(s).split(",") if p.strip() != ""]
    if len(parts) != 2:
        raise ValueError(f"'{name}' must be two numbers like '1, 3'")
    return [float(parts[0]), float(parts[1])]


def _force(s, default):
    """Parse an adjustText force: '' -> default, '1' -> 1.0, '1, 2' -> (1.0, 2.0)."""
    s = str(s).strip()
    if s == "":
        return default
    parts = [p.strip() for p in s.split(",") if p.strip() != ""]
    if len(parts) == 1:
        return float(parts[0])
    return (float(parts[0]), float(parts[1]))


def build_volcano_params(v):
    """Turn a dict of raw widget values into volcano_plot() keyword arguments."""
    xmin = _as_float(v["xlim_min"], "x-axis min")
    xmax = _as_float(v["xlim_max"], "x-axis max")

    ymin = str(v["ylim_min"]).strip()
    ymax = str(v["ylim_max"]).strip()
    if ymin == "" and ymax == "":
        ylim = []  # auto y-limits
    else:
        ylim = [_as_float(ymin, "y-axis min"), _as_float(ymax, "y-axis max")]

    mode_v = str(v["mode"]).strip().lower()
    mode = None if mode_v in ("", "auto") else int(float(mode_v))

    return dict(
        logFC_cutoff=_opt_float(v["logFC_cutoff"]),
        logFC_cutoff2=_opt_float(v["logFC_cutoff2"]),
        FDR_cutoff=_req_float(v["FDR_cutoff"], 0.05),
        use_empirical_fdr=bool(v["use_empirical_fdr"]),
        mode=mode,
        fdr_alpha=_req_float(v["fdr_alpha"], 0.05),
        kappa=_req_float(v["kappa"], 1e-6),
        p_value_cutoff=_req_float(v["p_value_cutoff"], 1),
        file_suffix=str(v["file_suffix"]),
        highlight_genes=_genes(v["highlight_genes"]),
        protein_level_cutoff=_opt_float(v["protein_level_cutoff"]),
        xlim=[xmin, xmax],
        ylim=ylim,
        x_interval=_req_float(v["x_interval"], 2),
        y_interval=_req_float(v["y_interval"], 1),
        top_buffer=_req_float(v["top_buffer"], 0.1),
        imputation_option=bool(v["imputation_option"]),
        PharosTCRD=bool(v["PharosTCRD"]),
        highlight_kinase=bool(v["highlight_kinase"]),
        highlight_ub=bool(v["highlight_ub"]),
        highlight_Gloops=bool(v["highlight_Gloops"]),
        highlight_RTloops=bool(v["highlight_RTloops"]),
        label_topX_mid_fc=_opt_int(v["label_topX_mid_fc"]),
        max_label=_req_int(v["max_label"], 100),
        label_most_extreme=_opt_int(v["label_most_extreme"]),
        label_up=bool(v["label_up"]),
        label_down=bool(v["label_down"]),
        label_imputed=bool(v["label_imputed"]),
        adjust_labels=bool(v["adjust_labels"]),
        adjust_arrows=bool(v["adjust_arrows"]),
        adjust_force_text=_force(v["adjust_force_text"], (1, 2)),
        adjust_force_static=_force(v["adjust_force_static"], (1, 2)),
        dpi=_req_int(v["dpi"], 300),
        dot_size=_req_float(v["dot_size"], 40),
        dot_alpha=_req_float(v["dot_alpha"], 0.5),
        color_bg=str(v["color_bg"]),
        color_up=str(v["color_up"]),
        color_down=str(v["color_down"]),
        color_imputed=str(v["color_imputed"]),
        color_highlight=str(v["color_highlight"]),
        color_kinase=str(v["color_kinase"]),
        color_ub=str(v["color_ub"]),
        color_gloops=str(v["color_gloops"]),
        color_rtloops=str(v["color_rtloops"]),
        color_tclin=str(v["color_tclin"]),
        color_tchem=str(v["color_tchem"]),
        color_tbio=str(v["color_tbio"]),
        color_tdark=str(v["color_tdark"]),
        title_fontsize=_req_float(v["title_fontsize"], 24),
        axis_label_fontsize=_req_float(v["axis_label_fontsize"], 20),
        tick_fontsize=_req_float(v["tick_fontsize"], 16),
        legend_fontsize=_req_float(v["legend_fontsize"], 12),
        gene_label_fontsize=_req_float(v["gene_label_fontsize"], 14),
    )


def build_bubble_params(v):
    """Turn raw widget values into (SAR dict, bubble_dendro_plot kwargs)."""
    sar_text = str(v["sar"]).strip()
    if not sar_text:
        raise ValueError("Enter at least one SAR group, e.g.:  : Positive_Control, IRAK1")
    SAR = {}
    for line in sar_text.splitlines():
        line = line.strip()
        if not line:
            continue
        if ":" in line:
            label, items = line.split(":", 1)
        else:
            label, items = "", line
        names = [x.strip() for x in items.split(",") if x.strip()]
        if names:
            SAR[label.strip()] = names
    if not SAR:
        raise ValueError("No valid SAR groups parsed.")

    legend_raw = str(v["legend_num"]).strip()
    legend_num = "auto" if legend_raw == "" or legend_raw.lower() == "auto" else int(float(legend_raw))

    kwargs = dict(
        SAR_suffix=str(v["sar_suffix"]),
        figure_filename=str(v["figure_filename"]).strip() or "bubble_plot.png",
        fig_title=str(v["fig_title"]),
        fig_width=_req_float(v["fig_width"], 50),
        fig_height=_req_float(v["fig_height"], 50),
        dendro_bubble_height_ratio=_num_pair(v["dendro_bubble_height_ratio"], "dendro/bubble height ratio"),
        bubble_legend_width_ratio=_num_pair(v["bubble_legend_width_ratio"], "bubble/legend width ratio"),
        compound_labelsize=_req_float(v["compound_labelsize"], 10),
        protein_labelsize=_req_float(v["protein_labelsize"], 10),
        colorFCrange=_num_pair(v["colorFCrange"], "color FC range"),
        highlight_G_loop=int(float(str(v["highlight_G_loop"]).strip() or 0)),
        highlight_RT_loop=int(float(str(v["highlight_RT_loop"]).strip() or 0)),
        rainbow_palette=1 if bool(v["rainbow_palette"]) else 0,
        invert_xy=bool(v["invert_xy"]),
        selected_genes=_genes(v["selected_genes"]),
        legend_num=legend_num,
        protein_label="gene" if str(v.get("protein_label_mode", "")).startswith("Gene") else "description_gene",
        bubble_size_scale=_req_float(v["bubble_size_scale"], 1.0),
        title_fontsize=_req_float(v["title_fontsize"], 14),
        axis_fontsize=_req_float(v["axis_fontsize"], 12),
        colorbar_label_fontsize=_req_float(v["colorbar_label_fontsize"], 12),
        colorbar_tick_fontsize=_req_float(v["colorbar_tick_fontsize"], 10),
        legend_fontsize=_req_float(v["legend_fontsize"], 10),
        dpi=_req_int(v["dpi"], 200),
    )
    return SAR, kwargs


def _stem_from_pg(path):
    """Derive the DIA-NN run stem from a *.pg_matrix.tsv path (or any file)."""
    name = os.path.basename(path)
    for suffix in (".pg_matrix.tsv", ".pr_matrix.tsv", ".tsv"):
        if name.lower().endswith(suffix):
            return name[: -len(suffix)]
    return os.path.splitext(name)[0]


def build_config(v):
    """Build an AnalysisConfig from the config-tab widget values.

    Inputs are the browsed pg/pr file paths; the working folder is wherever the
    pg file lives. Raises ValueError (with a human message) on bad input.
    """
    groups = [g.strip() for g in re.split(r"[,\n]", str(v["group_names"])) if g.strip()]
    # Optional exact-match patterns from the auto-picker; keep only those whose
    # label is still among the current group names (a manual edit can drop some).
    gp_raw = v.get("group_patterns") or {}
    group_patterns = {g: gp_raw[g] for g in groups if g in gp_raw} or None

    comp_raw = str(v["comparison_matrix"]).strip()
    if comp_raw == "":
        comparison_matrix = ()  # blank => all groups vs the reference group
    else:
        comparison_matrix = []
        for pair in comp_raw.split(","):
            if ":" not in pair:
                raise ValueError(f"Comparison '{pair.strip()}' must be 'treated:control'")
            t, c = pair.split(":", 1)
            comparison_matrix.append([t.strip(), c.strip()])
        comparison_matrix = tuple(comparison_matrix)

    pg = str(v.get("pg_path", "")).strip()
    pr = str(v.get("pr_path", "")).strip()
    if not pg:
        raise ValueError("Select a protein (pg_matrix.tsv) file with 'Browse...'.")
    if not pr:
        raise ValueError("Select a precursor (pr_matrix.tsv) file with 'Browse...'.")
    pg = os.path.abspath(pg)
    pr = os.path.abspath(pr)
    if not os.path.exists(pg):
        raise ValueError(f"pg_matrix file not found:\n{pg}")
    if not os.path.exists(pr):
        raise ValueError(f"pr_matrix file not found:\n{pr}")

    return AnalysisConfig(
        mode=int(float(str(v["mode"]).strip() or 0)),
        file=_stem_from_pg(pg),
        group_names=groups,
        reference_group=str(v["reference_group"]).strip(),
        comparison_matrix=comparison_matrix,
        control_group_detection_threshold=float(str(v["control_threshold"]).strip() or 0.5),
        imputation_option=bool(v["imputation_option"]),
        normalization_protein_id=str(v["normalization_protein_id"]).strip(),
        pharos_tcrd=False,  # volcano-only display option; controlled on the Volcano tab
        output_adjpval=bool(v["output_adjpval"]),
        pg_path=pg,
        pr_path=pr,
        group_patterns=group_patterns,
        drop_samples=[s for s in (v.get("drop_samples") or []) if str(s).strip()],
    )


# --------------------------------------------------------------------------- #
# Log plumbing
# --------------------------------------------------------------------------- #
class _QueueWriter:
    """File-like object that funnels writes into a queue (for stdout capture)."""
    def __init__(self, q):
        self.q = q

    def write(self, s):
        if s:
            self.q.put(s)

    def flush(self):
        pass


class _QueueLogHandler(logging.Handler):
    def __init__(self, q):
        super().__init__()
        self.q = q

    def emit(self, record):
        try:
            self.q.put(self.format(record) + "\n")
        except Exception:
            pass


# --------------------------------------------------------------------------- #
# Small widget helpers
# --------------------------------------------------------------------------- #
class Tooltip:
    """Show a small hint popup when the mouse hovers over a widget."""
    def __init__(self, widget, text, delay=450, wrap=340):
        self.widget = widget
        self.text = text
        self.delay = delay
        self.wrap = wrap
        self.tip = None
        self._after = None
        widget.bind("<Enter>", self._schedule, add="+")
        widget.bind("<Leave>", self._hide, add="+")
        widget.bind("<ButtonPress>", self._hide, add="+")

    def _schedule(self, _event=None):
        self._cancel()
        self._after = self.widget.after(self.delay, self._show)

    def _cancel(self):
        if self._after is not None:
            try:
                self.widget.after_cancel(self._after)
            except Exception:
                pass
            self._after = None

    def _show(self):
        if self.tip is not None:
            return
        try:
            x = self.widget.winfo_rootx() + 18
            y = self.widget.winfo_rooty() + self.widget.winfo_height() + 4
        except Exception:
            return
        self.tip = tk.Toplevel(self.widget)
        self.tip.wm_overrideredirect(True)
        self.tip.wm_geometry(f"+{x}+{y}")
        tk.Label(self.tip, text=self.text, justify="left", background="#ffffe0",
                 relief="solid", borderwidth=1, wraplength=self.wrap, padx=6, pady=4).pack()

    def _hide(self, _event=None):
        self._cancel()
        if self.tip is not None:
            self.tip.destroy()
            self.tip = None


def labeled_entry(parent, row, label, default="", width=14, tip=None, hint=None):
    lbl = ttk.Label(parent, text=label)
    lbl.grid(row=row, column=0, sticky="w", padx=4, pady=2)
    var = tk.StringVar(value=str(default))
    ent = ttk.Entry(parent, textvariable=var, width=width)
    ent.grid(row=row, column=1, sticky="w", padx=4, pady=2)
    if tip:
        ttk.Label(parent, text=tip, foreground="#666").grid(row=row, column=2, sticky="w", padx=4)
    if hint:
        Tooltip(lbl, hint)
        Tooltip(ent, hint)
    return var


def check(parent, row, label, default=False, hint=None):
    var = tk.BooleanVar(value=default)
    cb = ttk.Checkbutton(parent, text=label, variable=var)
    cb.grid(row=row, column=0, columnspan=3, sticky="w", padx=4, pady=1)
    if hint:
        Tooltip(cb, hint)
    return var


def _readable_fg(color):
    """Black or white text depending on how light the background colour is."""
    try:
        import tkinter as _tk
        # Use a hidden widget to resolve names to RGB if needed.
        r, g, b = (0, 0, 0)
        c = str(color)
        if c.startswith("#") and len(c) == 7:
            r, g, b = int(c[1:3], 16), int(c[3:5], 16), int(c[5:7], 16)
        else:
            named = {"black": (0, 0, 0), "white": (255, 255, 255), "grey": (128, 128, 128),
                     "gray": (128, 128, 128), "red": (255, 0, 0), "blue": (0, 0, 255),
                     "green": (0, 128, 0), "orange": (255, 165, 0), "yellow": (255, 255, 0),
                     "gold": (255, 215, 0), "cyan": (0, 255, 255), "magenta": (255, 0, 255)}
            r, g, b = named.get(c.lower(), (128, 128, 128))
        return "black" if (0.299 * r + 0.587 * g + 0.114 * b) > 150 else "white"
    except Exception:
        return "black"


def color_row(parent, row, label, default, hint=None):
    """A label + a swatch button that opens a colour picker; returns the colour var."""
    lbl = ttk.Label(parent, text=label)
    lbl.grid(row=row, column=0, sticky="w", padx=4, pady=2)
    var = tk.StringVar(value=str(default))
    btn = tk.Button(parent, textvariable=var, width=12, relief="groove")

    def _restyle(*_a):
        c = var.get()
        try:
            btn.configure(background=c, activebackground=c, foreground=_readable_fg(c))
        except Exception:
            pass

    def _pick():
        _rgb, hx = colorchooser.askcolor(color=var.get() or None, title=label, parent=parent)
        if hx:
            var.set(hx)
            _restyle()

    btn.configure(command=_pick)
    btn.grid(row=row, column=1, sticky="w", padx=4, pady=2)
    _restyle()
    if hint:
        Tooltip(lbl, hint)
        Tooltip(btn, hint)
    return var


def labeled_combo(parent, row, label, values, default, width=8, hint=None):
    lbl = ttk.Label(parent, text=label)
    lbl.grid(row=row, column=0, sticky="w", padx=4, pady=2)
    var = tk.StringVar(value=str(default))
    cb = ttk.Combobox(parent, textvariable=var, values=[str(v) for v in values],
                      state="readonly", width=width)
    cb.grid(row=row, column=1, sticky="w", padx=4, pady=2)
    if hint:
        Tooltip(lbl, hint)
        Tooltip(cb, hint)
    return var


def labeled_slider(parent, row, label, frm, to, default, resolution=0.1, length=160, hint=None):
    """A label + a horizontal slider that shows its current value; returns the DoubleVar."""
    lbl = ttk.Label(parent, text=label)
    lbl.grid(row=row, column=0, sticky="w", padx=4, pady=2)
    var = tk.DoubleVar(value=default)
    sld = tk.Scale(parent, from_=frm, to=to, resolution=resolution, orient="horizontal",
                   variable=var, showvalue=True, length=length)
    sld.grid(row=row, column=1, sticky="w", padx=4, pady=2)
    if hint:
        Tooltip(lbl, hint)
        Tooltip(sld, hint)
    return var


def path_row(parent, row, label, browse_cmd, width=32, hint=None):
    """A row with a label, a path entry, and a Browse button."""
    lbl = ttk.Label(parent, text=label)
    lbl.grid(row=row, column=0, sticky="w", padx=4, pady=2)
    var = tk.StringVar()
    ent = ttk.Entry(parent, textvariable=var, width=width)
    ent.grid(row=row, column=1, sticky="w", padx=4, pady=2)
    ttk.Button(parent, text="Browse...", command=browse_cmd).grid(row=row, column=2, sticky="w", padx=4)
    if hint:
        Tooltip(lbl, hint)
        Tooltip(ent, hint)
    return var


_DELIMS = set("_-./ :|\\\t")


def _uniq(vals):
    """Order-preserving unique, trimming whitespace and dropping empties."""
    seen, out = set(), []
    for v in vals:
        v = v.strip()
        if v and v not in seen:
            seen.add(v)
            out.append(v)
    return out


def _longest_common_substring(strings):
    """Longest substring (>= 2 chars) present in every string; '' if none."""
    if not strings:
        return ""
    shortest = min(strings, key=len)
    n = len(shortest)
    for length in range(n, 1, -1):
        for start in range(0, n - length + 1):
            sub = shortest[start:start + length]
            if all(sub in s for s in strings):
                return sub
    return ""


def _decompose(cores, depth=0):
    """Recursively split aligned strings on their consensus substrings.

    Returns an ordered list of tokens: ('anchor', str) for each constant region
    shared by all samples, and ('seg', [values]) for the variable segments in
    between. Handles MULTIPLE consensus regions (e.g. '..._Target_..._Tech_...'),
    so every variable field becomes its own segment.
    """
    anchor = _longest_common_substring(cores) if depth < 30 else ""
    if not anchor or len(anchor) < 2:
        return [("seg", list(cores))]
    left, right = [], []
    for c in cores:
        i = c.find(anchor)
        left.append(c[:i])
        right.append(c[i + len(anchor):])
    return _decompose(left, depth + 1) + [("anchor", anchor)] + _decompose(right, depth + 1)


def _make_part(values):
    """A decomposition slot: 'const' if every sample shares the value, else 'var'
    (a variable field the user can pick). ``values`` is aligned per sample."""
    uniq = _uniq(values)
    if len(uniq) <= 1:
        return {"kind": "const", "text": values[0] if values else ""}
    return {"kind": "var", "values": list(values), "uniq": uniq, "label": ""}


def _split_fields(values, delim):
    """Split each value on ``delim`` into positional sub-columns. Only fires when
    EVERY value has at least 2 parts -- this keeps multi-token names like
    'Positive_Control' whole (DMSO has just 1 part), while still splitting regular
    'NR-TDxdR_TDXd_1' names. Returns ``min(parts)`` columns; the last column
    absorbs any extra trailing tokens, so a ragged tail (an optional re-injection
    timestamp, or controls that drop a field) doesn't break the alignment. Returns
    None when there's nothing safe to split."""
    if not delim:
        return None
    pieces = [v.split(delim) for v in values]
    k = min(len(p) for p in pieces)
    if k < 2:
        return None
    cols = []
    for j in range(k):
        if j < k - 1:
            cols.append([p[j] for p in pieces])
        else:                                   # last kept field swallows the tail
            cols.append([delim.join(p[j:]) for p in pieces])
    return cols


def _label_fields(parts):
    """Give each variable field a human description from its neighbouring anchors."""
    n = 0
    for i, p in enumerate(parts):
        if p["kind"] != "var":
            continue
        n += 1
        def _anchor(rng):
            return next((parts[j]["text"].strip("".join(_DELIMS))
                         for j in rng if parts[j]["kind"] == "const"
                         and parts[j]["text"].strip("".join(_DELIMS))), "")
        nxt = _anchor(range(i + 1, len(parts)))
        prv = _anchor(range(i - 1, -1, -1))
        if nxt:
            p["label"] = f"Field {n} (before '{nxt}')"
        elif prv:
            p["label"] = f"Field {n} (after '{prv}')"
        else:
            p["label"] = f"Field {n}"


def analyze_sample_fields(sample_cols, delim="_", rep_suffix=""):
    """Decompose aligned sample names into ordered parts so the caller can pick
    which variable fields define a group.

    Returns ``(prefix, suffix, parts)`` where ``parts`` is an ordered list of:
      - ``{"kind": "const", "text": str}``   -- a region shared by every sample
      - ``{"kind": "var", "values": [...per sample...], "uniq": [...], "label": str}``
    For every sample i, ``column == prefix + concat(const text | values[i]) + REP + suffix``
    where REP is the (per-sample, possibly empty) replicate suffix matched by
    ``rep_suffix``.

    Three steps make this versatile:
      1. trim the shared prefix/suffix;
      2. if ``rep_suffix`` (a regex) is given, strip the trailing replicate run
         from each core first -- this collapses ragged tails like an optional
         re-injection timestamp AND controls that drop a field, so e.g. CDK2's
         'compound_dose_rep[_timestamp]' and 'DMSO_rep' both reduce to the group
         body ('22-11_1000nM' / 'DMSO');
      3. a consensus pass (longest-common-substring) splits on shared regions
         (keeps 'Positive_Control' whole, isolates a middle plate-well field),
         then any leftover region where every sample has >= 2 ``delim`` parts is
         split positionally (resistance / treatment fields).
    """
    cols = list(sample_cols)
    prefix = os.path.commonprefix(cols)
    suffix = os.path.commonprefix([c[::-1] for c in cols])[::-1]
    # Trim the common prefix/suffix back to a delimiter boundary so we never eat
    # into a shared name token (e.g. 'DrugA'/'DrugB' -> keep 'Drug', strip './').
    last_delim = max([i for i, ch in enumerate(prefix) if ch in _DELIMS], default=-1)
    prefix = prefix[:last_delim + 1]
    first_delim = next((i for i, ch in enumerate(suffix) if ch in _DELIMS), None)
    suffix = suffix[first_delim:] if first_delim is not None else ""
    cores = [c[len(prefix): len(c) - len(suffix)] if suffix else c[len(prefix):] for c in cols]

    if rep_suffix:
        try:
            rep_re = re.compile(rep_suffix + r"$")
            cores = [(lambda m, c: c[:m.start()] if m and m.start() else c)(rep_re.search(c), c)
                     for c in cores]
        except re.error:
            pass

    parts = []
    for kind, payload in _decompose(cores):
        if kind == "anchor":
            parts.append({"kind": "const", "text": payload})
            continue
        sub = _split_fields(payload, delim)         # payload is the per-sample list
        if sub is None:
            parts.append(_make_part(payload))
        else:
            for j, vals in enumerate(sub):
                if j:
                    parts.append({"kind": "const", "text": delim})
                parts.append(_make_part(vals))
    _label_fields(parts)
    return prefix, suffix, parts


def build_field_groups(prefix, suffix, parts, selected, delim="_", rep_suffix=""):
    """Given the variable-field indices the user selected as the group key, return
    an ordered list of ``(label, regex)``. The label joins the selected field
    values; the regex is anchored full-match with unselected fields wildcarded and
    the replicate suffix made optional, so matching is exact (no substring overlap
    between e.g. 'NR' and 'NR-TDxdR') while every replicate still matches."""
    sel = [i for i in selected if parts[i]["kind"] == "var"]
    var_vals = [p["values"] for p in parts if p["kind"] == "var"]
    n = len(var_vals[0]) if var_vals else 0
    sel_set = set(sel)
    tail = f"(?:{rep_suffix})?" if rep_suffix else ""
    out, seen = [], set()
    for s in range(n):
        key = tuple(parts[i]["values"][s] for i in sel)
        if key in seen:
            continue
        seen.add(key)
        label = delim.join(key) if key else "all"
        pat = ["^", re.escape(prefix)]
        for i, p in enumerate(parts):
            if p["kind"] == "const":
                pat.append(re.escape(p["text"]))
            elif i in sel_set:
                pat.append(re.escape(p["values"][s]))
            else:
                pat.append(".*?")
        pat.append(tail + re.escape(suffix) + "$")
        out.append((label, "".join(pat)))
    return out


class GroupPickerDialog(tk.Toplevel):
    """Popup to pick which sample-name field(s) define the groups. Multi-select:
    every checked field becomes part of the group key, everything else is treated
    as a replicate. Emits clean labels plus exact, anchored match patterns."""
    def __init__(self, parent, sample_cols, apply_cb):
        super().__init__(parent)
        self.title("Auto-pick group names")
        self.geometry("720x540")
        self.minsize(620, 420)
        self.apply_cb = apply_cb
        self.sample_cols = list(sample_cols)
        self.transient(parent)
        self.delim = tk.StringVar(value="_")
        # Trailing replicate run to strip before grouping (regex). The default
        # eats one or more '_<digits>' groups, which folds replicate numbers AND
        # an optional re-injection timestamp into the replicate dimension, so the
        # remaining body (e.g. '22-11_1000nM', 'DMSO') becomes the group. Clear it
        # to keep every field.
        self.rep = tk.StringVar(value=r"(_\d+)+")

        top = ttk.Frame(self)
        top.pack(side="top", fill="x", padx=10, pady=(8, 2))
        self.info = ttk.Label(top, justify="left", foreground="#444")
        self.info.pack(anchor="w")
        drow = ttk.Frame(self)
        drow.pack(side="top", fill="x", padx=10, pady=2)
        ttk.Label(drow, text="Field delimiter:").pack(side="left")
        ent = ttk.Entry(drow, textvariable=self.delim, width=4)
        ent.pack(side="left", padx=(4, 8))
        ttk.Label(drow, text="Replicate suffix to ignore (regex):").pack(side="left")
        rent = ttk.Entry(drow, textvariable=self.rep, width=12)
        rent.pack(side="left", padx=(4, 0))
        self.delim.trace_add("write", lambda *_: self._rebuild())
        self.rep.trace_add("write", lambda *_: self._rebuild())
        ttk.Label(self, text="Tick every field that defines a group (others = replicates):",
                  foreground="#444").pack(side="top", anchor="w", padx=10)

        btns = ttk.Frame(self)
        btns.pack(side="bottom", fill="x", padx=10, pady=8)
        ttk.Button(btns, text="Use selected", command=self._use).pack(side="right", padx=4)
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right", padx=4)

        body = ttk.Frame(self)
        body.pack(side="top", fill="both", expand=True, padx=10, pady=6)
        self.fields_box = ttk.LabelFrame(body, text="Fields")
        self.fields_box.pack(side="left", fill="both", expand=True)
        prev = ttk.LabelFrame(body, text="Resulting groups")
        prev.pack(side="left", fill="both", expand=True, padx=(8, 0))
        self.preview = tk.Text(prev, width=30, wrap="word", state="disabled")
        self.preview.pack(fill="both", expand=True)

        self._rebuild()
        self.grab_set()

    def _rep_suffix(self):
        """The replicate-suffix regex, or '' if blank/invalid."""
        rs = self.rep.get().strip()
        if not rs:
            return ""
        try:
            re.compile(rs)
            return rs
        except re.error:
            return ""

    def _rebuild(self):
        """Recompute the field decomposition for the current delimiter / replicate
        suffix and rebuild the checkbox list. The first variable field is pre-checked."""
        for w in self.fields_box.winfo_children():
            w.destroy()
        try:
            self.prefix, self.suffix, self.parts = analyze_sample_fields(
                self.sample_cols, self.delim.get(), self._rep_suffix())
        except Exception:
            self.parts = []
            self._set_preview("Could not analyze sample names.")
            return
        _delims = "".join(_DELIMS)
        anchors = [p["text"].strip(_delims) for p in self.parts
                   if p["kind"] == "const" and p["text"].strip(_delims)]
        self.info.configure(text=(
            f"{len(self.sample_cols)} sample columns.   "
            f"prefix '{self.prefix}'   suffix '{self.suffix}'\n"
            f"Consensus regions: " + ("  |  ".join(f"'{a}'" for a in anchors) if anchors else "(none)")))

        self.var_indices, self.vars = [], []
        first = True
        for i, p in enumerate(self.parts):
            if p["kind"] != "var":
                continue
            v = tk.IntVar(value=1 if first else 0)
            first = False
            self.var_indices.append(i)
            self.vars.append(v)
            sample = ", ".join(p["uniq"][:6]) + (" ..." if len(p["uniq"]) > 6 else "")
            cb = ttk.Checkbutton(self.fields_box, variable=v,
                                 text=f"{p['label']}  ({len(p['uniq'])} values: {sample})",
                                 command=self._refresh)
            cb.pack(anchor="w", padx=6, pady=2)
        if not self.vars:
            ttk.Label(self.fields_box, text="No varying fields found.").pack(anchor="w", padx=6, pady=6)
        self._refresh()

    def _selected_indices(self):
        return [idx for idx, v in zip(self.var_indices, self.vars) if v.get()]

    def _groups(self):
        if not getattr(self, "parts", None) or not self._selected_indices():
            return []
        return build_field_groups(self.prefix, self.suffix, self.parts,
                                  self._selected_indices(), self.delim.get(), self._rep_suffix())

    def _refresh(self, *_):
        groups = self._groups()
        counts = {}
        for label, pat in groups:
            counts[label] = sum(1 for c in self.sample_cols if re.search(pat, c))
        lines = [f"{lab}  ({counts.get(lab, 0)})" for lab, _ in groups]
        self._set_preview(f"{len(groups)} groups:\n\n" + "\n".join(lines) if groups
                          else "Tick at least one field.")

    def _set_preview(self, text):
        self.preview.configure(state="normal")
        self.preview.delete("1.0", "end")
        self.preview.insert("end", text)
        self.preview.configure(state="disabled")

    def _use(self):
        groups = self._groups()
        if not groups:
            messagebox.showinfo("Auto-pick groups", "Tick at least one field first.")
            return
        names = ", ".join(lab for lab, _ in groups)
        patterns = {lab: pat for lab, pat in groups}
        self.apply_cb(names, patterns)
        self.destroy()


class ScrollableFrame(ttk.Frame):
    """A vertically scrollable frame; put widgets in `.inner`."""
    def __init__(self, parent, width=460, height=480, **kw):
        super().__init__(parent, **kw)
        # A bounded canvas height keeps tall tab content from forcing the whole
        # window taller (which would push the Log box off-screen) -- it scrolls.
        self.canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0, width=width, height=height)
        self.vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        # Route yscrollcommand through _yscroll so the bar auto-hides when it fits.
        self.canvas.configure(yscrollcommand=self._yscroll)
        self.vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.inner = ttk.Frame(self.canvas)
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        # Mouse-wheel handling is done app-wide by VolcanoGUI._global_wheel, which
        # scrolls whichever ScrollableFrame the pointer is over (a per-instance
        # bind_all would let the last-created frame hijack the wheel everywhere).

    def _yscroll(self, first, last):
        """Show the scrollbar only when the content overflows the viewport."""
        if float(first) <= 0.0 and float(last) >= 1.0:
            if self.vsb.winfo_ismapped():
                self.vsb.pack_forget()
        elif not self.vsb.winfo_ismapped():
            self.vsb.pack(side="right", fill="y", before=self.canvas)
        self.vsb.set(first, last)

    def scrollable(self):
        """True when the content is taller than the viewport (wheel can move it)."""
        first, last = self.canvas.yview()
        return not (float(first) <= 0.0 and float(last) >= 1.0)


# --------------------------------------------------------------------------- #
# Main application
# --------------------------------------------------------------------------- #
class VolcanoGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DIA-NN Proteomics  -  Volcano / Bubble Explorer" + _version_suffix())
        self.geometry("1640x900")

        self.result = None      # AnalysisResult after a run
        self.cfg = None
        self.cfg_vars = {}
        self._group_patterns = {}   # {group label: exact regex} from the auto-picker
        self._drop_samples = []     # sample columns the user excluded from the analysis
        self.pca_vars = {}
        self.vol_vars = {}
        self.bub_vars = {}
        self._sar_text = None
        self._plot_canvas = None
        self._volcano_nb = None        # the volcano "plot all" notebook (if shown)
        self._current_volcano_comp = None   # last single volcano comparison shown
        self._pending_workspace = None      # settings to re-apply after a workspace-triggered run
        self.output_dir = None   # dedicated outputs folder under the data folder
        self._log_fh = None      # open file handle: <output_dir>/analysis_log.txt
        self._scrollables = []   # ScrollableFrames, for app-wide mouse-wheel routing

        self.log_q = queue.Queue()
        self.result_q = queue.Queue()

        self._build_ui()
        self.bind_all("<MouseWheel>", self._global_wheel, add="+")
        self._setup_logging()
        self._prefill_sample_data()

        self.after(120, self._drain_log)
        self.after(150, self._poll_result)

    # ----- logging -----
    def _setup_logging(self):
        handler = _QueueLogHandler(self.log_q)
        handler.setFormatter(logging.Formatter("%(levelname)s - %(message)s"))
        root = logging.getLogger()
        root.setLevel(logging.INFO)
        root.addHandler(handler)

    def _drain_log(self):
        try:
            while True:
                msg = self.log_q.get_nowait()
                self.log_text.configure(state="normal")
                self.log_text.insert("end", msg)
                self.log_text.see("end")
                self.log_text.configure(state="disabled")
                if self._log_fh is not None:
                    try:
                        self._log_fh.write(msg)
                        self._log_fh.flush()
                    except Exception:
                        pass
        except queue.Empty:
            pass
        self.after(120, self._drain_log)

    # ----- overall layout: notebook (left) + plot (right) + log (bottom) -----
    def _build_ui(self):
        # Reserve the Log box at the bottom FIRST so taller tab content can never
        # push it off-screen; then the main area fills the space above it.
        logf = ttk.LabelFrame(self, text="Log")
        logf.pack(side="bottom", fill="x", padx=8, pady=6)
        self.log_text = tk.Text(logf, height=7, wrap="word", state="disabled")
        sb = ttk.Scrollbar(logf, command=self.log_text.yview)
        self.log_text.configure(yscrollcommand=sb.set)
        sb.pack(side="right", fill="y")
        self.log_text.pack(side="left", fill="both", expand=True)

        main = ttk.PanedWindow(self, orient="horizontal")
        main.pack(side="top", fill="both", expand=True, padx=8, pady=(8, 4))
        self._paned = main

        self.nb = ttk.Notebook(main)
        nb = self.nb
        main.add(nb, weight=0)
        cfg_tab = ttk.Frame(nb)
        pca_tab = ttk.Frame(nb)
        vol_tab = ttk.Frame(nb)
        bub_tab = ttk.Frame(nb)
        lookup_tab = ttk.Frame(nb)
        nb.add(cfg_tab, text="1. Analysis configuration")
        nb.add(pca_tab, text="2. PCA")
        nb.add(vol_tab, text="3. Volcano settings")
        nb.add(bub_tab, text="4. Bubbleplot settings")
        nb.add(lookup_tab, text="5. Raw data lookup")
        self._build_config_tab(cfg_tab)
        self._build_pca_tab(pca_tab)
        self._build_volcano_tab(vol_tab)
        self._build_bubble_tab(bub_tab)
        self._build_lookup_tab(lookup_tab)

        # Right side: a persistent per-tab output area (cache). Each plotting tab
        # keeps its own last-rendered plot/table; switching tabs swaps which one
        # is shown (the notebook tab index maps to a key below).
        right = ttk.LabelFrame(main, text="Plot / table")
        main.add(right, weight=1)
        self.plot_host = ttk.Frame(right)
        self.plot_host.pack(side="top", fill="both", expand=True)
        self.tab_frames = {k: ttk.Frame(self.plot_host) for k in ("pca", "volcano", "bubble", "lookup")}
        for k, f in self.tab_frames.items():
            ttk.Label(f, text="Nothing rendered yet — use the button on this tab after a run.",
                      foreground="#888").pack(expand=True)
        self.placeholder_frame = ttk.Frame(self.plot_host)
        ttk.Label(self.placeholder_frame,
                  text="Run an analysis (tab 1), then plot from tabs 2-4 or look up raw data on tab 5.",
                  foreground="#888").pack(expand=True)
        self._current_area = None
        self._tab_index_key = {1: "pca", 2: "volcano", 3: "bubble", 4: "lookup"}
        self._show_tab_area(None)
        nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)

    def _show_tab_area(self, key):
        """Show the output frame for `key` (None -> the generic placeholder)."""
        target = self.tab_frames.get(key, self.placeholder_frame)
        if target is self._current_area:
            return
        if self._current_area is not None:
            self._current_area.pack_forget()
        target.pack(fill="both", expand=True)
        self._current_area = target

    def _tab_area(self, key):
        """Return the (cleared) output frame for a plotting tab."""
        frame = self.tab_frames[key]
        for w in frame.winfo_children():
            w.destroy()
        return frame

    def _on_tab_changed(self, _event=None):
        idx = self.nb.index(self.nb.select())
        key = self._tab_index_key.get(idx)
        self._show_tab_area(key)
        # Entering the lookup tab: default its comparison to the one currently shown
        # on the volcano panel (the active "plot all" sub-tab, or the last single plot).
        if key == "lookup":
            comp = self._current_volcano_comparison()
            if comp:
                try:
                    if comp in list(self.lookup_comp.cget("values")):
                        self.lookup_comp.set(comp)
                except Exception:
                    pass

    def _current_volcano_comparison(self):
        """The comparison currently displayed on the volcano panel, or None."""
        nb = self._volcano_nb
        try:
            if nb is not None and nb.winfo_exists():
                sel = nb.select()
                if sel:
                    return nb.tab(sel, "text")
        except Exception:
            pass
        return self._current_volcano_comp

    def _scroll_inner(self, tab):
        sf = ScrollableFrame(tab, width=470)
        sf.pack(side="top", fill="both", expand=True)
        self._scrollables.append(sf)
        return sf.inner

    def _global_wheel(self, event):
        """Scroll whichever ScrollableFrame the pointer is currently over -- unless
        a Listbox is under the pointer (an open combobox dropdown, or the drop-list),
        which keeps the wheel for itself instead of leaking it to the panel behind."""
        try:
            path = self.tk.call("winfo", "containing", event.x_root, event.y_root)
            if path and self.tk.call("winfo", "class", path) == "Listbox":
                self.tk.call(path, "yview", "scroll", int(-event.delta / 120), "units")
                return "break"
        except tk.TclError:
            pass
        for sf in self._scrollables:
            c = sf.canvas
            if not c.winfo_ismapped():
                continue
            rx, ry = c.winfo_rootx(), c.winfo_rooty()
            if rx <= event.x_root < rx + c.winfo_width() and ry <= event.y_root < ry + c.winfo_height():
                if sf.scrollable():               # don't move panels that already fit
                    c.yview_scroll(int(-event.delta / 120), "units")
                return "break"

    # ----- tab 1: analysis configuration -----
    def _build_config_tab(self, tab):
        bar = ttk.Frame(tab)
        bar.pack(side="bottom", fill="x", padx=4, pady=6)
        self.preview_btn = ttk.Button(bar, text="Preview groups", command=self._preview_groups)
        self.preview_btn.pack(side="left", padx=4)
        Tooltip(self.preview_btn, "Show how many samples match each group, computed from the file's "
                                  "column headers, without running the full analysis.")
        self.run_btn = ttk.Button(bar, text="Run analysis", command=self._on_run)
        self.run_btn.pack(side="left", padx=4)
        Tooltip(self.run_btn, "Load the data and run the pipeline (impute -> fold change -> t-test/limma) "
                              "in the background. Outputs and a log go to proteomics_GUI_output.")
        self.status = ttk.Label(bar, text="Not run yet", foreground="#a60")
        self.status.pack(side="left", padx=8)

        # Workspace recall. A small workspace.json (all settings + pointers to the
        # data files, not copies) is auto-saved to the output folder after each run,
        # so there's no Save button -- just reopen one here.
        self.load_ws_btn = ttk.Button(bar, text="Load workspace", command=self._on_load_workspace)
        self.load_ws_btn.pack(side="right", padx=4)
        Tooltip(self.load_ws_btn, "Reopen a saved workspace.json: restores every setting and rebuilds the "
                                  "plots/lookup from the saved results in the output folder -- no re-run "
                                  "(it only re-runs if those outputs are missing). "
                                  "A workspace is auto-saved after each run.")

        gframe = ttk.LabelFrame(tab, text="Group assignments (samples matched per group)")
        gframe.pack(side="bottom", fill="both", padx=4, pady=4)
        self._group_text = tk.Text(gframe, height=9, wrap="none", state="disabled")
        gsb = ttk.Scrollbar(gframe, command=self._group_text.yview)
        self._group_text.configure(yscrollcommand=gsb.set)
        gsb.pack(side="right", fill="y")
        self._group_text.pack(side="left", fill="both", expand=True)

        inner = self._scroll_inner(tab)
        grid = ttk.Frame(inner)
        grid.pack(fill="x", padx=4, pady=4)

        v = self.cfg_vars
        v["pg_path"] = path_row(grid, 0, "Protein file (.pg_matrix.tsv)", self._browse_pg,
                                hint="Your DIA-NN protein-group matrix (<stem>.pg_matrix.tsv). "
                                     "The folder it lives in becomes the working folder, and outputs "
                                     "go to a proteomics_GUI_output folder there.")
        v["pr_path"] = path_row(grid, 1, "Precursor file (.pr_matrix.tsv)", self._browse_pr,
                                hint="The matching <stem>.pr_matrix.tsv (precursor/peptide matrix). "
                                     "Auto-filled from the protein file. Used for the peptide-count "
                                     "check that decides imputation of high-missing proteins.")
        self.path_label = ttk.Label(grid, text="No data selected.", foreground="#555", wraplength=440, justify="left")
        self.path_label.grid(row=2, column=0, columnspan=3, sticky="w", padx=4, pady=(0, 6))
        lbl_groups = ttk.Label(grid, text="Groups (comma\nor one per line)")
        lbl_groups.grid(row=3, column=0, sticky="nw", padx=4, pady=2)
        v["group_names"] = tk.StringVar(value="DMSO, Positive_Control, IRAK1")
        # Multi-row, scrolling, drag-to-resize box (group lists can be long, e.g. a
        # 30-compound screen). It stays in sync with the group_names StringVar that
        # the rest of the app reads/writes, so nothing else has to change.
        gbox = ttk.Frame(grid)
        gbox.grid(row=3, column=1, sticky="we", padx=4, pady=2)
        gbox.columnconfigure(0, weight=1)
        self.groups_text = tk.Text(gbox, height=3, width=34, wrap="word", undo=True)
        gvsb = ttk.Scrollbar(gbox, orient="vertical", command=self.groups_text.yview)
        self.groups_text.configure(yscrollcommand=gvsb.set)
        self.groups_text.grid(row=0, column=0, sticky="we")
        gvsb.grid(row=0, column=1, sticky="ns")
        self.groups_text.insert("1.0", v["group_names"].get())
        grip = ttk.Separator(gbox, orient="horizontal")
        grip.grid(row=1, column=0, columnspan=2, sticky="we", pady=(2, 0))
        grip.configure(cursor="sb_v_double_arrow")

        def _grip_start(e):
            self._grip_y0 = e.y_root
            self._grip_h0 = int(float(self.groups_text.cget("height")))

        def _grip_drag(e):
            self.groups_text.configure(height=max(1, self._grip_h0 + (e.y_root - self._grip_y0) // 16))
        grip.bind("<Button-1>", _grip_start)
        grip.bind("<B1-Motion>", _grip_drag)

        # Two-way sync. Typing also drops any exact patterns the auto-picker set
        # (the picker writes via .set(), which fires the trace but not KeyRelease,
        # so its patterns survive a programmatic update).
        self._syncing_groups = False

        def _groups_text_to_var(*_):
            if self._syncing_groups:
                return
            self._syncing_groups = True
            v["group_names"].set(self.groups_text.get("1.0", "end").strip())
            self._group_patterns = {}
            self._syncing_groups = False

        def _groups_var_to_text(*_):
            if self._syncing_groups:
                return
            new = v["group_names"].get()
            if self.groups_text.get("1.0", "end").strip() != new:
                self._syncing_groups = True
                self.groups_text.delete("1.0", "end")
                self.groups_text.insert("1.0", new)
                self._syncing_groups = False
        self.groups_text.bind("<KeyRelease>", _groups_text_to_var)
        v["group_names"].trace_add("write", _groups_var_to_text)

        btn_auto = ttk.Button(grid, text="Auto-pick...", command=self._autopick_groups)
        btn_auto.grid(row=3, column=2, sticky="nw", padx=4)
        _gh = ("Group names, separated by commas or newlines. Each is matched as a regex against the "
               "sample-column headers (or exactly, if set via Auto-pick). Drag the bar below the box to "
               "resize it; it scrolls when the list is long.")
        Tooltip(lbl_groups, _gh)
        Tooltip(btn_auto, "Split the sample-column headers into fields and tick which field(s) define a "
                          "group (the rest are treated as replicates). Groups get clean labels but are "
                          "matched exactly, so 'NR' won't leak into 'NR-TDxdR' columns.")
        v["reference_group"] = labeled_entry(grid, 4, "Reference group", "DMSO", width=18,
                                             hint="The control group every treatment is compared against "
                                                  "(e.g. DMSO). Must be one of the group names above.")
        v["comparison_matrix"] = labeled_entry(grid, 5, "Comparisons", "", width=34,
                                               tip="blank = all vs reference",
                                               hint="treated:control pairs, comma-separated "
                                                    "(e.g. IRAK1:DMSO, Positive_Control:DMSO). "
                                                    "Leave blank to compare every group against the reference.")
        ttk.Label(grid, text="(format: treated:control, comma-separated)", foreground="#666").grid(
            row=6, column=1, columnspan=2, sticky="w", padx=4)
        v["mode"] = labeled_entry(grid, 7, "Mode (0/1)", "0", width=6, tip="0 = degradation, 1 = enrichment",
                                  hint="0 = degradation (impute the treated group); "
                                       "1 = enrichment/pulldown (impute the control group).")
        v["control_threshold"] = labeled_entry(grid, 8, "Ctrl detect thresh", "0.5", width=6,
                                               hint="Minimum fraction of reference samples a protein must be "
                                                    "detected in to be kept. Only applied when Comparisons is blank.")
        v["normalization_protein_id"] = labeled_entry(grid, 9, "Normalize to UniProt", "", width=16, tip="blank = none",
                                                      hint="UniProt ID to normalize every sample to (e.g. a loading "
                                                           "control / spike-in). Blank = no normalization.")

        # Drop samples: pick a sample column from the (live) header dropdown and Add
        # it to the exclusion list. Dropped columns are removed from both matrices
        # right after load, so grouping/PCA/imputation never see them.
        dropf = ttk.LabelFrame(inner, text="Drop samples from analysis (optional)")
        dropf.pack(fill="x", padx=4, pady=4)
        drow = ttk.Frame(dropf)
        drow.pack(fill="x", padx=4, pady=(4, 2))
        ttk.Label(drow, text="Sample:").pack(side="left")
        self.drop_combo = ttk.Combobox(drow, state="readonly", width=46, values=[],
                                       postcommand=self._refresh_sample_list)
        self.drop_combo.pack(side="left", padx=4)
        add_btn = ttk.Button(drow, text="Add", width=6, command=self._drop_add)
        add_btn.pack(side="left", padx=2)
        Tooltip(self.drop_combo, "Sample columns read live from the selected pg_matrix header. "
                                 "Pick one and click Add to exclude it from the analysis.")
        lrow = ttk.Frame(dropf)
        lrow.pack(fill="both", expand=True, padx=4, pady=(0, 2))
        self.drop_list = tk.Listbox(lrow, height=3, selectmode="extended", exportselection=False)
        dsb = ttk.Scrollbar(lrow, command=self.drop_list.yview)
        self.drop_list.configure(yscrollcommand=dsb.set)
        self.drop_list.pack(side="left", fill="both", expand=True)
        dsb.pack(side="left", fill="y")
        brow = ttk.Frame(dropf)
        brow.pack(fill="x", padx=4, pady=(0, 4))
        ttk.Button(brow, text="Remove selected", command=self._drop_remove).pack(side="left", padx=2)
        ttk.Button(brow, text="Clear all", command=self._drop_clear).pack(side="left", padx=2)
        self._drop_hint = ttk.Label(brow, text="0 dropped", foreground="#666")
        self._drop_hint.pack(side="left", padx=8)

        checks = ttk.Frame(inner)
        checks.pack(fill="x", padx=4, pady=6)
        v["imputation_option"] = check(checks, 0, "Imputation", True,
                                       hint="Impute missing values (treated and reference groups). "
                                            "Off = use raw values only; proteins with too few valid values are dropped.")
        v["output_adjpval"] = check(checks, 1, "Plot adjusted P (FDR)", True,
                                    hint="Use the adjusted P-value (BH/FDR) on the volcano y-axis. "
                                         "Off = use the raw limma p-value. Both are written to the CSV either way.")
        ttk.Label(checks, text="Statistics: R/limma moderated t-test (requires R + R_HOME + limma).",
                  foreground="#666").grid(row=2, column=0, columnspan=3, sticky="w", padx=4, pady=(4, 0))
        # Pharos TCRD is purely a volcano-plot coloring option -- it lives on the
        # Volcano tab (Highlights), not here, since it doesn't affect the analysis.

    def _set_group_text(self, group_columns):
        """Fill the dedicated group-assignments box with per-group sample counts."""
        self._group_text.configure(state="normal")
        self._group_text.delete("1.0", "end")
        if not group_columns:
            self._group_text.insert("end", "No groups matched any sample columns.\n")
        else:
            total = sum(len(c) for c in group_columns.values())
            self._group_text.insert("end", f"{len(group_columns)} groups, {total} samples matched:\n\n")
            for grp, cols in group_columns.items():
                self._group_text.insert("end", f"{grp}: {len(cols)}\n")
                for c in cols:
                    self._group_text.insert("end", f"    {c}\n")
        self._group_text.configure(state="disabled")

    _META_COLS = {"Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"}

    def _sample_columns(self):
        """Non-meta sample columns from the current pg file header ([] on failure)."""
        pg = self.cfg_vars["pg_path"].get().strip()
        if not pg or not os.path.exists(pg):
            return []
        try:
            import pandas as pd
            cols = list(pd.read_csv(pg, sep="\t", index_col=0, nrows=0).columns)
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            return []
        return [c for c in cols if c not in self._META_COLS]

    # ----- drop samples from the analysis -----
    def _refresh_sample_list(self):
        """Repopulate the drop-sample dropdown from the pg header, minus what's
        already dropped. Runs as the combobox's postcommand (on each open)."""
        avail = [c for c in self._sample_columns() if c not in self._drop_samples]
        self.drop_combo.configure(values=avail)

    def _update_drop_hint(self):
        self._drop_hint.configure(text=f"{len(self._drop_samples)} dropped")

    def _drop_add(self):
        s = self.drop_combo.get().strip()
        if s and s not in self._drop_samples:
            self._drop_samples.append(s)
            self.drop_list.insert("end", s)
            self.drop_combo.set("")
            self._update_drop_hint()

    def _drop_remove(self):
        for i in reversed(self.drop_list.curselection()):
            s = self.drop_list.get(i)
            self.drop_list.delete(i)
            if s in self._drop_samples:
                self._drop_samples.remove(s)
        self._update_drop_hint()

    def _drop_clear(self):
        self._drop_samples = []
        self.drop_list.delete(0, "end")
        self._update_drop_hint()

    def _set_drop_samples(self, samples):
        """Replace the dropped-samples list (used by workspace load)."""
        self._drop_samples = [str(s) for s in (samples or [])]
        self.drop_list.delete(0, "end")
        for s in self._drop_samples:
            self.drop_list.insert("end", s)
        self._update_drop_hint()

    def _autopick_groups(self):
        """Open a popup to pick group names from the pg matrix header columns."""
        pg = self.cfg_vars["pg_path"].get().strip()
        if not pg or not os.path.exists(pg):
            messagebox.showerror("Auto-pick groups", "Select a valid pg_matrix.tsv file first.")
            return
        samples = self._sample_columns()
        if len(samples) < 2:
            messagebox.showinfo("Auto-pick groups", "Not enough sample columns to analyze.")
            return
        def _apply(names, patterns):
            self.cfg_vars["group_names"].set(names)
            self._group_patterns = patterns
        GroupPickerDialog(self, samples, _apply)

    def _preview_groups(self):
        """Compute group assignments from the pg file's columns, without running."""
        pg = self.cfg_vars["pg_path"].get().strip()
        if not pg or not os.path.exists(pg):
            messagebox.showerror("Preview groups", "Select a valid pg_matrix.tsv file first.")
            return
        groups = [g.strip() for g in re.split(r"[,\n]", self.cfg_vars["group_names"].get()) if g.strip()]
        if not groups:
            messagebox.showerror("Preview groups", "Enter at least one group name.")
            return
        try:
            import pandas as pd
            from scripts.io import assign_groups
            df0 = pd.read_csv(pg, sep="\t", index_col=0, nrows=0)  # header only
            if self._drop_samples:                                # mirror the analysis-time drop
                df0 = df0.drop(columns=[c for c in self._drop_samples if c in df0.columns])
            patterns = {g: self._group_patterns[g] for g in groups if g in self._group_patterns} or None
            self._set_group_text(assign_groups(df0, groups, patterns))
            logging.getLogger().info("Previewed group assignments from %s", os.path.basename(pg))
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            messagebox.showerror("Preview groups", "Could not read columns. See the Log.")

    # ----- tab 2: PCA -----
    def _build_pca_tab(self, tab):
        bar = ttk.Frame(tab)
        bar.pack(side="bottom", fill="x", padx=4, pady=6)
        self.pca_btn = ttk.Button(bar, text="Plot PCA", command=self._on_plot_pca, state="disabled")
        self.pca_btn.pack(side="left", padx=4)
        Tooltip(self.pca_btn, "Render the PCA plot with the options below (enabled after a run).")
        ttk.Label(bar, text="  (hover a dot for its sample; click a dot to toggle it in the "
                            "drop-off list in step 1)", foreground="#888").pack(side="left", padx=4)

        inner = self._scroll_inner(tab)
        box = ttk.LabelFrame(inner, text="PCA options")
        box.pack(fill="x", padx=4, pady=4)
        p = self.pca_vars
        p["title"] = labeled_entry(box, 0, "Title", "PCA plot", width=26,
                                   hint="Plot title shown at the top of the PCA figure.")
        p["filename"] = labeled_entry(box, 1, "Output filename", "PCA_plot.png", width=26,
                                      hint="File name for the saved PCA image, written to the outputs folder.")
        p["text"] = check(box, 2, "Label samples on the plot", True,
                          hint="Annotate each point with its sample name (uses adjustText to avoid overlap).")
        p["dpi"] = labeled_combo(box, 3, "DPI", [100, 150, 200, 300, 600, 1200], 300,
                                 hint="Resolution of the saved PNG (dots per inch). Higher = sharper, bigger file.")

        fonts = ttk.LabelFrame(inner, text="Font sizes")
        fonts.pack(fill="x", padx=4, pady=4)
        p["title_fontsize"] = labeled_entry(fonts, 0, "Title", "20", hint="Font size of the plot title.")
        p["axis_fontsize"] = labeled_entry(fonts, 1, "Axis labels", "15",
                                           hint="Font size of the PC1 / PC2 axis labels.")
        p["tick_fontsize"] = labeled_entry(fonts, 2, "Tick labels", "", tip="blank = default",
                                           hint="Font size of the axis tick numbers. Blank = matplotlib default.")
        p["legend_fontsize"] = labeled_entry(fonts, 3, "Legend", "", tip="blank = default",
                                             hint="Font size of the group legend. Blank = matplotlib default.")
        p["point_label_fontsize"] = labeled_entry(fonts, 4, "Sample labels", "4",
                                                  hint="Font size of the per-sample text labels (when enabled above).")
        ttk.Label(inner,
                  text="PCA uses every sample across all groups; missing values are\n"
                       "mean-imputed per protein (as in the original analysis). It does\n"
                       "not depend on the volcano/bubble settings.",
                  foreground="#666", justify="left").pack(anchor="w", padx=8, pady=6)

    # ----- tab 3: volcano settings -----
    def _build_volcano_tab(self, tab):
        bar = ttk.Frame(tab)
        bar.pack(side="bottom", fill="x", padx=4, pady=6)
        self.plot_all_btn = ttk.Button(bar, text="Plot all comparisons", command=self._on_plot_all, state="disabled")
        self.plot_all_btn.pack(side="left", padx=3)
        Tooltip(self.plot_all_btn, "Plot a volcano for every comparison (one tab each) using the current settings.")
        self.plot_btn = ttk.Button(bar, text="Plot selected", command=self._on_plot_volcano, state="disabled")
        self.plot_btn.pack(side="left", padx=3)
        Tooltip(self.plot_btn, "Plot a volcano for just the treatment/control pair selected below.")
        ttk.Label(bar, text="  (hover a dot for its gene; click a dot to look it up in step 5)",
                  foreground="#888").pack(side="left", padx=4)

        sel = ttk.Frame(tab)
        sel.pack(side="top", fill="x", padx=4, pady=4)
        ttk.Label(sel, text="Single comparison (for 'Plot selected'):", foreground="#666").grid(
            row=0, column=0, columnspan=2, sticky="w", padx=4)
        ttk.Label(sel, text="Treatment").grid(row=1, column=0, sticky="w", padx=4)
        self.treat_cb = ttk.Combobox(sel, state="readonly", width=24, values=[])
        self.treat_cb.grid(row=1, column=1, padx=4, pady=2)
        ttk.Label(sel, text="Control").grid(row=2, column=0, sticky="w", padx=4)
        self.ctrl_cb = ttk.Combobox(sel, state="readonly", width=24, values=[])
        self.ctrl_cb.grid(row=2, column=1, padx=4, pady=2)

        self._build_volcano_widgets(self._scroll_inner(tab))

    @staticmethod
    def _bind_show(var, frame):
        """Show `frame` (grid) only while the boolean `var` is True."""
        def _upd(*_a):
            frame.grid() if var.get() else frame.grid_remove()
        var.trace_add("write", _upd)
        _upd()

    @staticmethod
    def _bind_enable(var, frame):
        """Grey out every widget inside `frame` while the boolean `var` is False."""
        def _apply(widget, on):
            for child in widget.winfo_children():
                try:                       # ttk widgets
                    child.state(["!disabled" if on else "disabled"])
                except Exception:
                    try:                   # classic tk widgets
                        child.configure(state="normal" if on else "disabled")
                    except Exception:
                        pass
                _apply(child, on)
        def _upd(*_a):
            _apply(frame, bool(var.get()))
        var.trace_add("write", _upd)
        _upd()

    def _build_volcano_widgets(self, parent):
        v = self.vol_vars

        thr = ttk.LabelFrame(parent, text="Thresholds")
        thr.pack(fill="x", padx=4, pady=4)
        v["logFC_cutoff"] = labeled_entry(thr, 0, "log2FC cutoff", "1", tip="blank = none",
                                          hint="Vertical cutoff: |log2FC| >= this (with FDR cutoff) marks up/down hits "
                                               "and draws the dashed lines. Blank = no fold-change cutoff.")
        v["logFC_cutoff2"] = labeled_entry(thr, 1, "log2FC cutoff 2", "", tip="secondary, blank = none",
                                           hint="Optional second, looser down-regulation cutoff (e.g. >35% down) shown "
                                                "as small blue points. Blank = off.")
        v["FDR_cutoff"] = labeled_entry(thr, 2, "FDR cutoff", "0.05",
                                        hint="Significance threshold on the adjusted P-value / FDR; also the horizontal "
                                             "dashed line.")
        v["protein_level_cutoff"] = labeled_entry(thr, 3, "Protein-level cutoff", "", tip="blank = off",
                                                  hint="Highlight down-regulated proteins whose mean control abundance "
                                                       "is below 1000 in a separate colour. Blank = off.")

        emp = ttk.LabelFrame(parent, text="Empirical FDR")
        emp.pack(fill="x", padx=4, pady=4)
        v["use_empirical_fdr"] = check(emp, 0, "Use empirical FDR curve", False,
                                       hint="Replace the straight cutoffs with a hyperbolic curve fitted by the "
                                            "wrong-side (decoy) method below, instead of fixed log2FC/FDR lines.")
        empf = ttk.Frame(emp)   # these only apply when the curve is enabled -> grey out otherwise
        empf.grid(row=1, column=0, columnspan=3, sticky="w")
        v["fdr_alpha"] = labeled_entry(empf, 0, "FDR alpha", "0.05",
                                       hint="Target empirical FDR for the fitted curve (e.g. 0.05).")
        v["kappa"] = labeled_entry(empf, 1, "kappa", "1e-6",
                                   hint="Small pseudocount used when estimating the empirical FDR ratio.")
        v["p_value_cutoff"] = labeled_entry(empf, 2, "p_value_cutoff", "1",
                                            hint="Vertical offset (in -log10 P) of the empirical-FDR curve's asymptote.")
        v["mode"] = labeled_entry(empf, 3, "Mode override", "Auto", tip="Auto / 0 / 1",
                                  hint="Side the curve is fit toward. Auto = use the analysis mode "
                                       "(0 degradation -> down side, 1 enrichment -> up side).")
        self._bind_enable(v["use_empirical_fdr"], empf)

        axes = ttk.LabelFrame(parent, text="Axes")
        axes.pack(fill="x", padx=4, pady=4)
        v["xlim_min"] = labeled_entry(axes, 0, "x min", "-8", hint="Left limit of the log2FC (x) axis.")
        v["xlim_max"] = labeled_entry(axes, 1, "x max", "8", hint="Right limit of the log2FC (x) axis.")
        v["ylim_min"] = labeled_entry(axes, 2, "y min", "", tip="blank = auto",
                                      hint="Bottom limit of the -log10(FDR) (y) axis. Blank = auto.")
        v["ylim_max"] = labeled_entry(axes, 3, "y max", "", tip="blank = auto",
                                      hint="Top limit of the y axis. Blank = auto (max + top buffer).")
        v["x_interval"] = labeled_entry(axes, 4, "x tick interval", "2", hint="Spacing between x-axis ticks.")
        v["y_interval"] = labeled_entry(axes, 5, "y tick interval", "1",
                                        hint="Spacing between y-axis ticks (used when y limits are set).")
        v["top_buffer"] = labeled_entry(axes, 6, "Top buffer", "0.1",
                                        hint="Extra headroom above the tallest point, as a fraction (auto y only).")

        hi = ttk.LabelFrame(parent, text="Highlights")
        hi.pack(fill="x", padx=4, pady=4)
        v["imputation_option"] = check(hi, 0, "Mark imputed proteins (orange)", True,
                                       hint="Colour proteins that were imputed (special high-missing protocol); "
                                            "the colour is the 'Imputed colour' in the Dots section.")

        # Each of these highlight sets reveals its colour picker(s) when ticked.
        v["PharosTCRD"] = check(hi, 1, "Pharos TCRD classes", False,
                                hint="Colour proteins by Pharos target-development level (Tclin/Tchem/Tbio/Tdark).")
        ph = ttk.Frame(hi)
        ph.grid(row=2, column=0, columnspan=3, sticky="w", padx=18)
        v["color_tclin"] = color_row(ph, 0, "Tclin colour", "#17becf")
        v["color_tchem"] = color_row(ph, 1, "Tchem colour", "#e377c2")
        v["color_tbio"] = color_row(ph, 2, "Tbio colour", "#8c564b")
        v["color_tdark"] = color_row(ph, 3, "Tdark colour", "#7f7f7f")
        self._bind_show(v["PharosTCRD"], ph)

        v["highlight_kinase"] = check(hi, 3, "Protein kinases", False,
                                      hint="Mark proteins in the protein-kinase reference list.")
        kf = ttk.Frame(hi)
        kf.grid(row=4, column=0, columnspan=3, sticky="w", padx=18)
        v["color_kinase"] = color_row(kf, 0, "Kinase colour", "#17becf")
        self._bind_show(v["highlight_kinase"], kf)

        v["highlight_ub"] = check(hi, 5, "Ubiquitin-related", False,
                                  hint="Mark ubiquitin-related proteins (and add the significant ones to the labels).")
        uf = ttk.Frame(hi)
        uf.grid(row=6, column=0, columnspan=3, sticky="w", padx=18)
        v["color_ub"] = color_row(uf, 0, "Ubiquitin colour", "#17becf")
        self._bind_show(v["highlight_ub"], uf)

        v["highlight_Gloops"] = check(hi, 7, "G-loop proteins", False,
                                      hint="Mark proteins in the G-loop reference list.")
        gf = ttk.Frame(hi)
        gf.grid(row=8, column=0, columnspan=3, sticky="w", padx=18)
        v["color_gloops"] = color_row(gf, 0, "G-loop colour", "#17becf")
        self._bind_show(v["highlight_Gloops"], gf)

        v["highlight_RTloops"] = check(hi, 9, "RT-loop proteins", False,
                                       hint="Mark proteins in the RT-loop reference list.")
        rf = ttk.Frame(hi)
        rf.grid(row=10, column=0, columnspan=3, sticky="w", padx=18)
        v["color_rtloops"] = color_row(rf, 0, "RT-loop colour", "#e377c2")
        self._bind_show(v["highlight_RTloops"], rf)

        v["highlight_genes"] = labeled_entry(hi, 11, "Highlight genes", "", width=28,
                                             tip="UniProt IDs or gene names, comma-separated",
                                             hint="Specific proteins to mark and always label (colour = 'Highlight "
                                                  "colour' in the Dots section). Comma-separated UniProt accessions "
                                                  "(e.g. P51617) or gene names (e.g. IRAK1); gene names are "
                                                  "case-insensitive.")

        dots = ttk.LabelFrame(parent, text="Dots (size, transparency, colours)")
        dots.pack(fill="x", padx=4, pady=4)
        v["dot_size"] = labeled_combo(dots, 0, "Dot size", [10, 20, 30, 40, 50, 60, 80, 100], 40,
                                      hint="Marker size of the main dots (background, up, down).")
        v["dot_alpha"] = labeled_combo(dots, 1, "Transparency (alpha)",
                                       ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"], "0.5",
                                       hint="Opacity of the main dots: 1.0 = solid, lower = more see-through.")
        v["color_up"] = color_row(dots, 2, "Up-regulated colour", "red",
                                  hint="Colour of significant up-regulated dots (click to pick).")
        v["color_down"] = color_row(dots, 3, "Down-regulated colour", "blue",
                                    hint="Colour of significant down-regulated dots (click to pick).")
        v["color_imputed"] = color_row(dots, 4, "Imputed colour", "orange",
                                       hint="Colour of the imputed-protein markers (click to pick).")
        v["color_highlight"] = color_row(dots, 5, "Highlight colour", "green",
                                         hint="Colour of the 'Highlight genes' markers (click to pick).")
        v["color_bg"] = color_row(dots, 6, "Background colour", "grey",
                                  hint="Colour of the non-significant background dots (click to pick).")

        gl = ttk.LabelFrame(parent, text="Gene labels (which to show)")
        gl.pack(fill="x", padx=4, pady=4)
        v["label_up"] = check(gl, 0, "Label up-regulated genes", True,
                              hint="Add the significant up-regulated (red) genes to the text labels.")
        v["label_down"] = check(gl, 1, "Label down-regulated genes", True,
                                hint="Add the significant down-regulated (blue) genes to the text labels.")
        v["label_imputed"] = check(gl, 2, "Label imputed genes", False,
                                   hint="Add the imputed (orange) proteins to the text labels.")
        v["label_topX_mid_fc"] = labeled_entry(gl, 3, "Label top-X mid FC", "", tip="blank = off",
                                               hint="Also label the X most significant 'mid' down-regulated proteins "
                                                    "(log2FC between -1 and -0.32). Blank/0 = off.")
        v["label_most_extreme"] = labeled_entry(gl, 4, "Label most extreme (per side)", "", tip="blank = off",
                                                hint="Label only the N points farthest from the origin on each side "
                                                     "(overrides the up/down label set). Blank = off.")
        v["max_label"] = labeled_entry(gl, 5, "Max labels", "100",
                                       hint="Skip drawing labels entirely if more than this many would be shown "
                                            "(prevents an unreadable, slow plot).")

        place = ttk.LabelFrame(parent, text="Label placement")
        place.pack(fill="x", padx=4, pady=4)
        v["adjust_labels"] = check(place, 0, "Auto-arrange labels (adjustText)", True,
                                   hint="Reposition labels to avoid overlap (adjustText). Off = place at the point, "
                                        "no arrows -- much faster with many labels.")
        placf = ttk.Frame(place)   # arrow/force options only apply when auto-arrange is on
        placf.grid(row=1, column=0, columnspan=3, sticky="w")
        v["adjust_arrows"] = check(placf, 0, "Draw arrows to labels", True,
                                   hint="Draw thin connector lines from each moved label back to its point.")
        v["adjust_force_text"] = labeled_entry(placf, 1, "Repel force (text)", "1, 2",
                                               tip="number or 'x, y'",
                                               hint="How strongly labels push apart from each other. "
                                                    "A number, or 'x, y' for separate horizontal/vertical force.")
        v["adjust_force_static"] = labeled_entry(placf, 2, "Repel force (points)", "1, 2",
                                                 tip="number or 'x, y'",
                                                 hint="How strongly labels are pushed away from the data points.")
        self._bind_enable(v["adjust_labels"], placf)

        fonts = ttk.LabelFrame(parent, text="Font sizes")
        fonts.pack(fill="x", padx=4, pady=4)
        v["title_fontsize"] = labeled_entry(fonts, 0, "Title", "24", hint="Font size of the plot title.")
        v["axis_label_fontsize"] = labeled_entry(fonts, 1, "Axis labels", "20",
                                                 hint="Font size of the x/y axis labels.")
        v["tick_fontsize"] = labeled_entry(fonts, 2, "Tick labels", "16", hint="Font size of the axis tick numbers.")
        v["legend_fontsize"] = labeled_entry(fonts, 3, "Legend", "12", hint="Font size of the legend.")
        v["gene_label_fontsize"] = labeled_entry(fonts, 4, "Gene labels", "14",
                                                 hint="Font size of the per-gene text labels on the plot.")

        out = ttk.LabelFrame(parent, text="Output (saved PNG)")
        out.pack(fill="x", padx=4, pady=4)
        v["dpi"] = labeled_combo(out, 0, "DPI", [100, 150, 200, 300, 600, 1200], 300,
                                 hint="Resolution of the saved PNG (dots per inch). Higher = sharper, bigger file.")
        v["file_suffix"] = labeled_entry(out, 1, "File suffix", "",
                                         hint="Extra text appended to the saved PNG file name.")

    # ----- tab 4: bubbleplot settings -----
    def _build_bubble_tab(self, tab):
        bar = ttk.Frame(tab)
        bar.pack(side="bottom", fill="x", padx=4, pady=6)
        self.bubble_btn = ttk.Button(bar, text="Plot bubble", command=self._on_plot_bubble, state="disabled")
        self.bubble_btn.pack(side="left", padx=4)
        Tooltip(self.bubble_btn, "Render the clustered bubble/dendrogram plot from the SAR groups and "
                                 "options below (enabled after a run).")

        self._build_bubble_widgets(self._scroll_inner(tab))

    def _build_bubble_widgets(self, parent):
        v = self.bub_vars

        sar = ttk.LabelFrame(parent, text="SAR groups (downregulated proteins clustered across these)")
        sar.pack(fill="x", padx=4, pady=4)
        ttk.Label(sar, text="One group per line:   label: treatmentA, treatmentB",
                  foreground="#666").pack(anchor="w", padx=4)
        ttk.Label(sar, text="(label may be empty; treatments are the column names, suffix added below)",
                  foreground="#888").pack(anchor="w", padx=4)
        self._sar_text = tk.Text(sar, height=5, width=46, wrap="word")
        self._sar_text.pack(fill="x", padx=4, pady=3)
        suf = ttk.Frame(sar)
        suf.pack(fill="x")
        v["sar_suffix"] = labeled_entry(suf, 0, "Suffix appended to each", "_vs_DMSO", width=18,
                                        hint="Appended to every treatment name above to form the comparison suffix, "
                                             "e.g. IRAK1 + '_vs_DMSO' -> log2FC_IRAK1_vs_DMSO.")

        fig = ttk.LabelFrame(parent, text="Figure")
        fig.pack(fill="x", padx=4, pady=4)
        v["fig_title"] = labeled_entry(fig, 0, "Title", "", width=28, hint="Title shown above the plot.")
        v["figure_filename"] = labeled_entry(fig, 1, "Output filename", "bubble_plot.png", width=24,
                                             hint="File name for the saved bubble image (in the outputs folder).")
        v["fig_width"] = labeled_entry(fig, 2, "Width (in)", "18",
                                       hint="Figure width in inches. Increase for many compounds/proteins.")
        v["fig_height"] = labeled_entry(fig, 3, "Height (in)", "12", hint="Figure height in inches.")
        v["dendro_bubble_height_ratio"] = labeled_entry(fig, 4, "Dendro:bubble height", "1, 3",
                                                        hint="Height ratio of the dendrogram to the bubble grid, "
                                                             "e.g. '1, 3'. Ignored when axes are inverted.")
        v["bubble_legend_width_ratio"] = labeled_entry(fig, 5, "Bubble:legend width", "20, 1",
                                                       hint="Width ratio of the bubble grid to the legend/colorbar column.")
        v["colorFCrange"] = labeled_entry(fig, 6, "Color FC range", "-4, 0",
                                          hint="log2FC range mapped to the colour scale, e.g. '-4, 0'. "
                                               "Values outside are clipped to the ends.")
        v["legend_num"] = labeled_entry(fig, 7, "Legend # entries", "auto", tip="auto or integer",
                                        hint="Number of size-legend entries (FDR). 'auto' or an integer.")
        v["bubble_size_scale"] = labeled_slider(fig, 8, "Bubble size scale", 0.1, 5.0, 1.0, resolution=0.1,
                                                hint="Multiplier for all circle sizes, legend included "
                                                     "(1 = default, 2 = twice as big, 0.5 = half). Does not change "
                                                     "the FDR values the sizes represent.")
        v["dpi"] = labeled_combo(fig, 9, "DPI", [100, 150, 200, 300, 600], 200,
                                 hint="Resolution of the saved PNG (dots per inch). Higher = sharper, bigger file.")

        opt = ttk.LabelFrame(parent, text="Options")
        opt.pack(fill="x", padx=4, pady=4)
        v["highlight_G_loop"] = labeled_entry(opt, 0, "Highlight G-loop", "0", tip="0 none / 1 5res / 2 8res",
                                              hint="Shade protein labels in the G-loop set: 0 = off, "
                                                   "1 = 5-residue list, 2 = 8-residue list.")
        v["highlight_RT_loop"] = labeled_entry(opt, 1, "Highlight RT-loop", "0", tip="0 none / 1 5res",
                                               hint="Shade protein labels in the RT-loop set: 0 = off, 1 = 5-residue list.")
        v["rainbow_palette"] = check(opt, 2, "Rainbow palette (else coolwarm)", False,
                                     hint="Use a reversed Spectral (rainbow) colour map instead of coolwarm for log2FC.")
        v["invert_xy"] = check(opt, 3, "Invert axes (no dendrogram)", False,
                               hint="Swap the x/y axes (treatments on x, proteins on y) and omit the dendrogram.")
        v["selected_genes"] = labeled_entry(opt, 4, "Selected genes only", "", width=28,
                                            tip="row labels, comma-separated; blank = all",
                                            hint="Restrict the plot to these proteins, written exactly as the row "
                                                 "labels (see 'Protein label' below), comma-separated. Blank = all.")
        v["protein_label_mode"] = labeled_combo(opt, 5, "Protein label",
                                                ["Description | Gene", "Gene only"], "Description | Gene", width=18,
                                                hint="Row label for each protein: the full 'Description | Gene', or just "
                                                     "the gene name.")

        fonts = ttk.LabelFrame(parent, text="Font sizes")
        fonts.pack(fill="x", padx=4, pady=4)
        v["title_fontsize"] = labeled_entry(fonts, 0, "Title", "14", hint="Font size of the plot title.")
        v["axis_fontsize"] = labeled_entry(fonts, 1, "Distance axis", "12",
                                           hint="Font size of the dendrogram 'Distance' axis label.")
        v["compound_labelsize"] = labeled_entry(fonts, 2, "Compound labels", "10",
                                                hint="Font size of the treatment/compound axis labels.")
        v["protein_labelsize"] = labeled_entry(fonts, 3, "Protein labels", "10",
                                               hint="Font size of the protein axis labels.")
        v["colorbar_label_fontsize"] = labeled_entry(fonts, 4, "Colorbar label", "12",
                                                     hint="Font size of the 'Log2FC Value' colorbar label.")
        v["colorbar_tick_fontsize"] = labeled_entry(fonts, 5, "Colorbar ticks", "10",
                                                    hint="Font size of the colorbar tick numbers.")
        v["legend_fontsize"] = labeled_entry(fonts, 6, "Legend", "10",
                                             hint="Font size of the FDR size-legend entries.")

    # ----- tab 5: raw data lookup -----
    def _build_lookup_tab(self, tab):
        bar = ttk.Frame(tab)
        bar.pack(side="bottom", fill="x", padx=4, pady=6)
        self.lookup_btn = ttk.Button(bar, text="Look up", command=self._on_lookup, state="disabled")
        self.lookup_btn.pack(side="left", padx=4)
        Tooltip(self.lookup_btn, "Pull the raw PG (protein-group) and PR (precursor/peptide) "
                                 "intensities for this protein into the table on the right.")

        inner = self._scroll_inner(tab)
        box = ttk.LabelFrame(inner, text="Raw data lookup")
        box.pack(fill="x", padx=4, pady=4)
        self.lookup_query = labeled_entry(box, 0, "UniProt or gene", "", width=22,
                                          hint="A UniProt accession (e.g. P51617) or a gene symbol "
                                               "(e.g. IRAK1). Genes are matched case-insensitively.")
        ttk.Label(box, text="Samples").grid(row=1, column=0, sticky="w", padx=4, pady=2)
        self.lookup_comp = ttk.Combobox(box, state="readonly", width=24, values=[])
        self.lookup_comp.grid(row=1, column=1, sticky="w", padx=4, pady=2)
        Tooltip(self.lookup_comp, "Which samples to show: a comparison (its control + treated samples) "
                                  "or 'All samples'.")
        self.lookup_imputed = check(box, 2, "Also show imputed PG values", False,
                                    hint="Add a 'PG (imputed)' column with the post-imputation protein-group "
                                         "values for the selected comparison. Imputation is per-comparison, so "
                                         "this needs a comparison (not 'All samples'); precursors are never imputed.")
        ttk.Label(inner,
                  text="Shows the raw (pre-imputation) protein-group intensity and each\n"
                       "precursor's intensity for the chosen protein, across the selected\n"
                       "samples. The table appears in the area on the right.",
                  foreground="#666", justify="left").pack(anchor="w", padx=8, pady=6)

    @staticmethod
    def _short_sample(name):
        base = str(name).split("/")[-1].split("\\")[-1]
        for ext in (".mzML", ".raw", ".d", ".dia"):
            if base.lower().endswith(ext.lower()):
                base = base[: -len(ext)]
        return base

    @staticmethod
    def _fmt_val(v):
        try:
            f = float(v)
            return "" if f != f else f"{f:.1f}"   # f != f detects NaN
        except (TypeError, ValueError):
            return "" if v is None else str(v)

    def _resolve_protein(self, query):
        """Resolve a UniProt accession or gene symbol to a Protein.Group index.
        Accepts a ';'/','-joined token (e.g. 'TRAV8-2;TRAV8-4') and matches any part."""
        df = self.result.df_original

        def _one(q):
            q = q.strip()
            if not q:
                return None
            if q in df.index:
                return q
            # UniProt token inside a ';'-joined Protein.Group (e.g. 'A0A0B4J237;P01737')
            for idx in df.index:
                if q in str(idx).split(";"):
                    return idx
            # Gene symbol match (Genes column may be ';'-separated), case-insensitive
            if "Genes" in df.columns:
                ql = q.lower()
                for idx, g in df["Genes"].items():
                    if ql in [t.strip().lower() for t in str(g).split(";")]:
                        return idx
            return None

        hit = _one(query)
        if hit is not None:
            return hit
        for part in str(query).replace(",", ";").split(";"):   # try each member
            hit = _one(part)
            if hit is not None:
                return hit
        return None

    def _on_lookup(self):
        if self.result is None:
            return
        query = self.lookup_query.get().strip()
        if not query:
            messagebox.showerror("Raw data lookup", "Enter a UniProt ID or gene name.")
            return
        pg = self._resolve_protein(query)
        if pg is None:
            messagebox.showinfo("Raw data lookup", f"'{query}' was not found in the protein matrix.")
            return

        gc = self.result.group_columns
        sel = self.lookup_comp.get()
        sample_info = []  # list of (group_name, sample_col)
        comp = None
        if sel and sel != "All samples" and "_vs_" in sel:
            comp = sel
            treated, control = sel.split("_vs_")
            for grp in (control, treated):
                for col in gc.get(grp, []):
                    sample_info.append((grp, col))
        else:
            for grp, cols in gc.items():
                for col in cols:
                    sample_info.append((grp, col))

        want_imputed = bool(self.lookup_imputed.get())
        if want_imputed and comp is None:
            messagebox.showinfo("Raw data lookup",
                                "Imputed values are per-comparison — pick a comparison (not 'All samples') "
                                "to also show the imputed PG column.")
        try:
            self._show_lookup(pg, sample_info, comp if want_imputed else None)
            logging.getLogger().info("Raw lookup: %s (%s) over %d samples", query, pg, len(sample_info))
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            messagebox.showerror("Raw data lookup", "Lookup failed. See the Log.")

    def _show_lookup(self, pg, sample_info, imputed_comp=None):
        df = self.result.df_original
        dpr = self.result.df_peptide

        # Optional imputed PG values (per-comparison; only if the protein survived).
        imp_df = None
        if imputed_comp and imputed_comp in self.result.imputed_dataframes:
            cand = self.result.imputed_dataframes[imputed_comp]
            if pg in cand.index:
                imp_df = cand

        genes = df.loc[pg].get("Genes") if "Genes" in df.columns else None
        desc = df.loc[pg].get("First.Protein.Description") if "First.Protein.Description" in df.columns else ""

        # Precursor rows for this protein (PR matrix is indexed by Protein.Group).
        if pg in dpr.index:
            prec = dpr.loc[[pg]]
        else:
            prec = dpr.iloc[0:0]
        label_col = next((c for c in ("Precursor.Id", "Modified.Sequence", "Stripped.Sequence")
                          if c in dpr.columns), None)
        prec_labels = [str(prec.iloc[i][label_col]) if label_col else f"precursor {i + 1}"
                       for i in range(len(prec))]

        # Render into the lookup tab's cached output frame.
        area = self._tab_area("lookup")
        self._plot_canvas = None
        ttk.Label(area,
                  text=f"Protein.Group: {pg}    Genes: {genes}\n{desc}    "
                       f"({len(prec)} precursor{'s' if len(prec) != 1 else ''})",
                  justify="left", foreground="#222").pack(anchor="w", padx=8, pady=(6, 2))

        container = ttk.Frame(area)
        container.pack(side="top", fill="both", expand=True, padx=4, pady=4)

        cols = ["group", "sample", "pg"] + (["pgimp"] if imp_df is not None else []) \
            + [f"p{i}" for i in range(len(prec))]
        tree = ttk.Treeview(container, columns=cols, show="headings")
        vsb = ttk.Scrollbar(container, orient="vertical", command=tree.yview)
        hsb = ttk.Scrollbar(container, orient="horizontal", command=tree.xview)
        tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        tree.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="ew")
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)

        tree.heading("group", text="Group"); tree.column("group", width=130, anchor="w")
        tree.heading("sample", text="Sample"); tree.column("sample", width=240, anchor="w")
        tree.heading("pg", text="PG (raw)" if imp_df is not None else "PG")
        tree.column("pg", width=100, anchor="e")
        if imp_df is not None:
            tree.heading("pgimp", text="PG (imputed)"); tree.column("pgimp", width=110, anchor="e")
        for i, lab in enumerate(prec_labels):
            tree.heading(f"p{i}", text=lab)
            tree.column(f"p{i}", width=130, anchor="e")

        # Give each group its own light background so the groups are easy to tell apart.
        palette = ["#eaf3ff", "#fff0e6", "#eafbea", "#f3eaff", "#fffbe6",
                   "#ffeaea", "#e6fffb", "#f2f2f2"]
        group_tag = {}
        for grp, _col in sample_info:
            if grp not in group_tag:
                tag = f"g{len(group_tag)}"
                group_tag[grp] = tag
                tree.tag_configure(tag, background=palette[(len(group_tag) - 1) % len(palette)])

        imp_row = imp_df.loc[pg] if imp_df is not None else None
        for grp, col in sample_info:
            pg_val = df.loc[pg].get(col)
            row = [grp, self._short_sample(col), self._fmt_val(pg_val)]
            if imp_df is not None:
                row.append(self._fmt_val(imp_row.get(col)))
            row += [self._fmt_val(prec.iloc[i].get(col)) for i in range(len(prec))]
            tree.insert("", "end", tags=(group_tag[grp],), values=row)

        self._show_tab_area("lookup")

    # ----- data file browsing -----
    def _browse_pg(self):
        path = filedialog.askopenfilename(
            title="Select the protein groups matrix (<stem>.pg_matrix.tsv)",
            filetypes=[("DIA-NN protein matrix", "*.pg_matrix.tsv"),
                       ("TSV files", "*.tsv"), ("All files", "*.*")])
        if not path:
            return
        path = os.path.abspath(path)
        self.cfg_vars["pg_path"].set(path)
        # Auto-fill the matching pr_matrix with the same stem in the same folder.
        candidate = os.path.join(os.path.dirname(path), _stem_from_pg(path) + ".pr_matrix.tsv")
        if os.path.exists(candidate):
            self.cfg_vars["pr_path"].set(candidate)
        self._update_path_label()

    def _browse_pr(self):
        path = filedialog.askopenfilename(
            title="Select the precursor/peptide matrix (<stem>.pr_matrix.tsv)",
            filetypes=[("DIA-NN precursor matrix", "*.pr_matrix.tsv"),
                       ("TSV files", "*.tsv"), ("All files", "*.*")])
        if not path:
            return
        self.cfg_vars["pr_path"].set(os.path.abspath(path))
        self._update_path_label()

    def _output_dir_for(self, pg_path):
        workdir = os.path.dirname(os.path.abspath(pg_path))
        return os.path.join(workdir, "proteomics_GUI_output")

    # ----- workspace: save / recall all settings -----
    def _collect_workspace(self):
        """Snapshot every input setting as a JSON-serializable dict, plus small
        structural info (group_columns, comparisons) so a recall can rebuild the
        result from the saved outputs without re-deriving anything."""
        def vals(d):
            return {k: var.get() for k, var in d.items()}
        ws = {
            "app": "proteomics_GUI",
            "version": _app_version(),
            "config": vals(self.cfg_vars),
            "group_patterns": self._group_patterns or {},
            "drop_samples": list(self._drop_samples),
            "pca": vals(self.pca_vars),
            "volcano": vals(self.vol_vars),
            "bubble": vals(self.bub_vars),
            "sar_text": self._sar_text.get("1.0", "end").rstrip("\n") if self._sar_text else "",
            "treat": self.treat_cb.get(),
            "control": self.ctrl_cb.get(),
            "lookup_query": self.lookup_query.get(),
            "lookup_comp": self.lookup_comp.get(),
        }
        if self.result is not None:
            ws["group_columns"] = {k: list(v) for k, v in self.result.group_columns.items()}
            ws["comparisons"] = list(self.result.imputed_dataframes.keys())
        return ws

    def _apply_workspace(self, ws):
        """Restore settings from a workspace dict (only keys we recognise)."""
        def setv(d, data):
            for k, value in (data or {}).items():
                if k in d:
                    try:
                        d[k].set(value)
                    except Exception:
                        pass
        setv(self.cfg_vars, ws.get("config"))
        self._group_patterns = dict(ws.get("group_patterns") or {})
        self._set_drop_samples(ws.get("drop_samples"))
        setv(self.pca_vars, ws.get("pca"))
        setv(self.vol_vars, ws.get("volcano"))
        setv(self.bub_vars, ws.get("bubble"))
        if self._sar_text is not None and "sar_text" in ws:
            self._sar_text.delete("1.0", "end")
            self._sar_text.insert("1.0", ws.get("sar_text", ""))
        for cb, key in ((self.treat_cb, "treat"), (self.ctrl_cb, "control"),
                        (self.lookup_comp, "lookup_comp")):
            val = ws.get(key)
            if val:
                try:
                    cb.set(val)
                except Exception:
                    pass
        if "lookup_query" in ws:
            self.lookup_query.set(ws.get("lookup_query", ""))
        self._update_path_label()

    def _write_workspace(self, path):
        """Write the workspace as small JSON: settings + pointers to the data files
        (no copies of the matrices). Atomic (temp + replace)."""
        path = os.path.abspath(path)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        tmp = f"{path}.{os.getpid()}.tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(self._collect_workspace(), f, indent=2)
        os.replace(tmp, path)

    def _autosave_workspace(self, path):
        """Auto-save the (small) workspace JSON; logged, best-effort."""
        try:
            self._write_workspace(path)
            logging.getLogger().info("Workspace saved -> %s", path)
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")

    def _workspace_dir(self):
        """Best guess for where workspaces live: the output folder."""
        if self.output_dir:
            return self.output_dir
        pg = self.cfg_vars["pg_path"].get().strip()
        if pg:
            return self._output_dir_for(pg)
        return os.path.dirname(os.path.abspath(__file__))

    @staticmethod
    def _results_xlsx_name(config):
        """The results-Excel filename _save_workflow_excel would have produced."""
        imp = bool(getattr(config, "imputation_option", True))
        norm = bool(getattr(config, "normalization_protein_id", ""))
        if imp and norm:
            return "final_analysis_results_imputed_normalized.xlsx"
        if imp:
            return "final_analysis_results_imputed.xlsx"
        if norm:
            return "final_analysis_results_normalized.xlsx"
        return "final_analysis_results.xlsx"

    def _reconstruct_from_outputs(self, ws):
        """Rebuild an AnalysisResult from the saved results Excel + the workspace
        config (no limma re-run). Returns the result, or None if it can't be done
        (then the caller falls back to re-running)."""
        import pandas as pd
        from scripts.pipeline import AnalysisResult
        try:
            cfg = build_config(self._cfg_values())
        except Exception:
            return None
        pg = (getattr(cfg, "pg_path", "") or "").strip()
        pr = (getattr(cfg, "pr_path", "") or "").strip()
        if not pg:
            return None
        xlsx = os.path.join(self._output_dir_for(pg), self._results_xlsx_name(cfg))
        if not os.path.exists(xlsx):
            return None

        def _unrename(c):  # adjPvalue_<comp> -> bh_FDR_<comp> (internal name)
            c = str(c)
            return ("bh_FDR_" + c[len("adjPvalue_"):]) if c.startswith("adjPvalue_") else c

        def _is_true(v):
            return v is True or v == 1 or str(v).strip().upper() == "TRUE"

        xls = pd.ExcelFile(xlsx)
        summary = xls.parse("Fold_Change_Summary").rename(columns=_unrename)
        if "Protein.Group" in summary.columns:
            summary = summary.set_index("Protein.Group")

        comparisons = ws.get("comparisons") or [str(c)[len("log2FC_"):]
                                                for c in summary.columns if str(c).startswith("log2FC_")]
        group_columns = ws.get("group_columns") or {}

        stat_prefixes = ("FC_", "log2FC_", "Pvalue_", "bh_FDR_", "adjPvalue_", "-log_P_adj_")
        keep = [c for c in summary.columns if not str(c).startswith(stat_prefixes)]
        df_original = summary[keep].copy()           # metadata + sample intensities

        # An older workspace may not have stored group_columns (saved before the
        # analysis finished). Rebuild it from the config's group patterns matched
        # against the sample-intensity columns so plotting/imputation can work.
        if not group_columns and getattr(cfg, "group_names", None):
            from scripts.io import assign_groups
            group_columns = assign_groups(df_original, cfg.group_names,
                                          getattr(cfg, "group_patterns", None))

        imputed_dataframes, imputation_dict = {}, {}
        for comp in comparisons:
            sheet = comp[:31]
            if sheet not in xls.sheet_names:
                continue
            d = xls.parse(sheet).rename(columns=_unrename)
            if "Protein.Group" in d.columns:
                d = d.set_index("Protein.Group")
            imputation_dict[comp] = [i for i, v in d["Imputed"].items() if _is_true(v)] \
                if "Imputed" in d.columns else []
            imputed_dataframes[comp] = d

        # Precursor matrix for the lookup (read from the pointed file; empty if absent).
        try:
            df_peptide = pd.read_csv(pr, sep="\t", index_col=0) if pr and os.path.exists(pr) else pd.DataFrame()
        except Exception:
            df_peptide = pd.DataFrame()

        return AnalysisResult(config=cfg, df_original=df_original, df_peptide=df_peptide,
                              group_columns=group_columns, imputed_dataframes=imputed_dataframes,
                              imputation_dict=imputation_dict, summary=summary)

    def _activate_result(self, result, status="Loaded from saved outputs"):
        """Put a (reconstructed or computed) result into the GUI and show its plots."""
        self.cfg = result.config
        self.result = result
        try:
            if getattr(result.config, "pg_path", ""):
                self.output_dir = self._output_dir_for(result.config.pg_path)
        except Exception:
            pass
        self._set_group_text(result.group_columns)
        groups = list(result.group_columns.keys())
        self.treat_cb.configure(values=groups)
        self.ctrl_cb.configure(values=groups)
        if not self.treat_cb.get():
            if result.config.comparison_matrix:
                self.treat_cb.set(result.config.comparison_matrix[0][0])
                self.ctrl_cb.set(result.config.comparison_matrix[0][1])
            elif groups:
                self.treat_cb.set(groups[-1])
                self.ctrl_cb.set(result.config.reference_group)
        comps = list(result.imputed_dataframes.keys())
        self.lookup_comp.configure(values=comps + ["All samples"])
        if not self.lookup_comp.get():
            self.lookup_comp.set(comps[0] if comps else "All samples")
        for b in (self.plot_all_btn, self.plot_btn, self.pca_btn, self.bubble_btn, self.lookup_btn):
            b.configure(state="normal")
        self.status.configure(text=status, foreground="#080")
        self.nb.select(2)                  # Volcano tab
        self.after(50, self._on_plot_all)

    def _on_load_workspace(self):
        path = filedialog.askopenfilename(
            title="Load workspace", initialdir=self._workspace_dir(),
            filetypes=[("Workspace", "*.json"), ("All files", "*.*")])
        if not path:
            return
        try:
            with open(path, "r", encoding="utf-8") as f:
                ws = json.load(f)
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            messagebox.showerror("Workspace", "Could not read the workspace file. See the Log.")
            return
        self._apply_workspace(ws)
        # Prefer rebuilding from the saved results Excel (instant, no R needed).
        result = None
        try:
            result = self._reconstruct_from_outputs(ws)
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            result = None
        if result is not None:
            self._activate_result(result)
            logging.getLogger().info("Workspace reopened from saved outputs: %s", path)
            return
        # Fall back: re-run from the data files the workspace points to.
        pg = self.cfg_vars["pg_path"].get().strip()
        if pg and os.path.exists(pg):
            self._pending_workspace = ws
            logging.getLogger().info("No saved outputs; re-running from %s", pg)
            self._on_run()
        else:
            messagebox.showwarning(
                "Workspace",
                "Settings loaded, but there are no saved outputs to reopen and the protein file path "
                f"doesn't exist:\n{pg or '(none)'}\n\nPick the file on this machine, then click 'Run analysis'.")

    def _update_path_label(self):
        pg = self.cfg_vars["pg_path"].get().strip()
        if pg:
            workdir = os.path.dirname(os.path.abspath(pg))
            self.path_label.configure(
                text=f"Working folder: {workdir}\nOutputs -> {self._output_dir_for(pg)}")
        else:
            self.path_label.configure(text="No data selected.")

    def _prefill_sample_data(self):
        """If sample *.pg_matrix.tsv files ship next to this script, preselect one."""
        here = os.path.dirname(os.path.abspath(__file__))
        pgs = sorted(glob.glob(os.path.join(here, "*.pg_matrix.tsv")))
        if not pgs:
            self._update_path_label()
            return
        pg = pgs[0]
        self.cfg_vars["pg_path"].set(pg)
        candidate = os.path.join(here, _stem_from_pg(pg) + ".pr_matrix.tsv")
        if os.path.exists(candidate):
            self.cfg_vars["pr_path"].set(candidate)
        self._update_path_label()

    def _ensure_outdir_cwd(self):
        """Make the dedicated outputs folder the current dir so all writers land there."""
        if self.output_dir:
            os.makedirs(self.output_dir, exist_ok=True)
            os.chdir(self.output_dir)

    # ----- collecting values -----
    def _cfg_values(self):
        v = {k: var.get() for k, var in self.cfg_vars.items()}
        v["group_patterns"] = self._group_patterns or None
        v["drop_samples"] = list(self._drop_samples)
        return v

    def _vol_values(self):
        return {k: var.get() for k, var in self.vol_vars.items()}

    def _bub_values(self):
        vals = {k: var.get() for k, var in self.bub_vars.items()}
        vals["sar"] = self._sar_text.get("1.0", "end")
        return vals

    # ----- run analysis (stage 1) -----
    def _on_run(self):
        try:
            cfg = build_config(self._cfg_values())
        except Exception as e:
            messagebox.showerror("Configuration error", str(e))
            return
        self.output_dir = self._output_dir_for(cfg.pg_path)
        try:
            os.makedirs(self.output_dir, exist_ok=True)
        except Exception as e:
            messagebox.showerror("Output folder error",
                                 f"Could not create outputs folder:\n{self.output_dir}\n\n{e}")
            return
        # Start a fresh log file for this run in the outputs folder; everything
        # shown in the Log box is mirrored to it (see _drain_log).
        log_path = os.path.join(self.output_dir, "analysis_log.txt")
        try:
            if self._log_fh is not None:
                self._log_fh.close()
            self._log_fh = open(log_path, "w", encoding="utf-8")
        except Exception:
            self._log_fh = None
        for b in (self.run_btn, self.preview_btn, self.plot_all_btn, self.plot_btn,
                  self.pca_btn, self.bubble_btn, self.lookup_btn):
            b.configure(state="disabled")
        self.status.configure(text="Running...", foreground="#a60")
        log = logging.getLogger()
        log.info("Working folder: %s", os.path.dirname(cfg.pg_path))
        log.info("Outputs folder: %s", self.output_dir)
        if self._log_fh is not None:
            log.info("Log file: %s", log_path)
        log.info("Starting analysis for '%s'...", cfg.file)
        # rpy2's R<->pandas conversion rules live in a contextvars Context that a
        # freshly spawned thread does NOT inherit -- which made limma work on the
        # first run but fail on the second (new thread, empty context, and the
        # one-time activate() is skipped). Initialise R here on the main thread so
        # the rules persist in this context, then run each worker inside a COPY of
        # it so every run sees them.
        try:
            from scripts.stats import _ensure_r
            _ensure_r()
        except Exception:
            pass  # a genuine R/limma error will be logged by the worker
        ctx = contextvars.copy_context()
        threading.Thread(target=ctx.run, args=(self._run_worker, cfg), daemon=True).start()

    def _run_worker(self, cfg):
        old_out, old_err = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = _QueueWriter(self.log_q)
        try:
            os.makedirs(self.output_dir, exist_ok=True)
            os.chdir(self.output_dir)
            result = run_core(cfg)
            self.result_q.put(("ok", cfg, result))
        except Exception:
            self.result_q.put(("error", cfg, traceback.format_exc()))
        finally:
            sys.stdout, sys.stderr = old_out, old_err

    def _poll_result(self):
        try:
            kind, cfg, payload = self.result_q.get_nowait()
        except queue.Empty:
            self.after(150, self._poll_result)
            return

        self.run_btn.configure(state="normal")
        self.preview_btn.configure(state="normal")
        if kind == "error":
            self.status.configure(text="Failed", foreground="#c00")
            self.log_q.put(payload + "\n")
            messagebox.showerror("Analysis failed", "See the log for details.")
        else:
            self.cfg = cfg
            self.result = payload
            self._set_group_text(payload.group_columns)
            groups = list(payload.group_columns.keys())
            self.treat_cb.configure(values=groups)
            self.ctrl_cb.configure(values=groups)
            if cfg.comparison_matrix:
                self.treat_cb.set(cfg.comparison_matrix[0][0])
                self.ctrl_cb.set(cfg.comparison_matrix[0][1])
            elif groups:
                self.treat_cb.set(groups[-1])
                self.ctrl_cb.set(cfg.reference_group)
            self.vol_vars["imputation_option"].set(cfg.imputation_option)
            # Prefill bubble SAR with the treatment groups vs the reference.
            treatments = [k for k in groups if k != cfg.reference_group]
            self._sar_text.delete("1.0", "end")
            self._sar_text.insert("1.0", ": " + ", ".join(treatments))
            self.bub_vars["sar_suffix"].set("_vs_" + cfg.reference_group)
            # Raw-data lookup: offer the comparisons (or all samples).
            comps = list(payload.imputed_dataframes.keys())
            self.lookup_comp.configure(values=comps + ["All samples"])
            self.lookup_comp.set(comps[0] if comps else "All samples")
            for b in (self.plot_all_btn, self.plot_btn, self.pca_btn, self.bubble_btn, self.lookup_btn):
                b.configure(state="normal")
            self.status.configure(text="Analysis complete", foreground="#080")
            # If this run came from loading a workspace, re-apply its settings now
            # (the wiring above reset a few to run-derived defaults).
            if getattr(self, "_pending_workspace", None):
                self._apply_workspace(self._pending_workspace)
                self._pending_workspace = None
                self.status.configure(text="Loaded from workspace", foreground="#080")
            # Auto-save the (small) workspace -- settings + pointers to the data
            # files -- so the session can be recalled later by re-running from them.
            self._autosave_workspace(os.path.join(self.output_dir, "workspace.json"))
            # By default, generate volcanoes for ALL comparisons and show that tab.
            self.nb.select(2)  # Volcano tab
            self.after(50, self._on_plot_all)
        self.after(150, self._poll_result)

    # ----- plotting (renders into the active tab's cached output frame) -----
    def _embed(self, fig, key):
        container = self._tab_area(key)
        canvas = FigureCanvasTkAgg(fig, master=container)
        canvas.draw()
        NavigationToolbar2Tk(canvas, container).update()
        canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        self._plot_canvas = canvas
        self._show_tab_area(key)

    def _embed_notebook(self, items, key="volcano"):
        """Embed several figures, one inner tab per (name, figure)."""
        container = self._tab_area(key)
        nb = ttk.Notebook(container)
        nb.pack(side="top", fill="both", expand=True)
        if key == "volcano":
            self._volcano_nb = nb            # remember it to read the active comparison
        self._tab_canvases = []
        for name, fig in items:
            tab = ttk.Frame(nb)
            nb.add(tab, text=name)
            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.draw()
            NavigationToolbar2Tk(canvas, tab).update()
            canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
            self._tab_canvases.append(canvas)
        self._plot_canvas = self._tab_canvases[-1] if self._tab_canvases else None
        self._show_tab_area(key)

    def _embed_image(self, path, key):
        """Show a saved PNG (scaled to fit the panel, aspect preserved) instead of
        a live re-render -- so the preview matches the exact saved figure."""
        container = self._tab_area(key)
        holder = ttk.Frame(container)
        holder.pack(side="top", fill="both", expand=True)
        lbl = ttk.Label(holder, anchor="center")
        lbl.pack(fill="both", expand=True)
        try:
            base = Image.open(path); base.load()
        except Exception:
            lbl.configure(text=f"Could not load image:\n{path}")
            self._plot_canvas = None
            self._show_tab_area(key)
            return
        last = [0, 0]

        def _fit(_evt=None):
            W, H = holder.winfo_width(), holder.winfo_height()
            if W <= 1 or H <= 1:
                return
            iw, ih = base.size
            s = min(W / iw, H / ih)
            nw, nh = max(1, int(iw * s)), max(1, int(ih * s))
            if abs(nw - last[0]) < 2 and abs(nh - last[1]) < 2:
                return
            last[0], last[1] = nw, nh
            photo = ImageTk.PhotoImage(base.resize((nw, nh), Image.LANCZOS))
            lbl.configure(image=photo); lbl.image = photo

        holder.bind("<Configure>", _fit)
        self.after(10, _fit)
        self._plot_canvas = None   # an image view has no matplotlib canvas
        self._show_tab_area(key)

    def _point_hover(self, canvas, xs, ys, labels, on_pick=None):
        """Generic: show labels[i] in a tooltip when hovering near point (xs[i], ys[i]).
        If on_pick is given, a left-click near a point calls on_pick(i)."""
        if canvas is None:
            return
        try:
            fig = canvas.figure
            if not fig.axes:
                return
            ax = fig.axes[0]
            xs = np.asarray(xs, dtype=float)
            ys = np.asarray(ys, dtype=float)
            if xs.size == 0:
                return
            labels = [str(l) for l in labels]
            pts = np.column_stack([xs, ys])
            annot = ax.annotate("", xy=(0, 0), xytext=(12, 12), textcoords="offset points",
                                bbox=dict(boxstyle="round", fc="#ffffe0", ec="0.5", alpha=0.95),
                                fontsize=9, zorder=20)
            annot.set_visible(False)

            def _nearest(event):
                disp = ax.transData.transform(pts)
                d = np.hypot(disp[:, 0] - event.x, disp[:, 1] - event.y)
                i = int(np.argmin(d))
                return (i, d[i])

            def _on_move(event):
                if event.inaxes is not ax:
                    if annot.get_visible():
                        annot.set_visible(False); canvas.draw_idle()
                    return
                i, dist = _nearest(event)
                if dist <= 12:                       # within ~12 px of a point
                    annot.xy = (xs[i], ys[i])
                    annot.set_text(labels[i])
                    annot.set_visible(True)
                    canvas.draw_idle()
                elif annot.get_visible():
                    annot.set_visible(False); canvas.draw_idle()

            canvas.mpl_connect("motion_notify_event", _on_move)

            if on_pick is not None:
                def _on_click(event):
                    if event.inaxes is not ax or getattr(event, "button", None) != 1:
                        return
                    tb = getattr(canvas, "toolbar", None)   # ignore clicks while zoom/pan is active
                    if tb is not None and getattr(tb, "mode", ""):
                        return
                    i, dist = _nearest(event)
                    if dist <= 12:
                        on_pick(i)
                canvas.mpl_connect("button_press_event", _on_click)
        except Exception:
            pass

    def _attach_hover(self, canvas, treated, control):
        """Hover a volcano dot -> show its gene name."""
        if self.result is None or canvas is None:
            return
        df = self.result.summary
        xcol = f"log2FC_{treated}_vs_{control}"
        fcol = (f"bh_FDR_{treated}_vs_{control}" if getattr(self.cfg, "output_adjpval", True)
                else f"Pvalue_{treated}_vs_{control}")
        if xcol not in df.columns or fcol not in df.columns:
            return
        sub = df[[xcol, fcol, "Genes"]].copy()
        sub = sub[sub[xcol].notna() & sub[fcol].notna() & (sub[fcol] > 0)]
        if sub.empty:
            return
        labels = [g if str(g).strip() else a
                  for g, a in zip(sub["Genes"].astype(str), sub.index.astype(str))]
        accs = sub.index.astype(str).tolist()

        def _on_pick(i):
            # Prefer the gene label for display; fall back to the accession.
            self._lookup_from_volcano(labels[i] or accs[i], treated, control)

        self._point_hover(canvas, sub[xcol].to_numpy(dtype=float),
                          -np.log10(sub[fcol].to_numpy(dtype=float)), labels, on_pick=_on_pick)

    def _lookup_from_volcano(self, query, treated, control):
        """Click a volcano dot -> jump to the lookup tab, prefilled, and run it."""
        comp = f"{treated}_vs_{control}"
        try:
            if comp in list(self.lookup_comp.cget("values")):
                self.lookup_comp.set(comp)
        except Exception:
            pass
        self.lookup_query.set(str(query))
        self.nb.select(4)                 # step 5: Raw data lookup
        self.update_idletasks()
        self._on_lookup()

    def _attach_hover_pca(self, canvas, pca_df):
        """Hover a PCA dot -> show its sample name; click -> toggle it in step 1's
        drop-off list (takes effect on the next run)."""
        if canvas is None or pca_df is None or "PC1" not in getattr(pca_df, "columns", []):
            return
        samples = [str(s) for s in pca_df.index]          # full column names (match the drop list)
        labels = [s.split("/")[-1] for s in samples]      # short names for the tooltip
        self._point_hover(canvas, pca_df["PC1"].to_numpy(dtype=float),
                          pca_df["PC2"].to_numpy(dtype=float), labels,
                          on_pick=lambda i: self._drop_toggle_from_pca(samples[i]))

    def _drop_toggle_from_pca(self, sample):
        """Click a PCA dot: add its sample to the drop-off list, or remove it if
        already there. Logged; excluded from the analysis on the next run."""
        sample = str(sample)
        if sample in self._drop_samples:
            idx = self._drop_samples.index(sample)
            self._drop_samples.pop(idx)
            self.drop_list.delete(idx)
            self._update_drop_hint()
            logging.getLogger().info("Removed '%s' from the drop-off list.", sample)
        else:
            self._drop_samples.append(sample)
            self.drop_list.insert("end", sample)
            self._update_drop_hint()
            logging.getLogger().info("Added '%s' to the drop-off list; it will be "
                                     "excluded on the next run.", sample)

    def _on_plot_all(self):
        """Generate a volcano for every comparison using the current settings."""
        if self.result is None:
            return
        try:
            params = build_volcano_params(self._vol_values())
        except Exception as e:
            messagebox.showerror("Invalid volcano setting", str(e))
            return

        comparisons = list(self.result.imputed_dataframes.keys())  # 'treated_vs_control'
        if not comparisons:
            messagebox.showinfo("Volcano", "No comparisons were produced by the analysis.")
            return

        self._ensure_outdir_cwd()
        plt.close("all")
        items, errors = [], []
        for name in comparisons:
            treated_name, control_name = name.split("_vs_")
            try:
                volcano_plot(
                    treated_name, control_name,
                    df=self.result.summary,
                    group_columns=self.result.group_columns,
                    imputation_dict=self.result.imputation_dict,
                    config=self.cfg,
                    **params,
                )
                items.append((name, plt.gcf()))
            except Exception:
                errors.append(name)
                self.log_q.put(f"[{name}] {traceback.format_exc()}\n")

        if items:
            self._embed_notebook(items)
            for (name, _fig), canvas in zip(items, self._tab_canvases):
                t, c = name.split("_vs_")
                self._attach_hover(canvas, t, c)
        logging.getLogger().info("Plotted %d/%d comparisons (PNGs saved).", len(items), len(comparisons))
        if errors:
            messagebox.showwarning("Some comparisons failed",
                                   "Failed: " + ", ".join(errors) + "\nSee the log for details.")

    def _on_plot_volcano(self):
        if self.result is None:
            return
        treated = self.treat_cb.get()
        control = self.ctrl_cb.get()
        if not treated or not control:
            messagebox.showerror("Volcano", "Select a treatment and control group.")
            return
        try:
            params = build_volcano_params(self._vol_values())
        except Exception as e:
            messagebox.showerror("Invalid volcano setting", str(e))
            return
        try:
            self._ensure_outdir_cwd()
            plt.close("all")
            volcano_plot(
                treated, control,
                df=self.result.summary,
                group_columns=self.result.group_columns,
                imputation_dict=self.result.imputation_dict,
                config=self.cfg,
                **params,
            )
            self._embed(plt.gcf(), "volcano")          # single plot replaces any notebook
            self._volcano_nb = None
            self._current_volcano_comp = f"{treated}_vs_{control}"
            self._attach_hover(self._plot_canvas, treated, control)
            logging.getLogger().info("Plotted volcano %s vs %s", treated, control)
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            messagebox.showerror("Plot failed", "See the log for details.")

    def _on_plot_pca(self):
        if self.result is None:
            return
        try:
            self._ensure_outdir_cwd()
            plt.close("all")
            pca_df = generate_pca_plot(
                self.result.df_original, self.result.group_columns,
                filename=str(self.pca_vars["filename"].get()).strip() or "PCA_plot.png",
                title=str(self.pca_vars["title"].get()),
                text=bool(self.pca_vars["text"].get()),
                title_fontsize=_req_float(self.pca_vars["title_fontsize"].get(), 20),
                axis_fontsize=_req_float(self.pca_vars["axis_fontsize"].get(), 15),
                tick_fontsize=_opt_float(self.pca_vars["tick_fontsize"].get()),
                legend_fontsize=_opt_float(self.pca_vars["legend_fontsize"].get()),
                point_label_fontsize=_req_float(self.pca_vars["point_label_fontsize"].get(), 4),
                dpi=_req_int(self.pca_vars["dpi"].get(), 300),
            )
            self._embed(plt.gcf(), "pca")
            self._attach_hover_pca(self._plot_canvas, pca_df)
            logging.getLogger().info("Plotted PCA")
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            messagebox.showerror("PCA failed", "See the log for details.")

    def _on_plot_bubble(self):
        if self.result is None:
            return
        try:
            SAR, kwargs = build_bubble_params(self._bub_values())
        except Exception as e:
            messagebox.showerror("Invalid bubble setting", str(e))
            return
        try:
            self._ensure_outdir_cwd()
            # Pass the in-memory summary straight to the bubble plot (no CSV round-trip).
            plt.close("all")
            bubble_dendro_plot(SAR, self.cfg, df=self.result.summary, **kwargs)
            # Show the exact saved PNG (preserves the requested figure aspect, unlike
            # a fit-to-panel re-render). Fall back to a live embed if Pillow is absent.
            png = os.path.abspath(kwargs["figure_filename"])
            if _HAS_PIL and os.path.exists(png):
                plt.close("all")
                self._embed_image(png, "bubble")
            else:
                self._embed(plt.gcf(), "bubble")
            logging.getLogger().info("Plotted bubble plot -> %s", kwargs["figure_filename"])
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            messagebox.showerror(
                "Bubble plot failed",
                "See the log for details.\n\nTip: the bubble plot needs at least two "
                "proteins that are significantly down-regulated (log2FC < -1, bh_FDR < 0.01) "
                "across the SAR treatments.")


def main():
    app = VolcanoGUI()
    app.mainloop()


if __name__ == "__main__":
    main()
