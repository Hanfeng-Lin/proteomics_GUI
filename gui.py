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
import sys
import glob
import queue
import logging
import threading
import traceback
import contextvars

import tkinter as tk
from tkinter import ttk, messagebox, filedialog

import matplotlib
matplotlib.use("Agg")  # figures are embedded manually; no stray windows
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Make the package importable when running this file directly.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scripts import AnalysisConfig, run_core
from scripts.plots.volcano import volcano_plot
from scripts.plots.pca import generate_pca_plot
from scripts.plots.bubble import bubble_dendro_plot


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
        compound_labelsize=_req_float(v["compound_labelsize"], 20),
        protein_labelsize=_req_float(v["protein_labelsize"], 20),
        colorFCrange=_num_pair(v["colorFCrange"], "color FC range"),
        highlight_G_loop=int(float(str(v["highlight_G_loop"]).strip() or 0)),
        highlight_RT_loop=int(float(str(v["highlight_RT_loop"]).strip() or 0)),
        rainbow_palette=1 if bool(v["rainbow_palette"]) else 0,
        invert_xy=bool(v["invert_xy"]),
        selected_genes=_genes(v["selected_genes"]),
        legend_num=legend_num,
        title_fontsize=_req_float(v["title_fontsize"], 30),
        axis_fontsize=_req_float(v["axis_fontsize"], 30),
        colorbar_label_fontsize=_req_float(v["colorbar_label_fontsize"], 30),
        colorbar_tick_fontsize=_req_float(v["colorbar_tick_fontsize"], 30),
        legend_fontsize=_req_float(v["legend_fontsize"], 30),
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
    groups = [g.strip() for g in str(v["group_names"]).split(",") if g.strip()]

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


def group_candidates(sample_cols):
    """Propose group-name schemes from sample column headers.

    Returns (prefix, suffix, anchors, candidates):
      - anchors    : every consensus region shared by all sample names.
      - candidates : list of (description, [group names]), variable fields in
        header order first (earlier fields, usually the compound, ranked first),
        then a delimiter fallback. Handles MULTIPLE consensus regions and group
        names that contain the delimiter (e.g. 'Positive_Control').
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

    tokens = _decompose(cores)
    anchors = [t[1] for t in tokens if t[0] == "anchor" and t[1].strip()]

    raw = []  # (priority, description, groups)
    pos = 0
    for idx, t in enumerate(tokens):
        if t[0] != "seg":
            continue
        vals = _uniq(t[1])
        if len(vals) <= 1:  # constant or empty segment -> part of the consensus
            continue
        nxt = next((tokens[j][1] for j in range(idx + 1, len(tokens)) if tokens[j][0] == "anchor"), "")
        prv = next((tokens[j][1] for j in range(idx - 1, -1, -1) if tokens[j][0] == "anchor"), "")
        if nxt:
            desc = f"Field {pos + 1} (before '{nxt}')"
        elif prv:
            desc = f"Field {pos + 1} (after '{prv}')"
        else:
            desc = f"Field {pos + 1}"
        raw.append((pos, desc, vals))
        pos += 1

    # Delimiter fallback (lower priority) for names with no usable consensus.
    for delim in ("_", "-", "."):
        if not any(delim in core for core in cores):
            continue
        first = _uniq([core.split(delim, 1)[0] for core in cores])
        if 1 < len(first):
            raw.append((100, f"First part before '{delim}'", first))

    # Dedup by group set, keeping the best (lowest-priority) description.
    best = {}
    for pri, desc, g in raw:
        key = tuple(sorted(g))
        if key not in best or pri < best[key][0]:
            best[key] = (pri, desc, g)
    ranked = sorted(best.values(), key=lambda x: x[0])
    return prefix, suffix, anchors, [(desc, g) for _, desc, g in ranked]


class GroupPickerDialog(tk.Toplevel):
    """Popup showing the header consensus and candidate group fields to choose."""
    def __init__(self, parent, sample_cols, apply_cb):
        super().__init__(parent)
        self.title("Auto-pick group names")
        self.geometry("680x500")
        self.minsize(560, 360)
        self.apply_cb = apply_cb
        self.transient(parent)

        prefix, suffix, anchors, self.cands = group_candidates(sample_cols)
        anchors_str = "   |   ".join(f"'{a}'" for a in anchors) if anchors else "(none)"
        info = (f"{len(sample_cols)} sample columns detected.\n"
                f"Common prefix:  '{prefix}'\n"
                f"Common suffix:  '{suffix}'\n"
                f"Consensus regions:  {anchors_str}")
        ttk.Label(self, text=info, justify="left", foreground="#444").pack(anchor="w", padx=10, pady=8)
        ttk.Label(self, text="Choose which field defines your groups:").pack(anchor="w", padx=10)

        # Reserve the button bar at the bottom FIRST so it is always visible,
        # then let the body fill the space above it.
        btns = ttk.Frame(self)
        btns.pack(side="bottom", fill="x", padx=10, pady=8)
        ttk.Button(btns, text="Use selected", command=self._use).pack(side="right", padx=4)
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right", padx=4)

        body = ttk.Frame(self)
        body.pack(side="top", fill="both", expand=True, padx=10, pady=6)
        left = ttk.Frame(body)
        left.pack(side="left", fill="both", expand=True)
        self.lb = tk.Listbox(left, height=10, exportselection=False)
        for desc, g in self.cands:
            self.lb.insert("end", f"{desc}   ->   {len(g)} groups")
        self.lb.pack(side="left", fill="both", expand=True)
        lsb = ttk.Scrollbar(left, command=self.lb.yview)
        lsb.pack(side="left", fill="y")
        self.lb.configure(yscrollcommand=lsb.set)
        self.lb.bind("<<ListboxSelect>>", self._on_select)

        prev = ttk.LabelFrame(body, text="Resulting groups")
        prev.pack(side="left", fill="both", expand=True, padx=(8, 0))
        self.preview = tk.Text(prev, width=28, wrap="word", state="disabled")
        self.preview.pack(fill="both", expand=True)

        if self.cands:
            self.lb.selection_set(0)
            self._on_select()
        else:
            self._set_preview("No varying fields found in the sample names.")
        self.grab_set()

    def _set_preview(self, text):
        self.preview.configure(state="normal")
        self.preview.delete("1.0", "end")
        self.preview.insert("end", text)
        self.preview.configure(state="disabled")

    def _current(self):
        sel = self.lb.curselection()
        return self.cands[sel[0]] if sel else None

    def _on_select(self, event=None):
        c = self._current()
        self._set_preview("\n".join(c[1]) if c else "")

    def _use(self):
        c = self._current()
        if c:
            self.apply_cb(", ".join(c[1]))
            self.destroy()


class ScrollableFrame(ttk.Frame):
    """A vertically scrollable frame; put widgets in `.inner`."""
    def __init__(self, parent, width=460, height=480, **kw):
        super().__init__(parent, **kw)
        # A bounded canvas height keeps tall tab content from forcing the whole
        # window taller (which would push the Log box off-screen) -- it scrolls.
        self.canvas = tk.Canvas(self, borderwidth=0, highlightthickness=0, width=width, height=height)
        vsb = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.canvas.configure(yscrollcommand=vsb.set)
        vsb.pack(side="right", fill="y")
        self.canvas.pack(side="left", fill="both", expand=True)
        self.inner = ttk.Frame(self.canvas)
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        self.canvas.bind_all("<MouseWheel>", self._on_wheel)

    def _on_wheel(self, event):
        self.canvas.yview_scroll(int(-event.delta / 120), "units")


# --------------------------------------------------------------------------- #
# Main application
# --------------------------------------------------------------------------- #
class VolcanoGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("DIA-NN Proteomics  -  Volcano / Bubble Explorer")
        self.geometry("1360x900")

        self.result = None      # AnalysisResult after a run
        self.cfg = None
        self.cfg_vars = {}
        self.pca_vars = {}
        self.vol_vars = {}
        self.bub_vars = {}
        self._sar_text = None
        self._plot_canvas = None
        self.output_dir = None   # dedicated outputs folder under the data folder
        self._log_fh = None      # open file handle: <output_dir>/analysis_log.txt

        self.log_q = queue.Queue()
        self.result_q = queue.Queue()

        self._build_ui()
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

        nb = ttk.Notebook(main)
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

        right = ttk.LabelFrame(main, text="Plot / table")
        main.add(right, weight=1)
        self.plot_container = ttk.Frame(right)
        self.plot_container.pack(side="top", fill="both", expand=True)
        ttk.Label(self.plot_container,
                  text="Run an analysis (tab 1), then plot from tabs 2-4 or look up raw data on tab 5.",
                  foreground="#888").pack(expand=True)

    def _scroll_inner(self, tab):
        sf = ScrollableFrame(tab, width=470)
        sf.pack(side="top", fill="both", expand=True)
        return sf.inner

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
        lbl_groups = ttk.Label(grid, text="Groups (comma)")
        lbl_groups.grid(row=3, column=0, sticky="w", padx=4, pady=2)
        v["group_names"] = tk.StringVar(value="DMSO, Positive_Control, IRAK1")
        ent_groups = ttk.Entry(grid, textvariable=v["group_names"], width=34)
        ent_groups.grid(row=3, column=1, sticky="w", padx=4, pady=2)
        btn_auto = ttk.Button(grid, text="Auto-pick...", command=self._autopick_groups)
        btn_auto.grid(row=3, column=2, sticky="w", padx=4)
        _gh = ("Comma-separated group names. Each is matched as a regex against the sample-column "
               "headers, so it just needs to be a unique substring (e.g. DMSO, IRAK1).")
        Tooltip(lbl_groups, _gh)
        Tooltip(ent_groups, _gh)
        Tooltip(btn_auto, "Detect candidate group names from the sample-column headers and pick one.")
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

    def _autopick_groups(self):
        """Open a popup to pick group names from the pg matrix header columns."""
        pg = self.cfg_vars["pg_path"].get().strip()
        if not pg or not os.path.exists(pg):
            messagebox.showerror("Auto-pick groups", "Select a valid pg_matrix.tsv file first.")
            return
        try:
            import pandas as pd
            cols = list(pd.read_csv(pg, sep="\t", index_col=0, nrows=0).columns)
        except Exception:
            self.log_q.put(traceback.format_exc() + "\n")
            messagebox.showerror("Auto-pick groups", "Could not read the file header. See the Log.")
            return
        meta = {"Protein.Group", "Protein.Ids", "Protein.Names", "Genes", "First.Protein.Description"}
        samples = [c for c in cols if c not in meta]
        if len(samples) < 2:
            messagebox.showinfo("Auto-pick groups", "Not enough sample columns to analyze.")
            return
        GroupPickerDialog(self, samples, lambda s: self.cfg_vars["group_names"].set(s))

    def _preview_groups(self):
        """Compute group assignments from the pg file's columns, without running."""
        pg = self.cfg_vars["pg_path"].get().strip()
        if not pg or not os.path.exists(pg):
            messagebox.showerror("Preview groups", "Select a valid pg_matrix.tsv file first.")
            return
        groups = [g.strip() for g in self.cfg_vars["group_names"].get().split(",") if g.strip()]
        if not groups:
            messagebox.showerror("Preview groups", "Enter at least one group name.")
            return
        try:
            import pandas as pd
            from scripts.io import assign_groups
            df0 = pd.read_csv(pg, sep="\t", index_col=0, nrows=0)  # header only
            self._set_group_text(assign_groups(df0, groups))
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
        v["fdr_alpha"] = labeled_entry(emp, 1, "FDR alpha", "0.05",
                                       hint="Target empirical FDR for the fitted curve (e.g. 0.05).")
        v["kappa"] = labeled_entry(emp, 2, "kappa", "1e-6",
                                   hint="Small pseudocount used when estimating the empirical FDR ratio.")
        v["p_value_cutoff"] = labeled_entry(emp, 3, "p_value_cutoff", "1",
                                            hint="Vertical offset (in -log10 P) of the empirical-FDR curve's asymptote.")
        v["mode"] = labeled_entry(emp, 4, "Mode override", "Auto", tip="Auto / 0 / 1",
                                  hint="Side the curve is fit toward. Auto = use the analysis mode "
                                       "(0 degradation -> down side, 1 enrichment -> up side).")

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
                                       hint="Colour proteins that were imputed (special high-missing protocol) orange.")
        v["PharosTCRD"] = check(hi, 1, "Pharos TCRD classes", False,
                                hint="Colour proteins by Pharos target-development level (Tclin/Tchem/Tbio/Tdark).")
        v["highlight_kinase"] = check(hi, 2, "Protein kinases", False,
                                      hint="Mark proteins in the protein-kinase reference list.")
        v["highlight_ub"] = check(hi, 3, "Ubiquitin-related", False,
                                  hint="Mark ubiquitin-related proteins (and add the significant ones to the labels).")
        v["highlight_Gloops"] = check(hi, 4, "G-loop proteins", False,
                                      hint="Mark proteins in the G-loop reference list.")
        v["highlight_RTloops"] = check(hi, 5, "RT-loop proteins", False,
                                       hint="Mark proteins in the RT-loop reference list.")
        v["highlight_genes"] = labeled_entry(hi, 6, "Highlight genes", "", width=28,
                                             tip="UniProt IDs, comma-separated",
                                             hint="Specific proteins to mark (green) and always label. "
                                                  "Comma-separated UniProt IDs, e.g. Q13546, P51617.")

        lab = ttk.LabelFrame(parent, text="Labels")
        lab.pack(fill="x", padx=4, pady=4)
        v["label_topX_mid_fc"] = labeled_entry(lab, 0, "Label top-X mid FC", "", tip="blank = off",
                                               hint="Also label the X most significant 'mid' down-regulated proteins "
                                                    "(log2FC between -1 and -0.32). Blank/0 = off.")
        v["label_most_extreme"] = labeled_entry(lab, 1, "Label most extreme (per side)", "", tip="blank = off",
                                                hint="Label only the N points farthest from the origin on each side "
                                                     "(overrides the up/down label set). Blank = off.")
        v["max_label"] = labeled_entry(lab, 2, "Max labels", "100",
                                       hint="Skip drawing labels entirely if more than this many would be shown "
                                            "(prevents an unreadable, slow plot).")
        v["file_suffix"] = labeled_entry(lab, 3, "File suffix", "",
                                         hint="Extra text appended to the saved PNG file name.")

        place = ttk.LabelFrame(parent, text="Label selection & placement")
        place.pack(fill="x", padx=4, pady=4)
        v["label_up"] = check(place, 0, "Label up-regulated genes", True,
                              hint="Add the significant up-regulated (red) genes to the text labels.")
        v["label_down"] = check(place, 1, "Label down-regulated genes", True,
                                hint="Add the significant down-regulated (blue) genes to the text labels.")
        v["label_imputed"] = check(place, 2, "Label imputed genes", False,
                                   hint="Add the imputed (orange) proteins to the text labels.")
        v["adjust_labels"] = check(place, 3, "Auto-arrange labels (adjustText)", True,
                                   hint="Reposition labels to avoid overlap (adjustText). Off = place at the point, "
                                        "no arrows -- much faster with many labels.")
        v["adjust_arrows"] = check(place, 4, "Draw arrows to labels", True,
                                   hint="Draw thin connector lines from each moved label back to its point.")
        v["adjust_force_text"] = labeled_entry(place, 5, "Repel force (text)", "1, 2",
                                               tip="number or 'x, y'",
                                               hint="How strongly labels push apart from each other. "
                                                    "A number, or 'x, y' for separate horizontal/vertical force.")
        v["adjust_force_static"] = labeled_entry(place, 6, "Repel force (points)", "1, 2",
                                                 tip="number or 'x, y'",
                                                 hint="How strongly labels are pushed away from the data points.")

        fonts = ttk.LabelFrame(parent, text="Font sizes")
        fonts.pack(fill="x", padx=4, pady=4)
        v["title_fontsize"] = labeled_entry(fonts, 0, "Title", "24", hint="Font size of the plot title.")
        v["axis_label_fontsize"] = labeled_entry(fonts, 1, "Axis labels", "20",
                                                 hint="Font size of the x/y axis labels.")
        v["tick_fontsize"] = labeled_entry(fonts, 2, "Tick labels", "16", hint="Font size of the axis tick numbers.")
        v["legend_fontsize"] = labeled_entry(fonts, 3, "Legend", "12", hint="Font size of the legend.")
        v["gene_label_fontsize"] = labeled_entry(fonts, 4, "Gene labels", "14",
                                                 hint="Font size of the per-gene text labels on the plot.")

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
        v["compound_labelsize"] = labeled_entry(fig, 6, "Compound label size", "20",
                                                hint="Font size of the treatment/compound axis labels.")
        v["protein_labelsize"] = labeled_entry(fig, 7, "Protein label size", "20",
                                               hint="Font size of the protein axis labels.")
        v["colorFCrange"] = labeled_entry(fig, 8, "Color FC range", "-4, 0",
                                          hint="log2FC range mapped to the colour scale, e.g. '-4, 0'. "
                                               "Values outside are clipped to the ends.")
        v["legend_num"] = labeled_entry(fig, 9, "Legend # entries", "auto", tip="auto or integer",
                                        hint="Number of size-legend entries (FDR). 'auto' or an integer.")

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
                                            tip="'Desc | GENE', comma-separated; blank = all",
                                            hint="Restrict the plot to these proteins, written as "
                                                 "'Description | GENE' (the row labels), comma-separated. Blank = all.")

        fonts = ttk.LabelFrame(parent, text="Font sizes")
        fonts.pack(fill="x", padx=4, pady=4)
        ttk.Label(fonts, text="(compound / protein label sizes are in the Figure box above)",
                  foreground="#888").grid(row=0, column=0, columnspan=3, sticky="w", padx=4)
        v["title_fontsize"] = labeled_entry(fonts, 1, "Title", "30", hint="Font size of the plot title.")
        v["axis_fontsize"] = labeled_entry(fonts, 2, "Distance axis", "30",
                                           hint="Font size of the dendrogram 'Distance' axis label.")
        v["colorbar_label_fontsize"] = labeled_entry(fonts, 3, "Colorbar label", "30",
                                                     hint="Font size of the 'Log2FC Value' colorbar label.")
        v["colorbar_tick_fontsize"] = labeled_entry(fonts, 4, "Colorbar ticks", "30",
                                                    hint="Font size of the colorbar tick numbers.")
        v["legend_fontsize"] = labeled_entry(fonts, 5, "Legend", "30",
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
        """Resolve a UniProt accession or gene symbol to a Protein.Group index."""
        df = self.result.df_original
        q = query.strip()
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

        # Clear the right area and add an info header + the table.
        for w in self.plot_container.winfo_children():
            w.destroy()
        self._plot_canvas = None
        ttk.Label(self.plot_container,
                  text=f"Protein.Group: {pg}    Genes: {genes}\n{desc}    "
                       f"({len(prec)} precursor{'s' if len(prec) != 1 else ''})",
                  justify="left", foreground="#222").pack(anchor="w", padx=8, pady=(6, 2))

        container = ttk.Frame(self.plot_container)
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
        return {k: var.get() for k, var in self.cfg_vars.items()}

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
            # By default, generate volcanoes for ALL comparisons right away.
            self.after(50, self._on_plot_all)
        self.after(150, self._poll_result)

    # ----- plotting (stages 2 & 3) -----
    def _embed(self, fig):
        for w in self.plot_container.winfo_children():
            w.destroy()
        canvas = FigureCanvasTkAgg(fig, master=self.plot_container)
        canvas.draw()
        NavigationToolbar2Tk(canvas, self.plot_container).update()
        canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        self._plot_canvas = canvas

    def _embed_notebook(self, items):
        """Embed several figures, one tab per (name, figure)."""
        for w in self.plot_container.winfo_children():
            w.destroy()
        nb = ttk.Notebook(self.plot_container)
        nb.pack(side="top", fill="both", expand=True)
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
            self._embed(plt.gcf())
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
            generate_pca_plot(
                self.result.df_original, self.result.group_columns,
                filename=str(self.pca_vars["filename"].get()).strip() or "PCA_plot.png",
                title=str(self.pca_vars["title"].get()),
                text=bool(self.pca_vars["text"].get()),
                title_fontsize=_req_float(self.pca_vars["title_fontsize"].get(), 20),
                axis_fontsize=_req_float(self.pca_vars["axis_fontsize"].get(), 15),
                tick_fontsize=_opt_float(self.pca_vars["tick_fontsize"].get()),
                legend_fontsize=_opt_float(self.pca_vars["legend_fontsize"].get()),
                point_label_fontsize=_req_float(self.pca_vars["point_label_fontsize"].get(), 4),
            )
            self._embed(plt.gcf())
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
            # bubble_dendro_plot reads "<stem>_analyzed.csv" from the working dir;
            # write the analysis summary there so it has data to plot.
            analyzed = self.cfg.file.split(".")[0] + "_analyzed.csv"
            self.result.summary.to_csv(analyzed)
            plt.close("all")
            bubble_dendro_plot(SAR, self.cfg, **kwargs)
            self._embed(plt.gcf())
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
