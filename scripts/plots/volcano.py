"""Volcano plot (notebook cell 10).

The body is copied verbatim from the notebook. Changes are confined to plumbing
that used to rely on globals:

* ``df``, ``group_columns``, ``imputation_dict`` and ``config`` are passed in.
* ``mode`` / ``imputation_option`` / ``PharosTCRD`` keep their override-able
  keyword form but default to the corresponding ``config`` value (the notebook
  defaulted them from the global settings); ``output_adjPval`` is read from
  ``config``.
* The annotation lists (Tclin/Tchem/..., kinase/ubiquitin/loop lists) are
  imported at module scope so the body can reference them by their bare names,
  exactly as before.
* The imputation-highlight block (which read the loop variable ``pair`` and the
  ``comparison_matrix``/``reference_group`` globals) is rewritten to its single
  equivalent form: both original branches resolve to
  ``imputation_dict[treatment_group + "_vs_" + control_group]`` with identical
  scatter calls, so the result is numerically identical.
"""

import logging

import numpy as np
import pandas as pd
import matplotlib.pylab as plt
import seaborn as sns
from adjustText import adjust_text

from ..stats import empirical_fdr_curve
from ..reference_data import (
    Tclin, Tchem, Tbio, Tdark,
    protein_kinase_list, protein_ubiquitination_list,
    g_loop_5res_noEC, RT_loop_5res,
)


def volcano_plot(treatment_group, control_group, *, df, group_columns, imputation_dict, config,
                 logFC_cutoff=None, logFC_cutoff2=None, FDR_cutoff=0.05, use_empirical_fdr=False, mode=None, fdr_alpha=0.05, kappa=1e-6, p_value_cutoff=1,
                 file_suffix="", highlight_genes=[], protein_level_cutoff=None,
                 xlim=[], ylim=[], x_interval=2, y_interval=1, top_buffer=0.1,
                 imputation_option=None, PharosTCRD=None,
                 highlight_kinase=False, highlight_ub=False, highlight_Gloops=False, highlight_RTloops=False, label_topX_mid_fc=None, max_label=100, label_most_extreme=None,
                 title_fontsize=24, axis_label_fontsize=20, tick_fontsize=16, legend_fontsize=12, gene_label_fontsize=14,
                 label_up=True, label_down=True, label_imputed=False,
                 adjust_labels=True, adjust_force_text=(1, 2), adjust_force_static=(1, 2), adjust_arrows=True,
                 dpi=None,
                 dot_size=40, dot_alpha=0.5,
                 color_bg="grey", color_up="red", color_down="blue",
                 color_imputed="orange", color_highlight="green"):
    # Resolve settings that the notebook bound from globals.
    if mode is None:
        mode = config.mode
    if imputation_option is None:
        imputation_option = config.imputation_option
    if PharosTCRD is None:
        PharosTCRD = config.pharos_tcrd
    output_adjPval = config.output_adjpval

    logFC="log2FC_"+treatment_group+"_vs_"+control_group
    # Plot the adjusted FDR (bh_FDR_) or the raw p (Pvalue_) on the y-axis; both
    # columns are always present in the summary now.
    if output_adjPval:
        FDR="bh_FDR_"+treatment_group+"_vs_"+control_group
    else:
        FDR="Pvalue_"+treatment_group+"_vs_"+control_group
    y_max = df[FDR].apply(lambda x:-np.log10(x)).max()

    plt.figure(figsize=(12, 9))
    plt.scatter(x=df[logFC],y=df[FDR].apply(lambda x:-np.log10(x)),s=dot_size,alpha=dot_alpha, color=color_bg)

    if logFC_cutoff2:
        slight_down = df[(df[logFC]<=-logFC_cutoff2)&(df[FDR]<=FDR_cutoff)]
        plt.scatter(x=slight_down[logFC],y=slight_down[FDR].apply(lambda x:-np.log10(x)),s=3,label=">35% Down-regulated",color=color_down)
    else:
        slight_down = None

    up, down = pd.DataFrame(), pd.DataFrame()

    if logFC_cutoff and not use_empirical_fdr:
        up = df[(df[logFC]>=logFC_cutoff)&(df[FDR]<=FDR_cutoff)]
        down = df[(df[logFC]<=-logFC_cutoff)&(df[FDR]<=FDR_cutoff)]
        plt.scatter(x=up[logFC],y=up[FDR].apply(lambda x:-np.log10(x)),s=dot_size,alpha=dot_alpha, label="Up-regulated",color=color_up)
        if not (logFC_cutoff2 or protein_level_cutoff):
            plt.scatter(x=down[logFC],y=down[FDR].apply(lambda x:-np.log10(x)),s=dot_size,alpha=dot_alpha,label="Down-regulated",color=color_down)
        plt.axvline(-logFC_cutoff,color="grey",linestyle="--")
        plt.axvline(logFC_cutoff,color="grey",linestyle="--")

    if use_empirical_fdr and mode == 0:
        res, best = empirical_fdr_curve(df, alpha=fdr_alpha, kappa=kappa, p_value_cutoff=p_value_cutoff, mode="neg", x_col=logFC, p_col=FDR)
    if use_empirical_fdr and mode == 1:
        res, best = empirical_fdr_curve(df, alpha=fdr_alpha, kappa=kappa, p_value_cutoff=p_value_cutoff, mode="pos", x_col=logFC, p_col=FDR)
    if use_empirical_fdr and best is not None:
        plt.axvline(0,color="grey",linestyle=":")

        x0, c = best['x0'], best['c']
        # Positive side
        x_pos = np.linspace(x0+1e-6, xlim[1], 200)
        y_pos = c/(x_pos-x0) + p_value_cutoff
        mask_pos = y_pos <= y_max
        plt.plot(x_pos[mask_pos], y_pos[mask_pos], 'k-', lw=0.5)
        # Negative side
        x_neg = np.linspace(xlim[0], -x0-1e-6, 200)
        y_neg = c/(-x_neg-x0) + p_value_cutoff
        mask_neg = y_neg <= y_max
        plt.plot(x_neg[mask_neg], y_neg[mask_neg], 'k-', lw=0.5)
        # Annotate best parameters
        txt = (f"Empirical FDR c/(|x|-x0)+{p_value_cutoff} optimized:\n"
               f"x0={x0:.2f}, c={c:.2f}\n"
               f"FDR={best['FDR']:.3f}, npos={best['npos']}, nneg={best['nneg']}")
        if best['FDR'] > fdr_alpha:
            txt = (f"No feasible solution for empirical FDR < {fdr_alpha} found.\n\n"
                   f"Reporting lowest empirical FDR:\n")+ txt

        plt.text(0.05, 0.95, txt, ha="left", va="top", fontsize=9, transform=plt.gca().transAxes,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.7))
        up = df[(df[FDR].apply(lambda x:-np.log10(x))>c/(df[logFC]-x0)+p_value_cutoff) & (df[logFC] > x0)]
        down = df[(df[FDR].apply(lambda x:-np.log10(x))>c/(-df[logFC]-x0)+p_value_cutoff) & (df[logFC] < -x0)]
        plt.scatter(x=up[logFC],y=up[FDR].apply(lambda x:-np.log10(x)),s=dot_size,alpha=dot_alpha, label="Up-regulated",color=color_up)
        if not (logFC_cutoff2 or protein_level_cutoff):
            plt.scatter(x=down[logFC],y=down[FDR].apply(lambda x:-np.log10(x)),s=dot_size,alpha=dot_alpha,label="Down-regulated",color=color_down)

    if protein_level_cutoff:
        lowabundance_down=df[(df[logFC]<=-logFC_cutoff) & (df[FDR]<=FDR_cutoff) & (df[group_columns[control_group]].mean(axis=1)<1000)]
        plt.scatter(x=down[logFC],y=down[FDR].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color=color_down)
        plt.scatter(x=lowabundance_down[logFC],y=lowabundance_down[FDR].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated, protein level<"+str(protein_level_cutoff),color="turquoise")

    if highlight_genes:
        # Find which of the requested genes are actually in the DataFrame index
        genes_found = [gene for gene in highlight_genes if gene in df.index and not pd.isna(df.loc[gene, logFC])]
        # Find which genes were not found so you can report them
        genes_not_found = [gene for gene in highlight_genes if gene not in df.index or pd.isna(df.loc[gene, logFC])]
        # If any of the requested genes were not found, log a warning
        if genes_not_found:
            print(f"Could not find the following genes to highlight: {genes_not_found}")

        # Proceed to plot only the genes that were actually found
        if genes_found:
            highlight = df.loc[genes_found]
            plt.scatter(x=highlight[logFC], y=highlight[FDR].apply(lambda x:-np.log10(x)), color=color_highlight, zorder=5) # zorder makes them plot on top
        else:
            logging.warning("None of the specified highlight genes were found in the data.")
            highlight = None
    else:
        highlight = None

    if imputation_option:
        imputation_proteins = df.loc[imputation_dict[treatment_group+"_vs_"+control_group]]
        plt.scatter(x=imputation_proteins[logFC], y=imputation_proteins[FDR].apply(lambda x:-np.log10(x)), s=30, color=color_imputed, label="Imputation from missing value")
    else:
        imputation_proteins = None

    # Build the set of genes to label, honoring the up/down/imputed toggles.
    # df.iloc[0:0] seeds an empty frame WITH the right columns so downstream
    # label logic works even when every toggle is off.
    label_parts = [df.iloc[0:0]]
    if label_up:
        label_parts.append(up)
    if label_down:
        label_parts.append(down)
        if slight_down is not None:
            label_parts.append(slight_down)
    if highlight is not None:                       # explicit highlight_genes always labeled
        label_parts.append(highlight)
    if label_imputed and imputation_proteins is not None:
        label_parts.append(imputation_proteins.dropna(subset=[logFC, FDR]))
    concat_df = pd.concat(label_parts).drop_duplicates()

    if PharosTCRD:
        Tclin_filtered = [item for item in Tclin if item in df.index]
        Tchem_filtered = [item for item in Tchem if item in df.index]
        Tbio_filtered = [item for item in Tbio if item in df.index]
        Tdark_filtered = [item for item in Tdark if item in df.index]
        Tclin_highlight = df.loc[Tclin_filtered]
        Tchem_highlight = df.loc[Tchem_filtered]
        Tbio_highlight = df.loc[Tbio_filtered]
        Tdark_highlight = df.loc[Tdark_filtered]
        plt.scatter(x=Tclin_highlight[logFC], y=Tclin_highlight[FDR].apply(lambda x:-np.log10(x)), s=30, alpha=1, color='#17becf', label="Tclin") # Cyan for Tclin
        plt.scatter(x=Tchem_highlight[logFC], y=Tchem_highlight[FDR].apply(lambda x:-np.log10(x)), s=30, alpha=1, color='#e377c2', label="Tchem") # Pink for Tchem
        plt.scatter(x=Tbio_highlight[logFC],   y=Tbio_highlight[FDR].apply(lambda x:-np.log10(x)), s=30, alpha=1, color='#8c564b', label="Tbio") # Brown for Tbio
        plt.scatter(x=Tdark_highlight[logFC], y=Tdark_highlight[FDR].apply(lambda x:-np.log10(x)), s=30, alpha=1, color='#7f7f7f', label="Tdark") # Grey for Tdark

    if highlight_kinase:
        kinase_filtered = [item for item in protein_kinase_list if item in df.index]
        kinase_highlight = df.loc[kinase_filtered]
        plt.scatter(x=kinase_highlight[logFC], y=kinase_highlight[FDR].apply(lambda x:-np.log10(x)), s=20, alpha=0.5, color='#17becf', label="Protein kinase") # Cyan for kinase

    if highlight_ub:
        ub_filtered = [item for item in protein_ubiquitination_list if item in df.index]
        ub_highlight = df.loc[ub_filtered]
        plt.scatter(x=ub_highlight[logFC], y=ub_highlight[FDR].apply(lambda x:-np.log10(x)), s=20, alpha=0.5, color='#17becf', label="Ubiquitin-related proteins") # Cyan for ubiquitin-related proteins
        # I want ub_highlight and logFC>1
        concat_df = pd.concat([concat_df, ub_highlight[(ub_highlight[logFC] > logFC_cutoff)|(ub_highlight[FDR]<=FDR_cutoff)]]).drop_duplicates()

    if highlight_Gloops:
        Gloops_filtered = [item for item in g_loop_5res_noEC if item in df.index]
        Gloops_highlight = df.loc[Gloops_filtered]
        plt.scatter(x=Gloops_highlight[logFC], y=Gloops_highlight[FDR].apply(lambda x:-np.log10(x)), s=20, alpha=0.5, color='#17becf', label="G-loop proteins") # Cyan for G-loops

    if highlight_RTloops:
        RTloops_filtered = [item for item in RT_loop_5res if item in df.index]
        RTloops_highlight = df.loc[RTloops_filtered]
        plt.scatter(x=RTloops_highlight[logFC], y=RTloops_highlight[FDR].apply(lambda x:-np.log10(x)), s=20, alpha=0.5, color='#e377c2', label="RT-loop proteins") # Pink for RT-loops

    if label_topX_mid_fc:
        # Identify mid FC candidates (logFC between -1 and -0.32 (20%deg))
        plt.axvline(-0.32,color="grey",linestyle=":")
        mid_fc_mask = (df[logFC] >= -1) & (df[logFC] <= -0.32) & (df[FDR] <= FDR_cutoff)
        mid_fc_candidates = df[mid_fc_mask].copy()
        # Add ranking for sorting
        mid_fc_candidates['_mid_fc_rank'] = mid_fc_candidates[FDR].rank(method='min')
        # Get top10 most significant mid FC proteins not already labeled
        existing_genes = concat_df['Genes'].unique()
        top10_mid_fc = mid_fc_candidates.sort_values('_mid_fc_rank').head(label_topX_mid_fc).loc[~mid_fc_candidates['Genes'].isin(existing_genes)]
        # Merge with main DataFrame while preserving order
        concat_df = pd.concat([concat_df, top10_mid_fc]).drop_duplicates()

    if label_most_extreme:
        extreme_df = concat_df[[logFC, FDR, 'Genes']].dropna().copy()
        extreme_df['_dist'] = np.sqrt(extreme_df[logFC]**2 + (-np.log10(extreme_df[FDR]))**2)
        left_top = extreme_df[extreme_df[logFC] < 0].nlargest(label_most_extreme, '_dist')
        right_top = extreme_df[extreme_df[logFC] > 0].nlargest(label_most_extreme, '_dist')
        concat_df = pd.concat([left_top, right_top]).drop_duplicates()

    # Generate texts for existing groups
    texts = []
    if len(concat_df) < max_label:
        for i, r in concat_df.iterrows():
            texts.append(plt.text(x=r[logFC], y=-np.log10(r[FDR]), s=r['Genes'], size=gene_label_fontsize))

    plt.xlim(xlim[0],xlim[1])
    plt.xticks(ticks=np.arange(round(xlim[0]),round(xlim[1])+x_interval,x_interval),fontsize=tick_fontsize)
    if ylim:
        plt.ylim(ylim[0],ylim[1])
        plt.yticks(ticks=np.arange(round(ylim[0]),round(ylim[1])+y_interval,y_interval),fontsize=tick_fontsize)
    else:
        plt.ylim(0, y_max*(1+top_buffer))
    plt.xlabel(logFC, labelpad=10, size=axis_label_fontsize)
    if output_adjPval:
        plt.ylabel("-log(adjusted P value)", labelpad=10, size=axis_label_fontsize)
    else:
        plt.ylabel("-log(P value)", labelpad=10, size=axis_label_fontsize)
    plt.title(logFC.split("_", maxsplit=1)[1]+"\nn="+str(len(df[FDR].dropna())),size=title_fontsize)
    if logFC_cutoff2:
        plt.axvline(-logFC_cutoff2,color="grey",linestyle="--")
    plt.axhline(-np.log10(FDR_cutoff),color="grey",linestyle="--")
    plt.legend(loc="upper right", fontsize=legend_fontsize)
    logFC=logFC[:3]+"₂"+logFC[4:]

    if adjust_labels and texts:
        arrowprops = dict(arrowstyle="-", color='black', lw=0.5) if adjust_arrows else None
        adjust_text(texts, force_text=adjust_force_text, force_static=adjust_force_static, arrowprops=arrowprops)

    imputation_suffix="_no_imputation"
    if imputation_option:
        imputation_suffix="_imputation"
    rep_suffix="_"+str(len(group_columns[control_group]))+"rep"
    plt.savefig(logFC+file_suffix+imputation_suffix+rep_suffix+'.png', transparent=True, dpi=dpi)
