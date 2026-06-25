"""Bubble / dendrogram SAR plot (notebook cell 12).

Body verbatim. Only changes: ``config`` is passed in so the input path uses
``config.file`` (it reads ``<file>_analyzed.csv``), and the loop annotation
lists are imported at module scope so the body keeps its bare-name references.

Note: this stage expects a ``<file>_analyzed.csv`` summary on disk (the notebook
never generated it automatically). Write ``AnalysisResult.summary`` to that path
first if you want to drive it from a real run.
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import get_cmap
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans  # Use KMeans for clustering
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram, set_link_color_palette
from matplotlib.colors import rgb2hex
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.ticker import FuncFormatter
from typing import Union

from ..reference_data import g_loop_5res_noEC, g_loop_8res_noEC, RT_loop_5res


def bubble_dendro_plot(SAR, config, SAR_suffix="", figure_filename="bubble_plot.png", fig_title="", fig_width=50, fig_height=50,
                       dendro_bubble_height_ratio=[1,3], bubble_legend_width_ratio=[20,1],
                       compound_labelsize=20, protein_labelsize=20,
                       colorFCrange=[-4,0], highlight_G_loop=0, highlight_RT_loop=0, rainbow_palette=0,
                       invert_xy=False, selected_genes=[], legend_num: Union[str, int] = "auto",
                       title_fontsize=30, axis_fontsize=30, colorbar_label_fontsize=30,
                       colorbar_tick_fontsize=30, legend_fontsize=30):
    # SAR_suffix="_1uM_vs_DMSO"
    # highlight_G_loop can be 0:none, 1: 5res-list, 2: 8res-list
    # invert_xy: if True, swap the x and y axes and do not draw the dendrogram

    # Add suffix for SAR list
    for key in SAR:
        SAR[key] = list(map(lambda s: s + SAR_suffix, SAR[key]))

    suffix_to_group = {suffix: cluster for cluster, suffixes in SAR.items() for suffix in suffixes}

    # Preparing downregulated proteins for meta-analysis
    df = pd.read_csv(config.file.split(".")[0]+"_analyzed.csv", sep=",", index_col=0)

    # Bubble plot table for x=treatment, y=protein, color=FC, size=FDR
    #df = df.drop(columns=df.filter(like="10uM").columns)

    # Identify log2FC and bh_FDR columns
    log2FC_columns = df.filter(regex="^log2FC_").columns

    # Extract suffixes for pairing
    log2FC_suffixes = [col.split("_", 1)[1] for col in log2FC_columns if col.split("_", 1)[1] in suffix_to_group]
    print(log2FC_suffixes)  # e.g. ['NGT20-11_8nM_vs_DMSO', 'NGT20-11_40nM_vs_DMSO',...]

    # Create a mask for filtering rows
    mask = pd.DataFrame(
        {
            suffix: (df[f"log2FC_{suffix}"] < -1) & (df[f"bh_FDR_{suffix}"] < 0.01)
            for suffix in log2FC_suffixes
        }
    ).any(axis=1)

    # Apply the mask to filter rows
    df_protein_downreg = df[mask]
    print(df_protein_downreg.shape)

    df_protein_downreg["Description_Genes"] = (
        df_protein_downreg["First.Protein.Description"].fillna("") + " | " + df_protein_downreg["Genes"].fillna("")
    )

    # Melt the DataFrame for plotting
    df_protein_downreg = df_protein_downreg.reset_index()
    melted_log2FC = pd.melt(
        df_protein_downreg,
        id_vars=["Description_Genes", "Protein.Group"],
        value_vars=[col for col in df_protein_downreg.columns if col.startswith("log2FC")],
        var_name="log2FC_Column",
        value_name="log2FC",
    )

    melted_log2FC["Suffix"] = melted_log2FC["log2FC_Column"].str.split("_", n=1).str[1]
    melted_log2FC["Group"] = melted_log2FC["Suffix"].map(suffix_to_group)
    melted_log2FC = melted_log2FC[melted_log2FC["Suffix"].isin(suffix_to_group.keys())]

    # Map bh_FDR values and calculate -log10 scale
    melted_log2FC["bh_FDR"] = melted_log2FC.apply(
        lambda row: df_protein_downreg[f"bh_FDR_{row['Suffix']}"].iloc[row.name % len(df_protein_downreg)], axis=1
    )
    melted_log2FC["bh_FDR_log10"] = -np.log10(melted_log2FC["bh_FDR"])
    melted_log2FC["bh_FDR_log10"] = np.nan_to_num(melted_log2FC["bh_FDR_log10"], nan=0)

    # Pivot to get a matrix of Description_Genes vs. log2FC values for clustering
    pivot_df = melted_log2FC.pivot_table(
        index="Description_Genes",
        columns="Suffix",
        values="log2FC"
    )
    pivot_df = pivot_df.fillna(1)

    # If selected_genes is provided, filter both pivot_df and melted_log2FC
    if selected_genes:
        pivot_df = pivot_df.loc[pivot_df.index.isin(selected_genes)]
        melted_log2FC = melted_log2FC[melted_log2FC["Description_Genes"].isin(selected_genes)]

    # Compute pairwise distances and hierarchical clustering
    distance_threshold = 3.5  # Set your cutoff distance here
    if pivot_df.shape[0] > 1:
        linkage_matrix = linkage(pivot_df, method="ward")  # Options: "ward", "average", etc.
        dendro = dendrogram(linkage_matrix, labels=pivot_df.index, no_plot=True)
        dendrogram_order = dendro["ivl"]
        # Define clusters using the distance threshold
        cluster_labels = fcluster(linkage_matrix, t=distance_threshold, criterion="distance")
        pivot_df["Cluster"] = cluster_labels
    else:
        # With a single protein, simply use its order and assign a default cluster
        dendrogram_order = list(pivot_df.index)
        pivot_df["Cluster"] = 1

    # Define clusters using a distance threshold
    cluster_labels = fcluster(linkage_matrix, t=distance_threshold, criterion="distance")
    pivot_df["Cluster"] = cluster_labels

    # Map cluster labels back to the original DataFrame
    cluster_mapping = pivot_df["Cluster"].to_dict()
    melted_log2FC["Cluster"] = melted_log2FC["Description_Genes"].map(cluster_mapping)

    # Convert cluster and Suffix to string (for categorical color coding)
    melted_log2FC["Cluster"] = melted_log2FC["Cluster"].astype(str)
    melted_log2FC["Suffix"] = melted_log2FC["Suffix"].astype(str)


    # Reorder categories using dendrogram order for protein descriptions
    melted_log2FC["Description_Genes"] = pd.Categorical(
        melted_log2FC["Description_Genes"],
        categories=dendrogram_order,
        ordered=True
    )
    # Order Suffix based on SAR dictionary (if needed)
    ordered_suffixes = [suffix for group in SAR.values() for suffix in group]
    melted_log2FC["Suffix"] = pd.Categorical(melted_log2FC["Suffix"], categories=ordered_suffixes, ordered=True)
    melted_log2FC.sort_values(by=["Suffix", "Description_Genes"], inplace=True)

    melted_log2FC.to_csv("temp.csv")

    # Define colors for groups (from SAR) and for clusters
    unique_groups = melted_log2FC["Group"].dropna().unique()
    group_colors = {group: get_cmap("tab20")(i / len(unique_groups)) for i, group in enumerate(unique_groups)}
    melted_log2FC["Group_Color"] = melted_log2FC["Group"].map(group_colors)
    unique_clusters = melted_log2FC["Cluster"].unique()
    cluster_colors = {cluster: get_cmap("tab20")(i / len(unique_clusters)) for i, cluster in enumerate(unique_clusters)}
    melted_log2FC["Cluster_Color"] = melted_log2FC["Cluster"].map(cluster_colors)

    Description_Genes_UniprotID_dict = df_protein_downreg.set_index("Description_Genes")["Protein.Group"].to_dict()  # {"proteinA|geneA": uniprotID, ...}

    # Set up figure and GridSpec depending on invert_xy flag
    if not invert_xy:
        # Layout with dendrogram (2 rows x 2 columns)
        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = GridSpec(2, 2, figure=fig, height_ratios=dendro_bubble_height_ratio, width_ratios=bubble_legend_width_ratio)
        # Dendrogram subplot
        ax_dendro = fig.add_subplot(gs[0, 0])
        dendro = dendrogram(
            linkage_matrix,
            labels=pivot_df.index,
            leaf_rotation=90,
            distance_sort="ascending",
            color_threshold=0,
            above_threshold_color="black",
            ax=ax_dendro
        )
        ax_dendro.set_title("", fontsize=16)
        ax_dendro.set_xticks([])
        ax_dendro.set_ylabel("Distance", fontsize=axis_fontsize)
        # Bubble plot subplot below dendrogram
        ax_bubble = fig.add_subplot(gs[1, 0])
    else:
        # Layout without dendrogram (1 row x 2 columns)
        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = GridSpec(1, 2, figure=fig, width_ratios=bubble_legend_width_ratio)
        ax_bubble = fig.add_subplot(gs[0, 0])

    # Plotting the bubble scatter plot with (x, y) depending on invert_xy.
    if not invert_xy:
        scatter = ax_bubble.scatter(
            melted_log2FC["Description_Genes"],  # x-axis: protein descriptions
            melted_log2FC["Suffix"],             # y-axis: treatment suffixes
            c=melted_log2FC["log2FC"],            # color by log2FC
            s=melted_log2FC["bh_FDR_log10"] * 500,  # size proportional to -log10(FDR)
            cmap=get_cmap("Spectral").reversed() if rainbow_palette else "coolwarm",
            alpha=0.7,
            edgecolors="w",
            clip_on=False,
        )
    else:
        scatter = ax_bubble.scatter(
            melted_log2FC["Suffix"],              # x-axis: treatment suffixes
            melted_log2FC["Description_Genes"],   # y-axis: protein descriptions
            c=melted_log2FC["log2FC"],
            s=melted_log2FC["bh_FDR_log10"] * 500,
            cmap=get_cmap("Spectral").reversed() if rainbow_palette else "coolwarm",
            alpha=0.7,
            edgecolors="w",
            clip_on=False,
        )

    # Set up legend/colorbar axes (adjust GridSpec nested location based on invert_xy)
    if not invert_xy:
        gs_nested = GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[1, 1], height_ratios=[2, 3])
        ax_empty = fig.add_subplot(gs_nested[0, 0])
        ax_empty.axis('off')
        ax_cbar = fig.add_subplot(gs_nested[1, :])
    else:
        gs_nested = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[0, 1], width_ratios=[2, 3])
        ax_empty = fig.add_subplot(gs_nested[0, 0])
        ax_empty.axis('off')
        ax_cbar = fig.add_subplot(gs_nested[:, 1])

    cbar = plt.colorbar(scatter, cax=ax_cbar, shrink=0.8, location="left", ticklocation="right")
    cbar.set_label("Log2FC Value", fontsize=colorbar_label_fontsize, rotation=270, labelpad=30)
    ax_cbar.tick_params(labelsize=colorbar_tick_fontsize)
    ax_cbar.set_aspect(5)
    scatter.set_clim(colorFCrange[0], colorFCrange[1])

    # Create size legend from scatter elements
    handles, labels = scatter.legend_elements(
        prop="sizes", alpha=0.6, num=legend_num+1 if legend_num is int else legend_num, func=lambda s: s / 500,
        fmt=FuncFormatter(lambda x, pos: f"{x:.2f}")
    )
    new_labels = [f"FDR = {10**(-float(l)):.2g}" for l in labels]
    ax_empty.legend(handles=handles, labels=new_labels, title="", fontsize=legend_fontsize,
                    labelspacing=1.5, borderaxespad=0.8, handleheight=0.2)

    # Color axis tick labels appropriately.
    if not invert_xy:
        # In standard mode: x-axis labels (proteins) colored by cluster, y-axis labels (suffixes) colored by SAR group.
        for label in ax_bubble.get_xticklabels():
            description = label.get_text()
            if description in melted_log2FC["Description_Genes"].values:
                cluster_id = melted_log2FC.loc[melted_log2FC["Description_Genes"] == description, "Cluster"].iloc[0]
                label.set_color(cluster_colors[cluster_id])
            # Highlight if needed (using external g_loop_5res_noEC / g_loop_8res_noEC)
            if highlight_G_loop == 1:
                if (description in Description_Genes_UniprotID_dict) and (Description_Genes_UniprotID_dict[description] in g_loop_5res_noEC):
                    label.set_bbox(dict(facecolor="grey", edgecolor="none", alpha=0.3))
            if highlight_G_loop == 2:
                if (description in Description_Genes_UniprotID_dict) and (Description_Genes_UniprotID_dict[description] in g_loop_8res_noEC):
                    label.set_bbox(dict(facecolor="grey", edgecolor="none", alpha=0.3))
            if highlight_RT_loop == 1:
                if (description in Description_Genes_UniprotID_dict) and (Description_Genes_UniprotID_dict[description] in RT_loop_5res):
                    label.set_bbox(dict(facecolor="pink", edgecolor="none", alpha=0.3))
        for label in ax_bubble.get_yticklabels():
            group = label.get_text()
            if suffix_to_group.get(group) in group_colors:
                label.set_color(group_colors[suffix_to_group[group]])
    else:
        # In invert mode: now x-axis labels are suffixes and y-axis labels are proteins.
        for label in ax_bubble.get_xticklabels():
            group = label.get_text()
            if suffix_to_group.get(group) in group_colors:
                label.set_color(group_colors[suffix_to_group[group]])
        for label in ax_bubble.get_yticklabels():
            description = label.get_text()
            if description in melted_log2FC["Description_Genes"].values:
                cluster_id = melted_log2FC.loc[melted_log2FC["Description_Genes"] == description, "Cluster"].iloc[0]
                label.set_color(cluster_colors[cluster_id])
            if highlight_G_loop == 1:
                if (description in Description_Genes_UniprotID_dict) and (Description_Genes_UniprotID_dict[description] in g_loop_5res_noEC):
                    label.set_bbox(dict(facecolor="grey", edgecolor="none", alpha=0.3))
            if highlight_G_loop == 2:
                if (description in Description_Genes_UniprotID_dict) and (Description_Genes_UniprotID_dict[description] in g_loop_8res_noEC):
                    label.set_bbox(dict(facecolor="grey", edgecolor="none", alpha=0.3))

    # Adjust axis labels, tick parameters, and grid.
    if not invert_xy:
        ax_bubble.set_xlabel("")
        ax_bubble.tick_params(axis='x', rotation=45, labelsize=protein_labelsize)
        for label in ax_bubble.get_xticklabels():
            label.set_ha('right')
        ax_bubble.set_xlim(-0.5, len(set(melted_log2FC["Description_Genes"])) - 0.5)
        ax_bubble.set_ylabel("")
        ax_bubble.invert_yaxis()
        ax_bubble.tick_params(axis='y', labelsize=compound_labelsize, pad=20)
    else:
        ax_bubble.set_xlabel("")
        ax_bubble.tick_params(axis='x', rotation=45, labelsize=protein_labelsize)
        for label in ax_bubble.get_xticklabels():
            label.set_ha('right')
        ax_bubble.set_ylabel("")
        ax_bubble.tick_params(axis='y', labelsize=compound_labelsize, pad=20)
        # Optionally, adjust limits if necessary:
        ax_bubble.set_xlim(-0.5, len(set(melted_log2FC["Suffix"])) - 0.5)

    # Title and grid settings.
    if not invert_xy:
        ax_dendro.set_title(fig_title, fontsize=title_fontsize, pad=10)
    ax_bubble.grid(True, linestyle="--", alpha=0.5)
    for spine in ax_bubble.spines.values():
        spine.set_visible(False)

    plt.subplots_adjust(hspace=0.01, wspace=0.01)
    plt.tight_layout()
    plt.savefig(figure_filename, dpi=200)
    plt.show()
