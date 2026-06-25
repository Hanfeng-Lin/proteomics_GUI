"""GSEA / Enrichr enrichment plot (notebook cell 16).

Body verbatim. ``df`` and ``config`` are passed in instead of read from globals
(``imputation_option`` is accepted for signature parity with the notebook but,
as in the original, is not used in the body).
"""

import gseapy
from gseapy import barplot, dotplot
import matplotlib.pyplot as plt


def gsea_analysis(treatment_group, control_group, df, config, logFC_cutoff=1, FDR_cutoff=0.01, up_regulation=False, down_regulation=False, organism="human", gsea_cutoff=0.05,
                 imputation_option=None):
    logFC="log2FC_"+treatment_group+"_vs_"+control_group
    FDR="bh_FDR_"+treatment_group+"_vs_"+control_group

    # Define your protein list with gene symbols
    if logFC_cutoff:
        if up_regulation and not down_regulation:
            updown="up"
            gene_list = df[(df[logFC] >= logFC_cutoff) & (df[FDR] <= FDR_cutoff)]["Genes"].dropna().tolist()
        if down_regulation and not up_regulation:
            updown="down"
            gene_list = df[(df[logFC]<=-logFC_cutoff)&(df[FDR]<=FDR_cutoff)]["Genes"].dropna().tolist()
        if up_regulation and down_regulation:
            updown="updown"
            gene_list = df[((df[logFC] >= logFC_cutoff) | (df[logFC] <= -logFC_cutoff)) & (df[FDR] <= FDR_cutoff)]["Genes"].dropna().tolist()


    print(f"{treatment_group}_vs_{control_group}: {updown}: {gene_list}")
    # Perform PFAM domain enrichment analysis using the 'gseapy' package
    pfam_results = gseapy.enrichr(gene_list=gene_list,
                                  organism=organism,
                                  cutoff=gsea_cutoff,
                                  gene_sets=['InterPro_Domains_2019',  'KEGG_2021_Human', 'GO_Molecular_Function_2021'],
                                  outdir=f"gsea_{treatment_group}_vs_{control_group}")


    # ax = dotplot(pfam_results.res2d, title='Pfam_Human',cmap='viridis_r', size=10, figsize=(3, 5), ofname="pfam.png")

    #print(pfam_results.res2d)

    ax = dotplot(pfam_results.results,
                 column="Adjusted P-value",
                 x='Gene_set',  # set x axis, so you could do a multi-sample/library comparsion
                 size=10,
                 top_term=5,
                 figsize=(8, 15),
                 title="Enrichment",
                 xticklabels_rot=45,  # rotate xtick labels
                 show_ring=False,  # set to False to revmove outer ring
                 marker='o',
                 ofname=f"./gsea_{updown}_{treatment_group}_vs_{control_group}.png"
                 )
    plt.show()
