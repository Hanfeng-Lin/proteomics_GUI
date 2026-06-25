"""PCA plot (notebook cell 4). Function body is verbatim -- it already took its
inputs (``df``, ``group_columns``) as arguments and used no globals."""

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from adjustText import adjust_text


def generate_pca_plot(df, group_columns, filename="PCA_plot.png", title="PCA of Samples", text=False):
    """
    Generates and saves a PCA plot for the given data frame and group assignments.

    Args:
        df (pd.DataFrame): DataFrame with proteins as rows and samples as columns.
        group_columns (dict): Dictionary mapping group names to lists of column names.
        filename (str): The name of the file to save the plot.
        title (str): The title of the plot.
    """
    # Select only the sample data columns for PCA
    all_sample_columns = [col for sublist in group_columns.values() for col in sublist]
    data = df[all_sample_columns].copy()

    empty_samples = [col for col in data.columns if data[col].isnull().all()]
    if empty_samples:
        for sample in empty_samples:
            # Set a warning for the user
            print(f"Warning: Sample '{sample}' contains all NaN values and will be dropped from the PCA.")
        data.drop(columns=empty_samples, inplace=True)
    # PCA requires no missing values. We fill NaN with the mean of the protein (row).
    data_imputed = data.fillna(data.mean())

    # Transpose the data so that samples are rows and proteins are columns
    data_transposed = data_imputed.T

    # Standardize the data before performing PCA
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data_transposed)

    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(data_scaled)

    # Create a DataFrame with the principal components and group information
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=data_transposed.index)

    # Map sample names to group names
    sample_to_group = {sample: group for group, samples in group_columns.items() for sample in samples}
    pca_df['group'] = pca_df.index.map(sample_to_group)

    # Create the plot
    plt.figure(figsize=(15, 10))
    sns.scatterplot(
        x='PC1',
        y='PC2',
        hue='group',
        palette="tab20",
        data=pca_df,
        s=100,
        alpha=0.8
    )
    if text:
        texts = [plt.text(row['PC1'], row['PC2'], i.split('/')[-1], fontsize=4) for i, row in pca_df.iterrows()]
        adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))


    # Add labels and title
    plt.title(title, fontsize=20)
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.2f}% variance)', fontsize=15)
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.2f}% variance)', fontsize=15)
    plt.legend(title='Groups', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    # Save and show the plot
    plt.savefig(filename, dpi=300)
    plt.show()
