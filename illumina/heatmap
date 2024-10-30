import matplotlib

matplotlib.use('TkAgg')  # Cambia il backend per evitare problemi di GUI

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Verifica la directory corrente
print("Directory corrente:", os.getcwd())


# Function to create heatmaps with extended comparisons between MOCK, ZIKV, and CMV
def create_heatmap_three_panel_extended_v2(df, title, timepoint, top_n=50):
    # Select the circRNAs to be displayed (top N by fold change)
    if timepoint == "24H":
        mock = df['A1']
        cmv = df['A2']
        zikv = df['A3']
        mock_label = 'MOCK 24H'
        zikv_label = 'ZIKV 24H'
        cmv_label = 'CMV 24H'
    elif timepoint == "48H":
        mock = df['A4']
        cmv = df['A5']
        zikv = df['A6']
        mock_label = 'MOCK 48H'
        zikv_label = 'ZIKV 48H'
        cmv_label = 'CMV 48H'
    elif timepoint == "72H":
        mock = df['A7']
        cmv = df['A8']
        zikv = df['A9']
        mock_label = 'MOCK 72H'
        zikv_label = 'ZIKV 72H'
        cmv_label = 'CMV 72H'
    else:
        raise ValueError("Timepoint not recognized. Expected '24H', '48H', or '72H'.")

    # Calculate fold changes for each comparison
    fold_change_zikv = (np.maximum(abs(zikv), abs(mock)) / (np.minimum(abs(zikv), abs(mock)) + 1e-6))
    fold_change_cmv = (np.maximum(abs(cmv), abs(mock)) / (np.minimum(abs(cmv), abs(mock)) + 1e-6))

    # Create DataFrame for sorting and selection
    df_fold_change = pd.DataFrame({
        'circRNA_ID': df['circRNA_ID'],
        'MOCK': mock,
        'ZIKV': zikv,
        'CMV': cmv,
        'Fold Change ZIKV': fold_change_zikv,
        'Fold Change CMV': fold_change_cmv
    })

    # Sort by fold change (ZIKV vs MOCK) in descending order and select the top N
    top_df_zikv = df_fold_change.sort_values(by='Fold Change ZIKV', ascending=False).head(top_n)
    top_df_cmv = df_fold_change.sort_values(by='Fold Change CMV', ascending=False).head(top_n)

    # Set the index for heatmap visualization
    top_df_zikv.set_index('circRNA_ID', inplace=True)
    top_df_cmv.set_index('circRNA_ID', inplace=True)

    # Calculate expression differences for ZIKV and CMV comparisons
    expression_diff_zikv = (top_df_zikv['ZIKV'] - top_df_zikv['MOCK']).round().astype(int)
    expression_diff_cmv = (top_df_cmv['CMV'] - top_df_cmv['MOCK']).round().astype(int)

    # Get min and max values for consistent color scaling across all heatmaps based on top N values
    vmin = min(top_df_zikv[['MOCK', 'ZIKV', 'CMV']].min().min(), top_df_cmv[['MOCK', 'ZIKV', 'CMV']].min().min())
    vmax = max(top_df_zikv[['MOCK', 'ZIKV', 'CMV']].max().max(), top_df_cmv[['MOCK', 'ZIKV', 'CMV']].max().max())

    # Plot the heatmap for MOCK vs ZIKV
    fig, axes = plt.subplots(1, 3, figsize=(65, 35), gridspec_kw={'width_ratios': [3, 1, 1]})
    fig.suptitle(f"{title} - MOCK vs ZIKV", fontsize=60)

    # Plot expression levels as a heatmap
    sns.heatmap(top_df_zikv[['MOCK', 'ZIKV']].round().astype(int), annot=True, fmt="d", cmap='viridis', ax=axes[0],
                cbar=True, annot_kws={"size": 36}, cbar_kws={'label': 'Expression Level', 'shrink': 0.8}, vmin=vmin, vmax=vmax)
    axes[0].set_title("Expression Levels", fontsize=48)
    axes[0].set_ylabel("circRNA_ID", fontsize=44)
    axes[0].set_xlabel("", fontsize=44)
    axes[0].tick_params(axis='y', labelsize=50)
    axes[0].tick_params(axis='x', labelsize=50)

    # Plot fold change as a separate heatmap
    sns.heatmap(top_df_zikv[['Fold Change ZIKV']].round().astype(int), annot=True, fmt="d", cmap='Blues', ax=axes[1],
                cbar=True, annot_kws={"size": 36}, cbar_kws={'label': 'Fold Change', 'shrink': 0.8}, vmin=vmin, vmax=vmax)
    axes[1].set_title("Fold Change", fontsize=48)
    axes[1].set_yticks([])
    axes[1].tick_params(axis='x', labelsize=50)
    axes[1].tick_params(axis='y', labelsize=50)

    # Plot expression difference between ZIKV and MOCK
    sns.heatmap(expression_diff_zikv.to_frame(name='Expression Difference (ZIKV vs MOCK)'), annot=True, fmt="d",
                cmap='Reds', ax=axes[2], cbar=True, annot_kws={"size": 36}, cbar_kws={'label': 'Expression Difference', 'shrink': 0.8}, vmin=vmin, vmax=vmax)
    axes[2].set_title("Expression Difference", fontsize=48)
    axes[2].set_yticks([])
    axes[2].tick_params(axis='x', labelsize=50)
    axes[2].tick_params(axis='y', labelsize=50)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"extended_circRNAs_heatmap_{timepoint}_MOCK_vs_ZIKV_illumina_v2.png")

    plt.close()

    # Plot the heatmap for MOCK vs CMV
    fig, axes = plt.subplots(1, 3, figsize=(65, 35), gridspec_kw={'width_ratios': [3, 1, 1]})
    fig.suptitle(f"{title} - MOCK vs CMV", fontsize=60)

    # Plot expression levels as a heatmap
    sns.heatmap(top_df_cmv[['MOCK', 'CMV']].round().astype(int), annot=True, fmt="d", cmap='viridis', ax=axes[0],
                cbar=True, annot_kws={"size": 36}, cbar_kws={'label': 'Expression Level', 'shrink': 0.8}, vmin=vmin, vmax=vmax)
    axes[0].set_title("Expression Levels", fontsize=48)
    axes[0].set_ylabel("circRNA_ID", fontsize=44)
    axes[0].set_xlabel("", fontsize=44)
    axes[0].tick_params(axis='y', labelsize=50)
    axes[0].tick_params(axis='x', labelsize=50)

    # Plot fold change as a separate heatmap
    sns.heatmap(top_df_cmv[['Fold Change CMV']].round().astype(int), annot=True, fmt="d", cmap='Blues', ax=axes[1],
                cbar=True, annot_kws={"size": 36}, cbar_kws={'label': 'Fold Change', 'shrink': 0.8}, vmin=vmin, vmax=vmax)
    axes[1].set_title("Fold Change", fontsize=48)
    axes[1].set_yticks([])
    axes[1].tick_params(axis='x', labelsize=50)
    axes[1].tick_params(axis='y', labelsize=50)

    # Plot expression difference between CMV and MOCK
    sns.heatmap(expression_diff_cmv.to_frame(name='Expression Difference (CMV vs MOCK)'), annot=True, fmt="d",
                cmap='Reds', ax=axes[2], cbar=True, annot_kws={"size": 36}, cbar_kws={'label': 'Expression Difference', 'shrink': 0.8}, vmin=vmin, vmax=vmax)
    axes[2].set_title("Expression Difference", fontsize=48)
    axes[2].set_yticks([])
    axes[2].tick_params(axis='x', labelsize=50)
    axes[2].tick_params(axis='y', labelsize=50)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"extended_circRNAs_heatmap_{timepoint}_MOCK_vs_CMV_illumina_v2.png")

    plt.close()


# Load the datasets
# The CSV files should contain columns like circRNA_ID, A1, A2, A3, etc., corresponding to different experimental conditions.
file_paths = {
    "24H": "top_circRNAs_illumina_24H.csv",
    "48H": "top_circRNAs_illumina_48H.csv",
    "72H": "top_circRNAs_illumina_72H.csv"
}

# Load each CSV into a DataFrame
data_frames = {timepoint: pd.read_csv(path) for timepoint, path in file_paths.items()}

# Create and save extended heatmaps for each timepoint dataset, showing more circRNAs for MOCK vs ZIKV and MOCK vs CMV
for timepoint, df in data_frames.items():
    if timepoint in ["24H", "48H", "72H"]:  # Only create for valid timepoints
        create_heatmap_three_panel_extended_v2(df, f"Extended circRNAs Heatmap - {timepoint}", timepoint, top_n=20)
