

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data (assuming data is already loaded as 'data')
# Set epsilon for calculations to prevent division by zero
epsilon = 1e-6

# Calculate fold changes for each condition, adding epsilon only in the denominator during calculations
data['FC_CMV_24H'] = data['A2'] / (data['A1'] + epsilon)
data['FC_ZIKV_24H'] = data['A3'] / (data['A1'] + epsilon)

data['FC_CMV_48H'] = data['A5'] / (data['A4'] + epsilon)
data['FC_ZIKV_48H'] = data['A6'] / (data['A4'] + epsilon)

data['FC_CMV_72H'] = data['A8'] / (data['A7'] + epsilon)
data['FC_ZIKV_72H'] = data['A9'] / (data['A7'] + epsilon)

# Find the top circRNAs with highest fold changes in both Control vs CMV and Control vs ZIKV at each timepoint
top_zikv_24h = data[['circRNA_ID', 'FC_ZIKV_24H']].nlargest(5, 'FC_ZIKV_24H')
top_cmv_24h = data[['circRNA_ID', 'FC_CMV_24H']].nlargest(5, 'FC_CMV_24H')

top_zikv_48h = data[['circRNA_ID', 'FC_ZIKV_48H']].nlargest(5, 'FC_ZIKV_48H')
top_cmv_48h = data[['circRNA_ID', 'FC_CMV_48H']].nlargest(5, 'FC_CMV_48H')

top_zikv_72h = data[['circRNA_ID', 'FC_ZIKV_72H']].nlargest(5, 'FC_ZIKV_72H')
top_cmv_72h = data[['circRNA_ID', 'FC_CMV_72H']].nlargest(5, 'FC_CMV_72H')

# Getting the intersection of the top circRNAs for both comparisons at each timepoint
top_combined_24h = pd.merge(top_zikv_24h, top_cmv_24h, on='circRNA_ID')
top_combined_48h = pd.merge(top_zikv_48h, top_cmv_48h, on='circRNA_ID')
top_combined_72h = pd.merge(top_zikv_72h, top_cmv_72h, on='circRNA_ID')

# Extracting the IDs of the circRNAs that are in the top for both comparisons
combined_circRNA_ids = pd.concat([top_combined_24h['circRNA_ID'], top_combined_48h['circRNA_ID'], top_combined_72h['circRNA_ID']]).unique()

# Filtering the original data to include only the top circRNAs for both comparisons
filtered_combined_data = data[data['circRNA_ID'].isin(combined_circRNA_ids)].copy()
filtered_combined_data.set_index('circRNA_ID', inplace=True)

# Selecting only the relevant columns (A1 to A9) for heatmap generation without modifying original expression values
combined_heatmap_data = filtered_combined_data[['A1', 'A2', 'A3', 'A4', 'A5', 'A6', 'A7', 'A8', 'A9']]

# Extracting the fold change values for each time point
fc_combined_24h_data = filtered_combined_data[['FC_CMV_24H', 'FC_ZIKV_24H']]
fc_combined_48h_data = filtered_combined_data[['FC_CMV_48H', 'FC_ZIKV_48H']]
fc_combined_72h_data = filtered_combined_data[['FC_CMV_72H', 'FC_ZIKV_72H']]

# Convert the numerical values in the heatmaps to integers for better readability
combined_heatmap_data_int = combined_heatmap_data.round().astype(int)
fc_combined_24h_data_int = fc_combined_24h_data.round().astype(int)
fc_combined_48h_data_int = fc_combined_48h_data.round().astype(int)
fc_combined_72h_data_int = fc_combined_72h_data.round().astype(int)

# Sorting the fold change data and corresponding expression levels in descending order of fold change for better visualization
fc_combined_24h_sorted = fc_combined_24h_data_int.sort_values(by='FC_ZIKV_24H', ascending=False)
combined_heatmap_24h_sorted = combined_heatmap_data_int.loc[fc_combined_24h_sorted.index, ['A1', 'A2', 'A3']]

fc_combined_48h_sorted = fc_combined_48h_data_int.sort_values(by='FC_ZIKV_48H', ascending=False)
combined_heatmap_48h_sorted = combined_heatmap_data_int.loc[fc_combined_48h_sorted.index, ['A4', 'A5', 'A6']]

fc_combined_72h_sorted = fc_combined_72h_data_int.sort_values(by='FC_ZIKV_72H', ascending=False)
combined_heatmap_72h_sorted = combined_heatmap_data_int.loc[fc_combined_72h_sorted.index, ['A7', 'A8', 'A9']]

# Plotting heatmap for 24H time point with fold change alongside (sorted by fold change, integers)
fig, ax = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={'width_ratios': [3, 1]})

# Heatmap for expression levels at 24H with integer values, sorted
sns.heatmap(combined_heatmap_24h_sorted, annot=True, fmt="d", cmap='viridis', linewidths=.5, cbar_kws={'label': 'Expression Level'}, ax=ax[0])
ax[0].set_title('Heatmap of circRNA Expression Levels for 24H (Sorted by Fold Changes for ZIKV and CMV)')
ax[0].set_xlabel('Samples (Conditions at 24H)')
ax[0].set_ylabel('circRNA IDs')

# Heatmap for fold changes at 24H with integer values, sorted
sns.heatmap(fc_combined_24h_sorted, annot=True, fmt="d", cmap='coolwarm', linewidths=.5, cbar_kws={'label': 'Fold Change'}, ax=ax[1])
ax[1].set_title('Fold Change Heatmap (24H, Sorted)')
ax[1].set_xlabel('Fold Changes')
ax[1].set_yticks([])  # Hide y-tick labels as they're the same as the main heatmap

plt.tight_layout()
plt.show()

# Plotting heatmap for 48H time point with fold change alongside (sorted by fold change, integers)
fig, ax = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={'width_ratios': [3, 1]})

# Heatmap for expression levels at 48H with integer values, sorted
sns.heatmap(combined_heatmap_48h_sorted, annot=True, fmt="d", cmap='viridis', linewidths=.5, cbar_kws={'label': 'Expression Level'}, ax=ax[0])
ax[0].set_title('Heatmap of circRNA Expression Levels for 48H (Sorted by Fold Changes for ZIKV and CMV)')
ax[0].set_xlabel('Samples (Conditions at 48H)')
ax[0].set_ylabel('circRNA IDs')

# Heatmap for fold changes at 48H with integer values, sorted
sns.heatmap(fc_combined_48h_sorted, annot=True, fmt="d", cmap='coolwarm', linewidths=.5, cbar_kws={'label': 'Fold Change'}, ax=ax[1])
ax[1].set_title('Fold Change Heatmap (48H, Sorted)')
ax[1].set_xlabel('Fold Changes')
ax[1].set_yticks([])  # Hide y-tick labels as they're the same as the main heatmap

plt.tight_layout()
plt.show()

# Plotting heatmap for 72H time point with fold change alongside (sorted by fold change, integers)
fig, ax = plt.subplots(1, 2, figsize=(18, 8), gridspec_kw={'width_ratios': [3, 1]})

# Heatmap for expression levels at 72H with integer values, sorted
sns.heatmap(combined_heatmap_72h_sorted, annot=True, fmt="d", cmap='viridis', linewidths=.5, cbar_kws={'label': 'Expression Level'}, ax=ax[0])
ax[0].set_title('Heatmap of circRNA Expression Levels for 72H (Sorted by Fold Changes for ZIKV and CMV)')
ax[0].set_xlabel('Samples (Conditions at 72H)')
ax[0].set_ylabel('circRNA IDs')

# Heatmap for fold changes at 72H with integer values, sorted
sns.heatmap(fc_combined_72h_sorted, annot=True, fmt="d", cmap='coolwarm', linewidths=.5, cbar_kws={'label': 'Fold Change'}, ax=ax[1])
ax[1].set_title('Fold Change Heatmap (72H, Sorted)')
ax[1].set_xlabel('Fold Changes')
ax[1].set_yticks([])  # Hide y-tick labels as they're the same as the main heatmap

plt.tight_layout()
plt.show()



