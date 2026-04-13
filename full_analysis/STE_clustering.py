# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:53:05 2026

@author: Tom
"""

#%% Import modules
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from scipy.stats import f_oneway
from sklearn.metrics import silhouette_score


#%% Load data
df = pd.read_csv("fits.csv")


#%% Clustering function



def run_computational_clustering(df, param_cols, symptom_cols, n_clusters=3):
    """
    Clusters participants based on parameters and compares symptom scores across clusters.
    """
    # 1. Prepare and Scale Data
    # We only cluster based on the computational parameters
    cluster_data = df[param_cols].dropna()
    scaler = StandardScaler()
    scaled_params = scaler.fit_transform(cluster_data)
    
    
    # 2. Fit K-Means
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(scaled_params)
    
    # Add labels back to a copy of the original dataframe
    results_df = df.loc[cluster_data.index].copy()
    results_df['Cluster'] = cluster_labels
    
    # 3. Visualize the Computational "Profiles" of each cluster
    # This shows what the clusters actually represent computationally
    cluster_means = results_df.groupby('Cluster')[param_cols].mean()
    
    plt.figure(figsize=(10, 5))
    sns.heatmap(cluster_means, cmap='RdBu_r', center=0)
    plt.title("Computational Profile of Each Cluster (Z-scored)")
    plt.show()

    # 4. Compare Clusters on Symptoms
    # We "melt" the data to make it easy for seaborn to plot
    melted_df = results_df.melt(id_vars='Cluster', value_vars=symptom_cols, 
                                var_name='Symptom', value_name='Score')
    
    plt.figure(figsize=(12, 6))
    sns.boxplot(x='Symptom', y='Score', hue='Cluster', data=melted_df)
    plt.title("Cluster differences in transdiagnsotic symptoms")
    plt.axhline(0, color='black', linestyle='--', alpha=0.3)
    plt.show()
    
    # 5. Statistical Check (One-way ANOVA)
    for symptom in symptom_cols:
        groups = [results_df[results_df['Cluster'] == i][symptom].dropna() for i in range(n_clusters)]
        f_stat, p_val = f_oneway(*groups)
        print(f"ANOVA for {symptom}: F={f_stat:.3f}, p={p_val:.3f}")
        
    return results_df






    


#%% Delta parameters
print('--- Delta parameters ---')
param_list = ['Delta_om2', 
              'Delta_logal', 
              'Delta_logzeta0', 
              'Delta_zeta1', 
              'Delta_beta0', 
              'Delta_beta1',
              'Delta_beta2',
              'Delta_beta3',
              'Delta_beta4',
              'Delta_logsa']
param_list = ['Delta_om2', 
              'Delta_logal', 
              'Delta_zeta1']
symptom_list = ['AD', 'Compul', 'SW']
clustered_df_delta = run_computational_clustering(df, param_list, symptom_list, n_clusters=3)


#%% Safe parameters
print('--- Safe parameters ---')
param_list = ['Safe_om2', 
              'Safe_logal', 
              'Safe_logzeta0', 
              'Safe_zeta1', 
              'Safe_beta0', 
              'Safe_beta1',
              'Safe_beta2',
              'Safe_beta3',
              'Safe_beta4',
              'Safe_logsa']

#param_list = ['Safe_om2', 
#              'Safe_logal', 
#              'Safe_zeta1']
symptom_list = ['AD', 'Compul', 'SW']
clustered_df_safe = run_computational_clustering(df, param_list, symptom_list, n_clusters=3)


#%% Threat parameters
print('--- Threat parameters ---')
param_list = ['Threat_om2', 
              'Threat_logal', 
              'Threat_logzeta0', 
              'Threat_zeta1', 
              'Threat_beta0', 
              'Threat_beta1',
              'Threat_beta2',
              'Threat_beta3',
              'Threat_beta4',
              'Threat_logsa']
#param_list = ['Threat_om2', 
#              'Threat_logal', 
#              'Threat_zeta1']
symptom_list = ['AD', 'Compul', 'SW']
clustered_df_threat = run_computational_clustering(df, param_list, symptom_list, n_clusters=3)


#%% Find people that transition

import plotly.graph_objects as go
import pandas as pd

import plotly.io as pio
pio.renderers.default = "browser"

def plot_cluster_transitions(df_safe, df_threat):
    # Combine labels into one dataframe
    transfer_df = pd.DataFrame({
        'Safe_Cluster': df_safe['Cluster'],
        'Threat_Cluster': df_threat['Cluster']
    })

    # Create a contingency table (who moved where)
    counts = transfer_df.groupby(['Safe_Cluster', 'Threat_Cluster']).size().reset_index(name='count')

    # Create the Sankey Diagram
    fig = go.Figure(data=[go.Sankey(
        node = dict(
          pad = 15, thickness = 20,
          line = dict(color = "black", width = 0.5),
          label = ["Safe 0", "Safe 1", "Safe 2", "Threat 0", "Threat 1", "Threat 2"],
          color = "blue"
        ),
        link = dict(
          source = counts['Safe_Cluster'], # 0, 1, 2
          target = counts['Threat_Cluster'] + 3, # 3, 4, 5
          value = counts['count']
      ))])

    fig.update_layout(title_text="Participant Transitions: Safe to Threat", font_size=10)
    fig.show()

#plot_cluster_transitions(clustered_df_safe, clustered_df_threat)



from scipy.stats import ttest_ind

def run_switching_analysis(df_safe, df_threat, symptom_cols):
    # 1. Create the Switching variable
    # We assume both DFs have the same index (Participant ID)
    switching_df = pd.DataFrame(index=df_safe.index)
    switching_df['Safe_Cluster'] = df_safe['Cluster']
    switching_df['Threat_Cluster'] = df_threat['Cluster']
    
    # Switcher = 1 if clusters are different, 0 if they are the same
    switching_df['Is_Switcher'] = (switching_df['Safe_Cluster'] != switching_df['Threat_Cluster']).astype(int)
    
    # Add symptoms back in
    for symp in symptom_cols:
        switching_df[symp] = df_safe[symp]
        
    # 2. Statistical Analysis: T-test for each symptom
    print("--- T-test: Switchers vs. Non-Switchers ---")
    for symp in symptom_cols:
        group_0 = switching_df[switching_df['Is_Switcher'] == 0][symp].dropna()
        group_1 = switching_df[switching_df['Is_Switcher'] == 1][symp].dropna()
        
        t_stat, p_val = ttest_ind(group_1, group_0) # Group 1 (Switchers) vs Group 0 (Stables)
        print(f"{symp}: t = {t_stat:.3f}, p = {p_val:.3f}")
        print(f"   Mean (Stables): {group_0.mean():.2f}")
        print(f"   Mean (Switchers): {group_1.mean():.2f}\n")

    # 3. Visualization
    melted_switch = switching_df.melt(id_vars='Is_Switcher', value_vars=symptom_cols)
    
    plt.figure(figsize=(10, 6))
    sns.barplot(x='variable', y='value', hue='Is_Switcher', data=melted_switch, palette='viridis')
    plt.title("Symptom Profiles: Computational 'Stables' vs. 'Switchers'")
    plt.ylabel("Z-scored Symptom Intensity")
    plt.xlabel("Symptom Dimension")
    plt.legend(title='Switched Cluster?', labels=['Stable', 'Switcher'])
    plt.show()

    return switching_df

#switch_results = run_switching_analysis(clustered_df_safe, clustered_df_threat, ['AD', 'Compul', 'SW'])