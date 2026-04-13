# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 14:36:36 2026

@author: Tom
"""

#%% Import modules
import pandas as pd
import numpy as np
from sklearn.cross_decomposition import CCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

#%% Load data
df = pd.read_csv("fits.csv")


#%% Write function
def run_phenotype_analysis(df, symptom, param_delta_cols):
    """
    Runs CCA to find the relationship between symptom profiles 
    and computational parameter changes.
    """
    # 1. Prepare Data (Handle NaNs and Scale)
    # CCA is sensitive to scale, so we use Z-scores
    analysis_df = df[symptom + param_delta_cols].dropna()
    
    X = StandardScaler().fit_transform(analysis_df[symptom])
    Y = StandardScaler().fit_transform(analysis_df[param_delta_cols])
    
    print(f"Running analysis on N={len(analysis_df)} participants.")

    # 2. Fit CCA
    # n_components=1 finds the strongest single "mode" of relationship
    cca = CCA(n_components=1)
    cca.fit(X, Y)
    
    # 3. Transform data to the new "Canonical Space"
    X_c, Y_c = cca.transform(X, Y)
    
    # Calculate the correlation between the two modes
    canonical_corr = np.corrcoef(X_c[:, 0], Y_c[:, 0])[0, 1]
    
    # 4. Extract Loadings (Which symptoms and params drive the effect?)
    symptom_loadings = pd.Series(cca.x_loadings_[:, 0], index=symptom)
    param_loadings = pd.Series(cca.y_loadings_[:, 0], index=param_delta_cols)

    # --- Visualization ---
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # Plot A: The Canonical Correlation (The Phenotype)
    sns.regplot(x=X_c[:, 0], y=Y_c[:, 0], ax=axes[0], scatter_kws={'alpha':0.5})
    axes[0].set_title(f"Canonical Correlation: r = {canonical_corr:.2f}")
    axes[0].set_xlabel("Symptom Profile (Latent)")
    axes[0].set_ylabel("Computational Profile (Latent)")

    # Plot B: Loadings
    loadings = pd.concat([symptom_loadings, param_loadings])
    loadings.plot(kind='barh', ax=axes[1], color='skyblue')
    axes[1].set_title("Feature Loadings (The Phenotype Signature)")
    
    plt.tight_layout()
    plt.show()

    return symptom_loadings, param_loadings, canonical_corr



#%%

symptom_list = ['AD', 'Compul', 'SW']
delta_params = ['Delta_om2', 'Delta_logal', 'Delta_zeta1']
s_load, p_load, r = run_phenotype_analysis(df, symptom_list, delta_params)




