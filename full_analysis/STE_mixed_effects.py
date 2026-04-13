# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 11:49:51 2026

@author: Tom
"""

#%% Import modules
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns

#%% Load data
df = pd.read_csv("fits.csv")

#%% Define model fit
def fit_lmm(df, symptom, param):
    # reshape
    safe_col = f"Safe_{param}"
    threat_col = f"Threat_{param}"

    
    df_long = pd.melt(df, 
                      id_vars=["ID", symptom],
                      value_vars=[safe_col, threat_col],
                      var_name="Condition",
                      value_name="Parameter_Value")
    
    # Drop NaNs
    df_long = df_long.dropna(subset=["Parameter_Value", symptom]).reset_index(drop=True)
    
    # Clean 'condition' to just 'Safe' or 'Threat'
    df_long["Condition"] = df_long["Condition"].str.replace(f"_{param}", "")

    # Define the formula (use C() to ensure 'safe' is the reference category)
    formula = f"Parameter_Value ~ C(Condition, Treatment('Safe')) * {symptom}"
    
    # Fit the model
    # groups = 'ID' accounts for the repeated measures (random intercept per subject)
    model = smf.mixedlm(formula, df_long, groups=df_long["ID"])
    result = model.fit()
    
    # --- 3. Visualization Logic ---
    # We create a prediction dataframe for plotting lines
    s_mean = df_long[symptom].mean()
    s_std = df_long[symptom].std()
    
    # Extract coefficients manually for the "High" and "Low" lines
    # (Handling the intercept, main effects, and interaction)
    b0 = result.params['Intercept']
    b_cond = result.params["C(Condition, Treatment('Safe'))[T.Threat]"]
    b_symp = result.params[symptom]
    b_inter = result.params[f"C(Condition, Treatment('Safe'))[T.Threat]:{symptom}"]

    def get_pred(cond_is_threat, symp_val):
        return b0 + (b_cond * cond_is_threat) + (b_symp * symp_val) + (b_inter * cond_is_threat * symp_val)

    # Build the lines
    plot_rows = []
    for label, s_val in zip(['High (+1 SD)', 'Low (-1 SD)'], [s_mean + s_std, s_mean - s_std]):
        plot_rows.append({'Condition': 'Safe', 'Value': get_pred(0, s_val), 'Group': f'{symptom} {label}'})
        plot_rows.append({'Condition': 'Threat', 'Value': get_pred(1, s_val), 'Group': f'{symptom} {label}'})
    
    viz_df = pd.DataFrame(plot_rows)

    # --- 4. Generate Plot ---
    plt.figure(figsize=(7, 5))
    
    # Plot the predicted interaction lines
    sns.lineplot(data=viz_df, x='Condition', y='Value', hue='Group', marker='o', linewidth=3)
    
    
    plt.title(f'Interaction: {symptom} x Condition on {param}')
    plt.ylabel(f'{param}')
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.show()
    
    return result




#%% AD / om2
result = fit_lmm(df, "AD", "om2")
print(result.summary())

#%% AD / alpha
result = fit_lmm(df, "AD", "logal")
print(result.summary())

#%% AD / zeta1
result = fit_lmm(df, "AD", "zeta1")
print(result.summary())

#%% SW / om2
result = fit_lmm(df, "SW", "om2")
print(result.summary())

#%% SW / alpha
result = fit_lmm(df, "SW", "logal")
print(result.summary())

#%% SW / zeta1
result = fit_lmm(df, "SW", "zeta1")
print(result.summary())

#%% Compul / om2
result = fit_lmm(df, "Compul", "om2")
print(result.summary())

#%% Compul / alpha
result = fit_lmm(df, "Compul", "logal")
print(result.summary())

#%% Compul / zeta1
result = fit_lmm(df, "Compul", "zeta1")
print(result.summary())






