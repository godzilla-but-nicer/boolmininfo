import pandas as pd
import numpy as np
import statsmodels.formula.api as smf

# open df get list for new df rows
dyn = pd.read_csv(snakemake.input.dyn_csv, index_col=0)
df_list = []

# one set of features per rule
for rule in range(2**8):
    new_row = {}
    new_row['rule'] = rule
    eca = dyn[dyn['rule'] == rule]

    # mean transient length
    results = smf.ols('np.log2(mean_transient) ~ n_cells', data=eca).fit()
    r2 = results.rsquared
    slope = results.params[1]
    new_row['mean_transient'] = slope

    # variance of transient length
    results = smf.ols('np.log2(var_transient) ~ n_cells', data=eca).fit()
    r2 = results.rsquared
    slope = results.params[1]
    new_row['var_transient'] = slope

    # mean attractor period
    results = smf.ols('np.log2(mean_period) ~ n_cells', data=eca).fit()
    r2 = results.rsquared
    slope = results.params[1]
    new_row['mean_period'] = slope

    # variance of attractor period
    # set zeros to the min value, this is probably a bad strategy
    eca['var_period'] = eca['var_period'].replace(
        0.0, eca[eca['var_period'] > 0]['var_period'].min())
    
    # if the new min value is not zero, do regression
    if eca['var_period'].min() > 0:
        results = smf.ols('np.log2(var_period) ~ n_cells', data=eca).fit()
        r2 = results.rsquared
        slope = results.params[1]
        new_row['var_period'] = slope
    # otherwise put in nan
    else:
        new_row['var_period'] = np.nan
    
    # number of attractors
    results = smf.ols('np.log2(m_attractors) ~ n_cells', data=eca).fit()
    r2 = results.rsquared
    slope = results.params[1]
    new_row['m_attractors'] = slope

    # derrida coefficient should be equal for all number of cells
    new_row['derrida_coeff'] = eca['derrida_coeff'].min()

    df_list.append(new_row)

out_df = pd.DataFrame(df_list)
out_df.to_csv(snakemake.output[0])
