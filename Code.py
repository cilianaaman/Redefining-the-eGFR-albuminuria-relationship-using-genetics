#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 12:06:50 2024

@author: cilianaaman
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:04:36 2024

@author: cilianaaman
"""
import openpyxl
import numpy as np
import pandas as pd
import statsmodels.api as sm
from sklearn import preprocessing
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text




#%% Load eGFR Data
#########################
#### Load eGFR Data #####
#########################


# Load Excel file
df = pd.read_excel('/Users/cilianaaman/Documents/Uni/Bachelor/6. sem/Bachelorprojekt/eGFR.xlsx')

df=df.drop(labels="LocusCoordinates", axis=1)

# Create a matrix

eGFR_matrix_overall = df.drop(range(424,453))

del df

#%% UKBB Data import (run one time)
#####################
####  UKBB Data #####
#####################

#########################
### Load the alb data ###
#########################

df_Genomic_alb = pd.read_excel('/Users/cilianaaman/Documents/Uni/Bachelor/6. sem/Bachelorprojekt/Data UKBB/Data_UKBB_alb.xlsx', sheet_name='Ark1')

# Select columns starting with "rs"
rs_columns = [col for col in df_Genomic_alb.columns if col.startswith('rs')]

# round the values 
df_Genomic_alb[rs_columns] = df_Genomic_alb[rs_columns].apply(np.round)


##########################
### Load the eGFR data ###
##########################


df_Genomic_eGFR = pd.read_excel('/Users/cilianaaman/Documents/Uni/Bachelor/6. sem/Bachelorprojekt/Data UKBB/Data_ukbb_eGFR.xlsx', sheet_name='Ark1')

# Select columns starting with "rs"
rs_columns = [col for col in df_Genomic_eGFR.columns if col.startswith('rs')]

# round the values 
df_Genomic_eGFR[rs_columns] = df_Genomic_eGFR[rs_columns].apply(np.round)




###############################
### Load the phenotype data ###
###############################


df_phenotype = pd.read_excel('/Users/cilianaaman/Documents/Uni/Bachelor/6. sem/Bachelorprojekt/Data UKBB/phenotype_data_alb_eGFR.xlsx', sheet_name='phenotype_data2', decimal=',')


import pickle
with open('variables.pkl', 'wb') as f:
    pickle.dump((df_phenotype, df_Genomic_eGFR, df_Genomic_alb), f)

#%% Getting the UKBB imported data fast + creating 2 matresis with the phenotype and genotype data
from openpyxl import load_workbook
import pickle


# Load variables
with open('variables.pkl', 'rb') as f:
    df_phenotype, df_Genomic_eGFR, df_Genomic_alb = pickle.load(f)



ethnicity = pd.read_excel('/Users/cilianaaman/Documents/Uni/Bachelor/6. sem/Bachelorprojekt/Data UKBB/ethnicity.xlsx', sheet_name='ethnicity', decimal=',')

################
### Withdraw ###
################


# Load the participants who have withdrawn
withdraw = pd.read_csv('/Users/cilianaaman/Desktop/withdraw.csv', header=None)

# Rename the column to 'Eid'
withdraw.columns = ['Eid']



################################################################
### creating 2 matresis with the phenotype and genotype data ###
################################################################


#####################
### the eGFR data ###
#####################

df_eGFR_data = pd.merge(df_phenotype, df_Genomic_eGFR, left_on='eid', right_on='IID')

# Specify the columns to delete
columns_to_delete = ['BMI', 'micro albumin in urin', 'Unnamed: 9', 'Unnamed: 10', 'Unnamed: 11', 'Unnamed: 12', 'Unnamed: 13','Unnamed: 0','SEX', 'IID']

# Drop the specified columns from the DataFrame
df_eGFR_data = df_eGFR_data.drop(columns_to_delete, axis=1)

# Filter out rows where 'Eid' in Overall_matrix matches 'Eid' in withdraw
df_eGFR_data = df_eGFR_data[~df_eGFR_data['eid'].isin(withdraw['Eid'])]

# Remove rows with NaN values
df_eGFR_data = df_eGFR_data.dropna()



######################
### overall matrix ###
######################

genomic_eGFR_columns = [col for col in df_Genomic_eGFR.columns if col.startswith('rs')]
rename_mapping = {col: col + '_eGFR' for col in genomic_eGFR_columns}
df_Genomic_eGFR = df_Genomic_eGFR.rename(columns=rename_mapping)


Overall_matrix = pd.merge(df_phenotype, df_Genomic_eGFR, left_on='eid', right_on='IID')
# Specify the columns to delete
columns_to_delete = ['Unnamed: 0', 'BMI', 'eGFR_Creatinine_2021_0','Unnamed: 9' ,'Unnamed: 10', 'Unnamed: 11', 'Unnamed: 12', 'Unnamed: 13']

# Drop the specified columns from the DataFrame
Overall_matrix = Overall_matrix.drop(columns_to_delete, axis=1)

Overall_matrix = pd.merge(Overall_matrix, ethnicity, left_on='eid', right_on='eid')

# Filter out rows where 'Eid' in Overall_matrix matches 'Eid' in withdraw
Overall_matrix = Overall_matrix[~Overall_matrix['eid'].isin(withdraw['Eid'])]


# Deleting variables that are not needed     
variables_to_delete = ['columns_to_delete','df_Genomic_eGFR','f','genomic_eGFR_columns','rename_mapping']


for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete

#%% Ethnicity european - How many ?

Overall_matrix['22006-0.0'] = Overall_matrix['22006-0.0'].fillna(0)

value_counts = Overall_matrix['21000-0.0'].value_counts()

total_count = value_counts.sum()

desired_values = [1001, 1002, 1003]
desired_counts = value_counts.loc[desired_values].sum()

# Calculate the overall percentage
overall_percentage = (desired_counts / total_count) * 100

# Print the overall percentage
print("Overall percentage: {:.2f}%".format(overall_percentage))

value_counts = Overall_matrix['22006-0.0'].value_counts()

# Get the count of 1 values
count_1 = value_counts.get(1, 0)

# Get the total number of values in the column
total_count = Overall_matrix['22006-0.0'].count()

# Calculate the percentage of 1s compared to the entire column
percentage_1 = (count_1 / total_count) * 100

print(f"Number of 1 values: {count_1}")
print(f"Percentage of 1 values: {percentage_1}%")

Overall_matrix = Overall_matrix.drop(['22006-0.0', '21000-0.0'], axis=1)


#%% Calculating Allele Frequency UK Biobank

allele_frequency = pd.DataFrame(columns=['rsID', 'Allele Frequency'])

for column in Overall_matrix.columns:
    if column.startswith('rs'):
        allele_sum = Overall_matrix[column].sum()
        total_count = len(Overall_matrix[column]) * 2
        allele_frequency.loc[len(allele_frequency)] = [column, allele_sum / total_count]




#%% Preprocessing data 

##########################
### Preprocessing data ###
##########################



###############################
### calculating eGFR values ###
###############################


# Convert creatinine from µmol/L to mg/dL
Overall_matrix['Creatinine'] = Overall_matrix['Creatinine'] * 0.0113123


sex_column = Overall_matrix['Sex']
creatinine_column = Overall_matrix['Creatinine']
age_column = Overall_matrix['Age at recruitment']

   
def calculate_eGFR_2021(sex, creatinine, age):
    if sex == 1:
        return 142*((min(creatinine/0.9, 1))**(-0.302))*((max(creatinine/0.9, 1))**(-1.200))*(0.9938**age)
    elif sex == 0:
        return 142*(min(creatinine/0.7, 1))**(-0.241)*(max(creatinine/0.7, 1))**(-1.200)*0.9938**age*1.012

egfr_2021 = []

for eid, sex, creatinine, age in zip(Overall_matrix['eid'], sex_column, creatinine_column, age_column):
    if not np.isnan(sex) and not np.isnan(creatinine) and not np.isnan(age):
        eGFR_2021 = calculate_eGFR_2021(sex, creatinine, age)
    else:
        eGFR_2021 = np.nan
    egfr_2021.append((eid, eGFR_2021))

egfr_2021 = pd.DataFrame(egfr_2021, columns=['eid', 'eGFR_2021'])

Overall_matrix = Overall_matrix.merge(egfr_2021, on='eid', how='left')


########################################
### Overall matrix removing outliers ###
########################################


eid_column = Overall_matrix['eid']

# Select the columns with phenotype data
phenotype_columns = Overall_matrix.iloc[:, 1:6]

# Calculate the mean and standard deviation for each phenotype column excluding NaN values
mean_values = np.nanmean(phenotype_columns, axis=0)
std_values = np.nanstd(phenotype_columns, axis=0)

# Define the lower and upper bounds for outlier removal
lower_bounds = mean_values - 4 * std_values
upper_bounds = mean_values + 4 * std_values

# Identify the rows where any phenotype value is outside the bounds, excluding NaN values
outlier_rows = np.any((phenotype_columns < lower_bounds) | (phenotype_columns > upper_bounds), axis=1) & np.all(~np.isnan(phenotype_columns), axis=1)

# Remove the outlier rows from the Overall_matrix matrix
Overall_matrix_filtered = Overall_matrix[~outlier_rows].copy()

Overall_matrix_filtered.loc[:, 'eid'] = eid_column.loc[~outlier_rows].values

# Convert UACR to ln(UACR) and eGFR to ln(eGFR)
Overall_matrix_filtered['micro albumin in urin'] = np.log(Overall_matrix_filtered['micro albumin in urin'])

Overall_matrix_filtered['eGFR_2021'] = np.log(Overall_matrix_filtered['eGFR_2021'])


#############################
### creating eGFR matrix  ###
#############################

columns_to_keep = ['eid', 'Sex', 'Systolic blood pressure', 'Age at recruitment','eGFR_2021']

# Columns to filter and rename
columns_to_filter = [col for col in Overall_matrix_filtered.columns if col.endswith('_eGFR')]
filtered_columns = [col.replace('_eGFR', '') for col in columns_to_filter]

# Create the new DataFrame with filtered columns
df_eGFR_data_filtered = Overall_matrix_filtered[columns_to_keep + columns_to_filter].copy()
df_eGFR_data_filtered.columns = columns_to_keep + filtered_columns


# Drop rows with NaN values
df_eGFR_data_filtered = df_eGFR_data_filtered.dropna()






# Deleting variables that are not needed     
variables_to_delete = ['eid_column','lower_bounds','mean_values','outlier_rows','std_values',
                       'upper_bounds']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete




#%% preprocessing scale


##########################
### Preprocessing data ###
##########################



columns_to_scale = ['Sex', 'Systolic blood pressure', 'Age at recruitment']

# Select all columns that start with 'rs'
rs_columns = [col for col in Overall_matrix_filtered.columns if col.startswith('rs')]

# Combine the columns to scale and the rs columns
columns_to_scale.extend(rs_columns)

df_selected_columns = Overall_matrix_filtered[columns_to_scale]

# Scale the selected columns
scaled_columns = preprocessing.scale(df_selected_columns)

Overall_matrix_stand = pd.DataFrame(index=Overall_matrix_filtered.index, data=scaled_columns, columns=columns_to_scale)

# Include the columns that were not scaled
non_scaled_columns = [col for col in Overall_matrix_filtered.columns if col not in columns_to_scale]
Overall_matrix_stand[non_scaled_columns] = Overall_matrix_filtered[non_scaled_columns]


columns_to_drop = ['IID', 'SEX']
Overall_matrix_stand.drop(columns_to_drop, axis=1, inplace=True)



#%% MLR 1 2021 eGFR 



###########################
### MLR 1 for eGFR 2021 ###
###########################


# Select the columns with '_eGFR' at the end and include 'eid', 'Sex', and 'age' and SBP
columns_to_select = [col for col in Overall_matrix_stand.columns if col.endswith('_eGFR')] + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure','eGFR_2021']

# Create a new DataFrame with the selected columns
new_matrix = Overall_matrix_stand[columns_to_select].copy()

# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)


y_log = new_matrix['eGFR_2021']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex

# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
#bonferroni_threshold = 0.05 / len(rs_columns)  # Bonferroni corrected threshold

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    p_value = result.pvalues[rsID]
    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
MLR_eGFR_2021_all = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
#MLR_eGFR_2021_significant_results_bonferroni_all = MLR_eGFR_2021_all[MLR_eGFR_2021_all['P-value'] <= bonferroni_threshold]

# Filter non-significant rsIDs based on a significance level of 0.05
MLR_eGFR_2021_significant_results_all = MLR_eGFR_2021_all[MLR_eGFR_2021_all['P-value'] <= 0.05]



variables_to_delete = ['model','p_value','result','results','rs_columns',
                       'rsID','X','X_with_rsID','y_log']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete



#%% MLR 2


###########################################
### eGFR 2021 rsID associated with UACR ### 
###########################################

rsIDs_to_select = MLR_eGFR_2021_significant_results_all['rsID'].tolist()

new_matrix = Overall_matrix_stand[rsIDs_to_select + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure','micro albumin in urin']].copy()

# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)

y_log = new_matrix['micro albumin in urin']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex


# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
bonferroni_threshold_eGFR_2021_associated_with_UACR = 0.05 / len(rs_columns)  # Bonferroni corrected threshold

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    p_value = result.pvalues[rsID]
    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
eGFR_2021_rsid_associated_with_UACR_model2 = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value','CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
eGFR_2021_rsid_associated_with_UACR_bonferroni_model2 = eGFR_2021_rsid_associated_with_UACR_model2[eGFR_2021_rsid_associated_with_UACR_model2['P-value'] <= bonferroni_threshold_eGFR_2021_associated_with_UACR]

# Filter non-significant rsIDs based on a significance level of 0.05
eGFR_2021_rsid_associated_with_UACR_significant_model2 = eGFR_2021_rsid_associated_with_UACR_model2[eGFR_2021_rsid_associated_with_UACR_model2['P-value'] <= 0.05]


#%% Preprocessing - creating matrix for scatterplots


##############################################################
### eGFR_2021_rsid_associated_with_UACR_significant_model2 ###
##############################################################

rsID = eGFR_2021_rsid_associated_with_UACR_model2['rsID'].str.split('_', expand=True)[0]
EA = eGFR_2021_rsid_associated_with_UACR_model2['rsID'].str.split('_', expand=True)[1]
chr_values = []
gene_values = []


for rs in rsID:
    row = eGFR_matrix_overall[eGFR_matrix_overall['rsID'].str.startswith(rs)]
    chr_values.append(row['chr'].values[0])
    gene_values.append(row['Gene'].values[0])

data = {
    'rsID_eGFR_2021': rsID,
    'chr': chr_values,
    'Gene': gene_values,
    'EA': EA,
    'Beta_UACR': eGFR_2021_rsid_associated_with_UACR_model2['Effect Size'],
    'p_val_UACR': eGFR_2021_rsid_associated_with_UACR_model2['P-value'],
    'UACR CI Lower': eGFR_2021_rsid_associated_with_UACR_model2['CI Lower'],
    'UACR CI Upper': eGFR_2021_rsid_associated_with_UACR_model2['CI Upper']
}

MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2 = pd.DataFrame(data)



#########################
### UACR X eGFR 2021  ### eGFR
#########################


Scatter_plot_eGFR_2021_associated_with_UACR = pd.DataFrame()

Scatter_plot_eGFR_2021_associated_with_UACR['rsID_eGFR_2021'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['rsID_eGFR_2021']
Scatter_plot_eGFR_2021_associated_with_UACR['chr'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['chr']
Scatter_plot_eGFR_2021_associated_with_UACR['Gene'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['Gene']
Scatter_plot_eGFR_2021_associated_with_UACR['EA'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['EA']
Scatter_plot_eGFR_2021_associated_with_UACR['Beta_UACR'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['Beta_UACR']
Scatter_plot_eGFR_2021_associated_with_UACR['p_val_UACR'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['p_val_UACR']
Scatter_plot_eGFR_2021_associated_with_UACR['UACR CI Lower'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['UACR CI Lower']
Scatter_plot_eGFR_2021_associated_with_UACR['UACR CI Upper'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['UACR CI Upper']

Scatter_plot_eGFR_2021_associated_with_UACR['Beta_eGFR_2021'] = np.nan
Scatter_plot_eGFR_2021_associated_with_UACR['p_val_eGFR_2021'] = np.nan
Scatter_plot_eGFR_2021_associated_with_UACR['eGFR CI Lower'] = np.nan
Scatter_plot_eGFR_2021_associated_with_UACR['eGFR CI Upper'] = np.nan


for idx, row in MLR_eGFR_2021_all.iterrows():
    rsID_prefix = row['rsID'].split('_')[0]
    if rsID_prefix in MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2['rsID_eGFR_2021'].values:
        matching_row = Scatter_plot_eGFR_2021_associated_with_UACR[Scatter_plot_eGFR_2021_associated_with_UACR['rsID_eGFR_2021'] == rsID_prefix]
        Scatter_plot_eGFR_2021_associated_with_UACR.loc[matching_row.index, 'Beta_eGFR_2021'] = row['Effect Size']
        Scatter_plot_eGFR_2021_associated_with_UACR.loc[matching_row.index, 'p_val_eGFR_2021'] = row['P-value']
        Scatter_plot_eGFR_2021_associated_with_UACR.loc[matching_row.index, 'eGFR CI Lower'] = row['CI Lower']
        Scatter_plot_eGFR_2021_associated_with_UACR.loc[matching_row.index, 'eGFR CI Upper'] = row['CI Upper']






#%% Scatter plot 2021 


############################################
### filtereing only the significant data ###
############################################



Scatter_plot_eGFR_2021_associated_with_UACR_significant = Scatter_plot_eGFR_2021_associated_with_UACR[
    (Scatter_plot_eGFR_2021_associated_with_UACR['p_val_eGFR_2021'] <=  0.05) &
    (Scatter_plot_eGFR_2021_associated_with_UACR['p_val_UACR'] <=  0.05)
].reset_index(drop=True)



############################################################
####  Scatterplot eGFR 2021 rsID association with UACR ##### all without names 
############################################################


Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'], errors='coerce')
Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'], errors='coerce')

Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = 'A'
# Extract base rsIDs from the Bonferroni dataset
uacr_bonferroni_rsIDs = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2['rsID'].str.split('_').str[0]

# Check if rsID_eGFR is in the UACR Bonferroni list and assign 'B' to plotting_groups
Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = np.where(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant['rsID_eGFR_2021'].isin(uacr_bonferroni_rsIDs),
    'B',
    Scatter_plot_eGFR_2021_associated_with_UACR_significant.get('plotting_groups', np.nan)  # Preserve existing values if any
)



# Create the scatter plot 
g = sns.relplot(data=Scatter_plot_eGFR_2021_associated_with_UACR_significant,
                x='Beta_eGFR_2021',
                y='Beta_UACR',
                aspect=2,
                hue='plotting_groups',
                palette=['black','violet'],
                linewidth=0,
                s=45,
                legend=None)



# Add labels and title to the plot
g.ax.set_xlabel('ln(eGFR) beta')
g.ax.set_ylabel('ln(UACR) beta')
g.ax.set_title('The associated loci in the general population')
g.ax.axhline(0, color='black', lw=0.5)
g.ax.axvline(0, color='black', lw=0.5)

# Set x and y axis limits
g.ax.set_xlim(-0.0085, 0.0075)
g.ax.set_ylim(-0.025, 0.025)
g.ax.grid(True, which='both', linestyle='--', linewidth=0.5)

# Define legend labels and colors
legend_labels = ['Significance level: 0.05', 'Bonferroni corrected: 0.000137']
legend_markers = ['o', 'o']  # Specify the markers for each label
legend_colors = ['black', 'm']

handles = [plt.Line2D([], [], color=color, marker=marker, linestyle='None') for color, marker in zip(legend_colors, legend_markers)]


g.ax.legend(handles, legend_labels, loc='upper right', markerscale=1, fontsize=12, frameon=True, edgecolor='black')

# save the plot
save_path = '/Users/cilianaaman/Documents/Article/images/scatterplot_eGFR_2021_associated_with_UACR_a.jpg'
plt.savefig(save_path, dpi=300, bbox_inches='tight')

plt.show()


############################################################ 
####  Scatterplot eGFR 2021 rsID association with UACR ##### only bonferroni corrected
############################################################

Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'], errors='coerce')
Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'], errors='coerce')


Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = 'A'

# Iterate over the rows and check for matching rsIDs
# Extract base rsIDs from the Bonferroni dataset
uacr_bonferroni_rsIDs = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2['rsID'].str.split('_').str[0]

# Assign 'B' to rows where rsID_eGFR matches a Bonferroni rsID
Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] = np.where(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant['rsID_eGFR_2021'].isin(uacr_bonferroni_rsIDs),
    'B',
    Scatter_plot_eGFR_2021_associated_with_UACR_significant.get('plotting_groups', 'A')  # Default to 'A'
)
filtered_df = Scatter_plot_eGFR_2021_associated_with_UACR_significant[Scatter_plot_eGFR_2021_associated_with_UACR_significant['plotting_groups'] == 'B'].copy()



# Create the scatter plot
g = sns.relplot(data=filtered_df,
                x='Beta_eGFR_2021',
                y='Beta_UACR',
                aspect=2,
                hue='plotting_groups',
                palette=['violet'],
                linewidth=0,
                s=45,
                legend=None)



# Add labels and title to the plot
g.ax.set_xlabel('ln(eGFR) beta')
g.ax.set_ylabel('ln(UACR) beta')
g.ax.set_title('The associated loci in the general population')
g.ax.axhline(0, color='black', lw=0.5)
g.ax.axvline(0, color='black', lw=0.5)
g.ax.grid(True, which='both', linestyle='--', linewidth=0.5)
# Set x and y axis limits
g.ax.set_xlim(-0.0085, 0.0075)
g.ax.set_ylim(-0.025, 0.025)


# Create the annotations
Anno = filtered_df.apply(
    lambda p: g.ax.annotate(p['Gene'], (p['Beta_eGFR_2021'], p['Beta_UACR']), fontsize=11),
    axis=1
).to_list()

adjust_text(Anno, arrowprops=dict(arrowstyle="-", color='k', lw=0.4), force_points=0.1, force_text=0.3)




# Define legend labels and colors
legend_labels = ['Bonferroni corrected: 0.000137']
legend_markers = ['o']  # Specify the markers for each label
legend_colors = ['m']

handles = [plt.Line2D([], [], color=color, marker=marker, linestyle='None') for color, marker in zip(legend_colors, legend_markers)]


g.ax.legend(handles, legend_labels, loc='upper right', markerscale=1, fontsize=12, frameon=True, edgecolor='black')

# save the plot
save_path = '/Users/cilianaaman/Documents/Article/images/scatterplot_eGFR_2021_associated_with_UACR_b.jpg'
plt.savefig(save_path, dpi=300, bbox_inches='tight')


plt.show()

#%% Doing it only for diabetes 

#####################################
### creating diabetes matrix data ###
#####################################

# Load the diabetes csv data using the Pandas library
df_T2D = pd.read_csv('/Users/cilianaaman/Desktop/df_preprocessed.csv')

# Select only the 'eid' and 'T2D' columns
df_T2D = df_T2D[['eid', 'T2D']]

Overall_matrix_T2D = pd.merge(Overall_matrix, df_T2D, on='eid', how='inner')

# Filter out rows where 'Eid' in Overall_matrix matches 'Eid' in withdraw
Overall_matrix_T2D = Overall_matrix_T2D[~Overall_matrix_T2D['eid'].isin(withdraw['Eid'])]

# only the people with T2D = 1

Overall_matrix_T2D = Overall_matrix_T2D[Overall_matrix_T2D['T2D'] == 1]




##########################
### Preprocessing data ###
##########################


eid_column = Overall_matrix_T2D['eid']

# Select the columns with phenotype data
phenotype_columns = Overall_matrix_T2D.iloc[:, 1:6]

# Calculate the mean and standard deviation for each phenotype column excluding NaN values
mean_values = np.nanmean(phenotype_columns, axis=0)
std_values = np.nanstd(phenotype_columns, axis=0)

# Define the lower and upper bounds for outlier removal
lower_bounds = mean_values - 4 * std_values
upper_bounds = mean_values + 4 * std_values

# Identify the rows where any phenotype value is outside the bounds, excluding NaN values
outlier_rows = np.any((phenotype_columns < lower_bounds) | (phenotype_columns > upper_bounds), axis=1) & np.all(~np.isnan(phenotype_columns), axis=1)

# Remove the outlier rows from the Overall_matrix matrix
Overall_matrix_T2D_filtered = Overall_matrix_T2D[~outlier_rows].copy()

Overall_matrix_T2D_filtered.loc[:, 'eid'] = eid_column.loc[~outlier_rows].values

# Convert UACR to ln(UACR) and eGFR to ln(eGFR)
Overall_matrix_T2D_filtered['micro albumin in urin'] = np.log(Overall_matrix_T2D_filtered['micro albumin in urin'])
Overall_matrix_T2D_filtered['eGFR_2021'] = np.log(Overall_matrix_T2D_filtered['eGFR_2021'])




columns_to_scale = ['Sex', 'Systolic blood pressure', 'Age at recruitment']

# Select all columns that start with 'rs'
rs_columns = [col for col in Overall_matrix_T2D_filtered.columns if col.startswith('rs')]

# Combine the columns to scale and the rs columns
columns_to_scale.extend(rs_columns)

df_selected_columns = Overall_matrix_T2D_filtered[columns_to_scale]

# Scale the selected columns
scaled_columns = preprocessing.scale(df_selected_columns)

Overall_matrix_stand_T2D = pd.DataFrame(index=Overall_matrix_T2D_filtered.index, data=scaled_columns, columns=columns_to_scale)

# Include the columns that were not scaled
non_scaled_columns = [col for col in Overall_matrix_T2D_filtered.columns if col not in columns_to_scale]
Overall_matrix_stand_T2D[non_scaled_columns] = Overall_matrix_T2D_filtered[non_scaled_columns]


columns_to_drop = ['IID', 'SEX']
Overall_matrix_stand_T2D.drop(columns_to_drop, axis=1, inplace=True)





###############################
### MLR 1 for eGFR 2021 T2D ###
###############################


#####################
### for all pheno ###
#####################





# Select the columns with '_eGFR' at the end and include 'eid', 'Sex', 'age' and SBP
columns_to_select = [col for col in Overall_matrix_stand_T2D.columns if col.endswith('_eGFR')] + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure','eGFR_2021']



# Create a new DataFrame with the selected columns
new_matrix = Overall_matrix_stand_T2D[columns_to_select].copy()

# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)


y_log = new_matrix['eGFR_2021']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex ans SBP

# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
#bonferroni_threshold_T2D = 0.05 / len(rs_columns)  # Bonferroni corrected threshold

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    p_value = result.pvalues[rsID]
    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
MLR_eGFR_2021_all_T2D = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value','CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
#MLR_eGFR_2021_significant_results_bonferroni_all_T2D = MLR_eGFR_2021_all_T2D[MLR_eGFR_2021_all_T2D['P-value'] <= bonferroni_threshold_T2D]

# Filter non-significant rsIDs based on a significance level of 0.05
MLR_eGFR_2021_significant_results_all_T2D = MLR_eGFR_2021_all_T2D[MLR_eGFR_2021_all_T2D['P-value'] <= 0.05]



variables_to_delete = ['model','p_value','result','results','rs_columns',
                       'rsID','X','X_with_rsID','y_log']

for var_name in variables_to_delete:
    del locals()[var_name]

del var_name
del variables_to_delete


###########################################
### eGFR 2021 rsID associated with UACR ### 
###########################################

rsIDs_to_select = MLR_eGFR_2021_significant_results_all_T2D['rsID'].tolist()

new_matrix = Overall_matrix_stand_T2D[rsIDs_to_select + ['eid', 'Sex', 'Age at recruitment','Systolic blood pressure','micro albumin in urin']].copy()

# Remove rows with NaN values
new_matrix.dropna(inplace=True)

new_matrix.reset_index(drop=True, inplace=True)

y_log = new_matrix['micro albumin in urin']

# Define the independent variables (X) as the phenotype data and the rsIDs
X = new_matrix[['Sex', 'Age at recruitment','Systolic blood pressure']]  # Selecting Age and Sex


# Perform individual regressions for each rsID
rs_columns = [col for col in new_matrix.columns if col.startswith('rs')]
results = []
bonferroni_threshold_eGFR_2021_associated_with_UACR_T2D = 0.05 / len(rs_columns)  # Bonferroni corrected threshold

for rsID in rs_columns:
    # Add the current rsID column to X
    X_with_rsID = X.join(new_matrix[rsID])

    # Add a constant column to X_with_rsID (for the intercept)
    X_with_rsID = sm.add_constant(X_with_rsID)

    y_log.reset_index(drop=True, inplace=True)

    # Fit the regression model
    model = sm.OLS(y_log, X_with_rsID)
    result = model.fit()
    p_value = result.pvalues[rsID]
    conf_int = result.conf_int().loc[rsID]  # Confidence interval for this rsID


    results.append((rsID, result.params[rsID], p_value, conf_int[0], conf_int[1]))

# Create a DataFrame to store the results
eGFR_2021_rsid_associated_with_UACR_model2_T2D = pd.DataFrame(results, columns=['rsID', 'Effect Size', 'P-value', 'CI Lower', 'CI Upper'])

# Filter non-significant rsIDs based on Bonferroni correction
eGFR_2021_rsid_associated_with_UACR_bonferroni_model2_T2D = eGFR_2021_rsid_associated_with_UACR_model2_T2D[eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'] <= bonferroni_threshold_eGFR_2021_associated_with_UACR_T2D]

# Filter non-significant rsIDs based on a significance level of 0.05
eGFR_2021_rsid_associated_with_UACR_significant_model2_T2D = eGFR_2021_rsid_associated_with_UACR_model2_T2D[eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'] <= 0.05]



##############################################################
### eGFR_2021_rsid_associated_with_UACR_significant_model2_T2D ###
##############################################################

rsID = eGFR_2021_rsid_associated_with_UACR_model2_T2D['rsID'].str.split('_', expand=True)[0]
EA = eGFR_2021_rsid_associated_with_UACR_model2_T2D['rsID'].str.split('_', expand=True)[1]
chr_values = []
gene_values = []


for rs in rsID:
    row = eGFR_matrix_overall[eGFR_matrix_overall['rsID'].str.startswith(rs)]
    chr_values.append(row['chr'].values[0])
    gene_values.append(row['Gene'].values[0])

data = {
    'rsID_eGFR_2021': rsID,
    'chr': chr_values,
    'Gene': gene_values,
    'EA': EA,
    'Beta_UACR': eGFR_2021_rsid_associated_with_UACR_model2_T2D['Effect Size'],
    'p_val_UACR': eGFR_2021_rsid_associated_with_UACR_model2_T2D['P-value'],
    'UACR CI Lower': eGFR_2021_rsid_associated_with_UACR_model2_T2D['CI Lower'],
    'UACR CI Upper': eGFR_2021_rsid_associated_with_UACR_model2_T2D['CI Upper']
}

MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D = pd.DataFrame(data)



#########################
### UACR X eGFR 2021  ### eGFR
#########################


Scatter_plot_eGFR_2021_associated_with_UACR_T2D = pd.DataFrame()

Scatter_plot_eGFR_2021_associated_with_UACR_T2D['rsID_eGFR_2021'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['rsID_eGFR_2021']
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['chr'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['chr']
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['Gene'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['Gene']
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['EA'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['EA']
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['Beta_UACR'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['Beta_UACR']
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['p_val_UACR'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['p_val_UACR']
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['UACR CI Lower'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['UACR CI Lower']
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['UACR CI Upper'] = MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['UACR CI Upper']



Scatter_plot_eGFR_2021_associated_with_UACR_T2D['Beta_eGFR_2021'] = np.nan
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['p_val_eGFR_2021'] = np.nan
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['eGFR CI Lower'] = np.nan
Scatter_plot_eGFR_2021_associated_with_UACR_T2D['eGFR CI Upper'] = np.nan



for idx, row in MLR_eGFR_2021_all_T2D.iterrows():
    rsID_prefix = row['rsID'].split('_')[0]
    if rsID_prefix in MLR2_matrix_eGFR_2021_rsid_associated_with_UACR_model2_T2D['rsID_eGFR_2021'].values:
        matching_row = Scatter_plot_eGFR_2021_associated_with_UACR_T2D[Scatter_plot_eGFR_2021_associated_with_UACR_T2D['rsID_eGFR_2021'] == rsID_prefix]
        Scatter_plot_eGFR_2021_associated_with_UACR_T2D.loc[matching_row.index, 'Beta_eGFR_2021'] = row['Effect Size']
        Scatter_plot_eGFR_2021_associated_with_UACR_T2D.loc[matching_row.index, 'p_val_eGFR_2021'] = row['P-value']
        Scatter_plot_eGFR_2021_associated_with_UACR_T2D.loc[matching_row.index, 'eGFR CI Lower'] = row['CI Lower']
        Scatter_plot_eGFR_2021_associated_with_UACR_T2D.loc[matching_row.index, 'eGFR CI Upper'] = row['CI Upper']






############################################
### filtereing only the significant data ###
############################################



Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D = Scatter_plot_eGFR_2021_associated_with_UACR_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_T2D['p_val_eGFR_2021'] <=  0.05) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_T2D['p_val_UACR'] <=  0.05)
].reset_index(drop=True)


#%%

############################################################
####  Scatterplot eGFR 2021 rsID association with UACR ##### all without names 
############################################################

Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'], errors='coerce')
Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] = pd.to_numeric(Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'], errors='coerce')

Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['plotting_groups'] = 'A'
# Extract base rsIDs from the Bonferroni dataset
uacr_bonferroni_rsIDs_T2D = eGFR_2021_rsid_associated_with_UACR_bonferroni_model2_T2D['rsID'].str.split('_').str[0]

# Assign 'B' to rows where rsID_eGFR matches a Bonferroni rsID
Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['plotting_groups'] = np.where(
    Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['rsID_eGFR_2021'].isin(uacr_bonferroni_rsIDs_T2D),
    'B',
    Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.get('plotting_groups', 'A')  # Default to 'A'
)


# Create the scatter plot 
g = sns.relplot(data=Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D,
                x='Beta_eGFR_2021',
                y='Beta_UACR',
                aspect=2,
                hue='plotting_groups',
                palette=['black', 'violet'],
                linewidth=0,
                s=45,
                legend=None)


# Add labels and title to the plot
g.ax.set_xlabel('ln(eGFR) beta')
g.ax.set_ylabel('ln(UACR) beta')
g.ax.set_title('The associated loci in the T2D population')
g.ax.axhline(0, color='black', lw=0.5)
g.ax.axvline(0, color='black', lw=0.5)

# Set x and y axis limits
g.ax.set_xlim(-0.0075, 0.0075)
g.ax.set_ylim(-0.05, 0.05)

# Create the annotations
Anno = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.apply(
    lambda p: g.ax.annotate(p['Gene'], (p['Beta_eGFR_2021'], p['Beta_UACR']), fontsize=11),
    axis=1
).to_list()

adjust_text(Anno, arrowprops=dict(arrowstyle="-", color='k', lw=0.4), force_points=0.1, force_text=0.3)





# Define legend labels and colors
legend_labels = ['Significance level: 0.05',
                 'Bonferroni corrected: 0.00031']
legend_markers = ['o','o']  # Specify the markers for each label
legend_colors = ['black','violet']

handles = [plt.Line2D([], [], color=color, marker=marker, linestyle='None') for color, marker in zip(legend_colors, legend_markers)]

g.ax.grid(True, which='both', linestyle='--', linewidth=0.5)
g.ax.legend(handles, legend_labels, loc='upper right', markerscale=1, fontsize=12, frameon=True, edgecolor='black')

# save the plot
save_path = '/Users/cilianaaman/Documents/Article/images/scatterplot_eGFR_2021_associated_with_UACR_T2D.jpg'
plt.savefig(save_path, dpi=300, bbox_inches='tight')


plt.show()


#%% count how many points are in each quadrent 

###############
### For T2D ###
###############

# Counting points in the 1. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] > 0)
].shape[0]

print(f"Number of points in the 1. quadrant {count}")

#  Counting points in the 2. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] < 0)
].shape[0]

print(f"Number of points in the 2. quadrant: {count}")

# Counting points in the 3. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] < 0)
].shape[0]

print(f"Number of points in the 3. quadrant: {count}")

# Counting points in the 4. quadrant 
count = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D['Beta_eGFR_2021'] > 0)
].shape[0]

print(f"Number of points in the 4. quadrant: {count}")


###################
### For overall ###
###################

# Counting points in the 1st quadrant and their occurrences in plotting_groups 'A' and 'B'
count_1st_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] > 0)
]
count_1st_quadrant_A = count_1st_quadrant[
    count_1st_quadrant['plotting_groups'] == 'A'
].shape[0]
count_1st_quadrant_B = count_1st_quadrant[
    count_1st_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 1st quadrant: {count_1st_quadrant.shape[0]}")
print(f"Number of points in the 1st quadrant with 'A' in plotting_groups: {count_1st_quadrant_A}")
print(f"Number of points in the 1st quadrant with 'B' in plotting_groups: {count_1st_quadrant_B}")

# Counting points in the 2nd quadrant and their occurrences in plotting_groups 'A' and 'B'
count_2nd_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] > 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] < 0)
]
count_2nd_quadrant_A = count_2nd_quadrant[
    count_2nd_quadrant['plotting_groups'] == 'A'
].shape[0]
count_2nd_quadrant_B = count_2nd_quadrant[
    count_2nd_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 2nd quadrant: {count_2nd_quadrant.shape[0]}")
print(f"Number of points in the 2nd quadrant with 'A' in plotting_groups: {count_2nd_quadrant_A}")
print(f"Number of points in the 2nd quadrant with 'B' in plotting_groups: {count_2nd_quadrant_B}")

# Counting points in the 3rd quadrant and their occurrences in plotting_groups 'A' and 'B'
count_3rd_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] < 0)
]
count_3rd_quadrant_A = count_3rd_quadrant[
    count_3rd_quadrant['plotting_groups'] == 'A'
].shape[0]
count_3rd_quadrant_B = count_3rd_quadrant[
    count_3rd_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 3rd quadrant: {count_3rd_quadrant.shape[0]}")
print(f"Number of points in the 3rd quadrant with 'A' in plotting_groups: {count_3rd_quadrant_A}")
print(f"Number of points in the 3rd quadrant with 'B' in plotting_groups: {count_3rd_quadrant_B}")

# Counting points in the 4th quadrant and their occurrences in plotting_groups 'A' and 'B'
count_4th_quadrant = Scatter_plot_eGFR_2021_associated_with_UACR_significant[
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_UACR'] < 0) &
    (Scatter_plot_eGFR_2021_associated_with_UACR_significant['Beta_eGFR_2021'] > 0)
]
count_4th_quadrant_A = count_4th_quadrant[
    count_4th_quadrant['plotting_groups'] == 'A'
].shape[0]
count_4th_quadrant_B = count_4th_quadrant[
    count_4th_quadrant['plotting_groups'] == 'B'
].shape[0]
print(f"Number of points in the 4th quadrant: {count_4th_quadrant.shape[0]}")
print(f"Number of points in the 4th quadrant with 'A' in plotting_groups: {count_4th_quadrant_A}")
print(f"Number of points in the 4th quadrant with 'B' in plotting_groups: {count_4th_quadrant_B}")

#%%
# Sort the DataFrame by both p_val_UACR and p_val_eGFR_2021 in ascending order
sorted_df = Scatter_plot_eGFR_2021_associated_with_UACR_significant.sort_values(by=['p_val_UACR', 'p_val_eGFR_2021'], ascending=True)

# Select the top 5 rows
top_5_rows = sorted_df.head(20)

# Print the values in the rsID_eGFR_2021, p_val_UACR, and p_val_eGFR_2021 columns for the top 5 rows
print(top_5_rows[['rsID_eGFR_2021', 'p_val_UACR', 'p_val_eGFR_2021']])

# Sort the DataFrame by both p_val_UACR and p_val_eGFR_2021 in ascending order
sorted_df_T2D = Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.sort_values(by=['p_val_UACR', 'p_val_eGFR_2021'], ascending=True)

# Select the top 5 rows
top_5_rows = sorted_df_T2D.head(10)

# Print the values in the rsID_eGFR_2021, p_val_UACR, and p_val_eGFR_2021 columns for the top 5 rows
print(top_5_rows[['rsID_eGFR_2021', 'p_val_UACR', 'p_val_eGFR_2021']])

#%%


# Save to a specific folder (replace the path with your desired folder)
#output_path = '/Users/cilianaaman/Documents/Article/Data/Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.xlsx'

# Export to Excel
#Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D.to_excel(output_path, index=False, header=True)

# Save to a specific folder (replace the path with your desired folder)
#output_path2 = '/Users/cilianaaman/Documents/Article/Data/Scatter_plot_eGFR_2021_associated_with_UACR_significant.xlsx'

# Export to Excel
#Scatter_plot_eGFR_2021_associated_with_UACR_significant.to_excel(output_path2, index=False, header=True)

#%% Result table export 




###################################################
### result table significant overall population ###
###################################################

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    eGFR_matrix_overall,
    Scatter_plot_eGFR_2021_associated_with_UACR_significant,
    left_on='rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table = result_table[[
    'rsID',             # SNP ID
    'chr_x',
    'Gene_x',
    'EA_x',               # Effect Allele
    'OA',               # Other Allele
    'EAF',              # MAF (minor allele frequency)
    'Beta_eGFR_2021',   # Beta for eGFR
    'p_val_eGFR_2021',  # P-value for eGFR
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',        # Beta for UACR
    'p_val_UACR',        # P-value for UACR
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'

})

# Export the final table to an Excel file
output_file = '/Users/cilianaaman/Documents/Article/result/significant_result_table_overall_population.xlsx'
result_table.to_excel(output_file, index=False)


#############################################
### result table MLR 1 overall population ###
#############################################
# Extract the clean rsID from the MLR_eGFR_2021_all DataFrame
MLR_eGFR_2021_all['clean_rsID'] = MLR_eGFR_2021_all['rsID'].str.split('_').str[0]

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    eGFR_matrix_overall,
    MLR_eGFR_2021_all,
    left_on='rsID',
    right_on='clean_rsID',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table = result_table[[
    'rsID_x',             # SNP ID
    'chr',
    'Gene',
    'EA',               # Effect Allele
    'OA',               # Other Allele
    'EAF',              # MAF (minor allele frequency)
    'Effect Size',   # Beta for eGFR
    'P-value',  # P-value for eGFR
    'CI Lower',
    'CI Upper'
]].rename(columns={
    'rsID_x': 'rsID',
    'chr': 'Chr',
    'EA': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Effect Size': 'Beta ln(eGFR)',
    'P-value': 'P-value eGFR',
    'CI Upper': 'CI Upper eGFR',
    'CI Lower': 'CI Lower eGFR'

})

# Export the final table to an Excel file
output_file = '/Users/cilianaaman/Documents/Article/result/result_table_MLR1_overall_population.xlsx'
result_table.to_excel(output_file, index=False)


#############################################
### result table MLR 2 overall population ###
#############################################

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    eGFR_matrix_overall,
    Scatter_plot_eGFR_2021_associated_with_UACR,
    left_on='rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)


# Select and rename columns for the result table
result_table = result_table[[
    'rsID',             # SNP ID
    'chr_x',
    'Gene_x',
    'EA_x',               # Effect Allele
    'OA',               # Other Allele
    'EAF',              # MAF (minor allele frequency)
    'Beta_eGFR_2021',   # Beta for eGFR
    'p_val_eGFR_2021',  # P-value for eGFR
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',        # Beta for UACR
    'p_val_UACR',        # P-value for UACR
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'
    
})

# Export the final table to an Excel file
output_file = '/Users/cilianaaman/Documents/Article/result/result_table_overall_population.xlsx'
result_table.to_excel(output_file, index=False)



###############################################
### result table significant T2D population ###
###############################################

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    eGFR_matrix_overall,
    Scatter_plot_eGFR_2021_associated_with_UACR_significant_T2D,
    left_on='rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table = result_table[[
    'rsID',             # SNP ID
    'chr_x',
    'Gene_x',
    'EA_x',               # Effect Allele
    'OA',               # Other Allele
    'EAF',              # MAF (minor allele frequency)
    'Beta_eGFR_2021',   # Beta for eGFR
    'p_val_eGFR_2021',  # P-value for eGFR
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',        # Beta for UACR
    'p_val_UACR',        # P-value for UACR
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'
    
})

# Export the final table to an Excel file
output_file = '/Users/cilianaaman/Documents/Article/result/significant_result_table_T2D_population.xlsx'
result_table.to_excel(output_file, index=False)




#############################################
### result table MLR 1 overall population ###
#############################################
MLR_eGFR_2021_all_T2D['clean_rsID'] = MLR_eGFR_2021_all_T2D['rsID'].str.split('_').str[0]

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    eGFR_matrix_overall,
    MLR_eGFR_2021_all_T2D,
    left_on='rsID',
    right_on='clean_rsID',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table = result_table[[
    'rsID_x',            # SNP ID
    'chr',
    'Gene',
    'EA',               # Effect Allele
    'OA',               # Other Allele
    'EAF',              # MAF (minor allele frequency)
    'Effect Size',   # Beta for eGFR
    'P-value',  # P-value for eGFR
    'CI Lower',
    'CI Upper'
]].rename(columns={
    'rsID_x': 'rsID',
    'chr': 'Chr',
    'EA': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Effect Size': 'Beta ln(eGFR)',
    'P-value': 'P-value eGFR',
    'CI Lower': 'CI Lower eGFR',
    'CI Upper': 'CI Upper eGFR'

})

# Export the final table to an Excel file
output_file = '/Users/cilianaaman/Documents/Article/result/result_table_MLR1_T2D_population.xlsx'
result_table.to_excel(output_file, index=False)


#############################################
### result table MLR 2 overall population ###
#############################################

# Merge dataframes on 'rsID' (or the equivalent shared identifier)
result_table = pd.merge(
    eGFR_matrix_overall,
    Scatter_plot_eGFR_2021_associated_with_UACR_T2D,
    left_on='rsID',
    right_on='rsID_eGFR_2021',
    how='inner'  # Use 'inner' to keep only matching rows
)

# Select and rename columns for the result table
result_table = result_table[[
    'rsID',             # SNP ID
    'chr_x',
    'Gene_x',
    'EA_x',               # Effect Allele
    'OA',               # Other Allele
    'EAF',              # MAF (minor allele frequency)
    'Beta_eGFR_2021',   # Beta for eGFR
    'p_val_eGFR_2021',  # P-value for eGFR
    'eGFR CI Lower',
    'eGFR CI Upper',
    'Beta_UACR',        # Beta for UACR
    'p_val_UACR',        # P-value for UACR
    'UACR CI Lower',
    'UACR CI Upper'
]].rename(columns={
    'chr_x': 'Chr',
    'Gene_x': 'Gene',
    'EA_x': 'Effect Allele',
    'OA': 'Other Allele',
    'EAF': 'EAF',
    'Beta_eGFR_2021': 'Beta ln(eGFR)',
    'p_val_eGFR_2021': 'P-value eGFR',
    'Beta_UACR': 'Beta ln(UACR)',
    'p_val_UACR': 'P-value UACR',
    'eGFR CI Lower': 'CI Lower eGFR',
    'eGFR CI Upper': 'CI Upper eGFR',
    'UACR CI Lower': 'CI Lower UACR',
    'UACR CI Upper': 'CI Upper UACR'
    
})

# Export the final table to an Excel file
output_file = '/Users/cilianaaman/Documents/Article/result/result_table_T2D_population.xlsx'
result_table.to_excel(output_file, index=False)







