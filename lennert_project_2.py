import pandas as pd
import scipy.stats
from statistics import mean
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from statsmodels.stats.multitest import multipletests
from lifelines.plotting import add_at_risk_counts
from sklearn.decomposition import PCA
from umap import UMAP
import statsmodels.api as sm
import forestplot as fp


if __name__ == '__main__':
    # Load clinical data of IMmotion151 and Braunetal
    signature = ['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'AdenoSig',
                 'JAVELIN_26', 'Angio_Score', 'Teff_Score', 'Myeloid_Score']
    results_dir = 'revision/results_lennert_project'

    # IMmotion151
    len_samples = 823
    data = pd.read_csv(f'{results_dir}/IMmotion151.ClinicalAndDeconvolution.csv', header=0, index_col='RNASEQ_SAMPLE_ID')
    data = data.rename({'PFS_MONTHS': 'PFS_MO'}, axis=1)
    data.loc[data['PFS_CENSOR'] == 0, 'PFS_EVENT'] = 1
    data.loc[data['PFS_CENSOR'] == 1, 'PFS_EVENT'] = 0
    cols = []
    cols.extend(signature)
    cols.extend(['AGE', 'SEX', 'PFS_MO', 'PFS_EVENT', 'ARM', 'PRIMARY_VS_METASTATIC', 'NMF_GROUP', 'MSKCC_RISK_SCORE', 'PDL1_IHC', 'OBJECTIVE_RESPONSE'])
    survival_df = data[cols]
    survival_df['Dataset'] = 'IMmotion151'
    survival_df = survival_df[survival_df['Myeloid_Score'].notna()]
    survival_df['NMF_GROUP'] = survival_df['NMF_GROUP'].astype(int)
    survival_df_T_P = survival_df.loc[(survival_df['ARM']=='atezo_bev') & (survival_df['PRIMARY_VS_METASTATIC']=='PRIMARY')]
    survival_df_T_M = survival_df.loc[(survival_df['ARM']=='atezo_bev') & (survival_df['PRIMARY_VS_METASTATIC']=='METASTATIC')]
    survival_df_C_P = survival_df.loc[(survival_df['ARM']=='sunitinib') & (survival_df['PRIMARY_VS_METASTATIC']=='PRIMARY')]
    survival_df_C_M = survival_df.loc[(survival_df['ARM']=='sunitinib') & (survival_df['PRIMARY_VS_METASTATIC']=='METASTATIC')]
    survival_df_T = survival_df.loc[survival_df['ARM']=='atezo_bev']
    survival_df_C = survival_df.loc[survival_df['ARM']=='sunitinib']
    survival_df_P = survival_df.loc[survival_df['PRIMARY_VS_METASTATIC']=='PRIMARY']
    survival_df_M = survival_df.loc[survival_df['PRIMARY_VS_METASTATIC']=='METASTATIC']

    # BraunEtAl
    len_samples = 311
    data = pd.read_csv(f'{results_dir}/Braun_et_al.ClinicalAndDeconvolution.csv', header=0, index_col='RNA_ID')
    # Unify datasets
    cols = []
    cols.extend(signature)
    cols.extend(['Age', 'Sex', 'OS', 'OS_CNSR', 'PFS', 'PFS_CNSR','Arm', 'Tumor_Sample_Primary_or_Metastas', 'ORR'])
    survival_df = data[cols]
    survival_df = survival_df.rename({'Age': 'AGE'}, axis=1)
    survival_df = survival_df.rename({'Sex': 'SEX'}, axis=1)
    survival_df.loc[survival_df['SEX'] == 'Female', 'SEX'] = 'F'
    survival_df.loc[survival_df['SEX'] == 'FEMALE', 'SEX'] = 'F'
    survival_df.loc[survival_df['SEX'] == 'MALE', 'SEX'] = 'M'
    survival_df.loc[survival_df['SEX'] == 'Male', 'SEX'] = 'M'
    survival_df.loc[survival_df['PFS_CNSR'] == 0, 'PFS_EVENT'] = 1
    survival_df.loc[survival_df['PFS_CNSR'] == 1, 'PFS_EVENT'] = 0
    survival_df.loc[survival_df['OS_CNSR'] == 0, 'OS_EVENT'] = 1
    survival_df.loc[survival_df['OS_CNSR'] == 1, 'OS_EVENT'] = 0
    survival_df = survival_df.rename({'OS': 'OS_MO'}, axis=1)
    survival_df = survival_df.rename({'PFS': 'PFS_MO'}, axis=1)
    # survival_df = survival_df.drop(columns= ['PFS_CNSR','OS_CNSR'])
    survival_df['Dataset'] = 'BraunEtAl'
    survival_df = survival_df[survival_df['Myeloid_Score'].notna()]
    survival_df_T_P = survival_df.loc[(survival_df['Arm']=='NIVOLUMAB') & (survival_df['Tumor_Sample_Primary_or_Metastas']=='PRIMARY')]
    survival_df_T_M = survival_df.loc[(survival_df['Arm']=='NIVOLUMAB') & (survival_df['Tumor_Sample_Primary_or_Metastas']=='METASTASIS')]
    survival_df_C_P = survival_df.loc[(survival_df['Arm']=='EVEROLIMUS') & (survival_df['Tumor_Sample_Primary_or_Metastas']=='PRIMARY')]
    survival_df_C_M = survival_df.loc[(survival_df['Arm']=='EVEROLIMUS') & (survival_df['Tumor_Sample_Primary_or_Metastas']=='METASTASIS')]
    survival_df_T = survival_df.loc[survival_df['Arm']=='NIVOLUMAB' ]
    survival_df_C = survival_df.loc[survival_df['Arm']=='EVEROLIMUS']
    survival_df_P = survival_df.loc[survival_df['Tumor_Sample_Primary_or_Metastas']=='PRIMARY']
    survival_df_M = survival_df.loc[survival_df['Tumor_Sample_Primary_or_Metastas']=='METASTASIS']

    # (1) ---------------------PCA and clustering analysis of RNA profiles (Jul 28)-------------------------
    # IMmotion151
    sub_dir = 'immotion151'
    plot_dir = f'{results_dir}/{sub_dir}/umap'
    metastatic_colname = 'PRIMARY_VS_METASTATIC'
    rna = pd.read_csv('/juno/work/reznik/xiea1/MIRTH/RNA_raw_imputation/IMmotion151.csv', header=0, index_col=0)
    rna = np.log(rna+1)
    # Step 1: Perform PCA on gene expression data
    pca = PCA(n_components=4)
    pca_components = pca.fit_transform(rna)

    # Step 2: Perform UMAP on the PCA components
    #umap = UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean')
    #umap_components = umap.fit_transform(pca_components)

    # Step 3: Create a combined dataframe with UMAP components and metadata
    umap_df = pd.DataFrame(pca_components, columns=['UMAP1', 'UMAP2', 'UMAP3', 'UMAP4'])
    umap_df.index = rna.index
    umap_df = pd.concat([umap_df, survival_df], axis=1)
    umap_df['color'] = umap_df[metastatic_colname].map({'PRIMARY': 0, 'METASTATIC': 1})

    # Step 4: Plot the UMAP results with color-coded metastatic vs. primary samples
    plt.figure(figsize=(8, 8))
    plt.scatter(umap_df['UMAP3'], umap_df['UMAP4'], c=umap_df['color'], cmap='coolwarm', edgecolors='k', s=50)
    plt.xlabel('PCA Component 3')
    plt.ylabel('PCA Component 4')
    plt.title(sub_dir)
    cbar = plt.colorbar(ticks=[0, 1], label='Sample Origin')
    cbar.ax.set_yticklabels(['PRIMARY', 'METASTATIC'])
    plt.savefig(f'{plot_dir}/{sub_dir}_log_rna_pca_34.pdf')
    plt.close()

    # Braunetal
    sub_dir = 'braunetal'
    plot_dir = f'{results_dir}/{sub_dir}/umap'
    metastatic_colname = 'Tumor_Sample_Primary_or_Metastas'
    rna = pd.read_csv('/juno/work/reznik/xiea1/MIRTH/RNA_raw_imputation/BraunEtAl.csv', header=0, index_col=0)
    rna = np.log(rna+1)
    # Step 1: Perform PCA on gene expression data
    pca = PCA(n_components=4)
    pca_components = pca.fit_transform(rna)

    # Step 2: Perform UMAP on the PCA components
    #umap = UMAP(n_neighbors=15, min_dist=0.1, metric='euclidean')
    #umap_components = umap.fit_transform(pca_components)

    # Step 3: Create a combined dataframe with UMAP components and metadata
    umap_df = pd.DataFrame(pca_components, columns=['UMAP1', 'UMAP2', 'UMAP3', 'UMAP4'])
    umap_df.index = rna.index
    umap_df = pd.concat([umap_df, survival_df], axis=1)
    umap_df['color'] = umap_df[metastatic_colname].map({'PRIMARY': 0, 'METASTASIS': 1})

    # Step 4: Plot the UMAP results with color-coded metastatic vs. primary samples
    plt.figure(figsize=(8, 8))
    plt.scatter(umap_df['UMAP3'], umap_df['UMAP4'], c=umap_df['color'], cmap='coolwarm', edgecolors='k', s=50)
    plt.xlabel('PCA Component 3')
    plt.ylabel('PCA Component 4')
    plt.title(sub_dir)
    cbar = plt.colorbar(ticks=[0, 1], label='Sample Origin')
    cbar.ax.set_yticklabels(['PRIMARY', 'METASTATIC'])
    plt.savefig(f'{plot_dir}/{sub_dir}_log_rna_pca_34.pdf')
    plt.close()

    # (2) ---------------------Signature ~ Response Chi Square test (Jul 29)-------------------------
    # IMmotion151
    # treatment arm
    result_df = pd.DataFrame({'Log Odds Ratio': [], 'Lower Bound': [], 'Upper Bound': [], 'Group': []})
    for sig in signature:
        for group, df in {'all': survival_df_T, 'primary': survival_df_T_P, 'meta': survival_df_T_M}.items():
            df = df.loc[df['OBJECTIVE_RESPONSE'] != 'NE']
            df['response'] = np.nan
            df.loc[(df['OBJECTIVE_RESPONSE'] == 'PR') | (df['OBJECTIVE_RESPONSE'] == 'CR'), 'response'] = 1
            df.loc[(df['OBJECTIVE_RESPONSE'] == 'SD') | (df['OBJECTIVE_RESPONSE'] == 'PD'), 'response'] = 0
            median_score = df[sig].median()
            df['signature'] = df[sig].apply(lambda x: 1 if x > median_score else 0)
            df.loc[(df['response']==1) & (df['signature']==1)].shape
            contingency_table = pd.crosstab(df['response'], df['signature'])
            # rows are response (0, 1), columns are signature (0, 1)
            result_df.loc[f'{sig}_{group}', 'Log Odds Ratio'] = sm.stats.Table2x2(contingency_table.values).log_oddsratio
            # log(odds ratio) > 0: higher score in responders
            lower_bound, upper_bound = sm.stats.Table2x2(contingency_table.values).log_oddsratio_confint()
            result_df.loc[f'{sig}_{group}', 'Lower Bound'] = lower_bound
            result_df.loc[f'{sig}_{group}', 'Upper Bound'] = upper_bound
            result_df.loc[f'{sig}_{group}', 'P value'] = sm.stats.Table2x2(contingency_table.values).oddsratio_pvalue()
            result_df.loc[f'{sig}_{group}', 'Group'] = sig

    result_df['Signature'] = result_df.index
    result_df.to_csv(f'{results_dir}/immotion151/umap/chi_square_test_response_immotion151_T.csv')
    fp.forestplot(result_df,
                  estimate="Log Odds Ratio",
                  ll="Lower Bound", hl="Upper Bound",
                  varlabel="Signature",
                  groupvar="Group",
                  pval="P value",
                  ylabel="Confidence interval",
                  xlabel="Log(Odds Ratio)")
    plt.title("IMmotion151: atezo_bev")
    plt.savefig(f'{results_dir}/immotion151/umap/forestplot_immotion151_T.pdf', bbox_inches="tight")

    # control arm
    result_df = pd.DataFrame({'Log Odds Ratio': [], 'Lower Bound': [], 'Upper Bound': [], 'Group': []})
    for sig in signature:
        for group, df in {'all': survival_df_C, 'primary': survival_df_C_P, 'meta': survival_df_C_M}.items():
            df = df.loc[df['OBJECTIVE_RESPONSE'] != 'NE']
            df['response'] = np.nan
            df.loc[(df['OBJECTIVE_RESPONSE'] == 'PR') | (df['OBJECTIVE_RESPONSE'] == 'CR'), 'response'] = 1
            df.loc[(df['OBJECTIVE_RESPONSE'] == 'SD') | (df['OBJECTIVE_RESPONSE'] == 'PD'), 'response'] = 0
            median_score = df[sig].median()
            df['signature'] = df[sig].apply(lambda x: 1 if x > median_score else 0)
            df.loc[(df['response'] == 1) & (df['signature'] == 1)].shape
            contingency_table = pd.crosstab(df['response'], df['signature'])
            # rows are response (0, 1), columns are signature (0, 1)
            result_df.loc[f'{sig}_{group}', 'Log Odds Ratio'] = sm.stats.Table2x2(
                contingency_table.values).log_oddsratio
            # log(odds ratio) > 0: higher score in responders
            lower_bound, upper_bound = sm.stats.Table2x2(contingency_table.values).log_oddsratio_confint()
            result_df.loc[f'{sig}_{group}', 'Lower Bound'] = lower_bound
            result_df.loc[f'{sig}_{group}', 'Upper Bound'] = upper_bound
            result_df.loc[f'{sig}_{group}', 'P value'] = sm.stats.Table2x2(contingency_table.values).oddsratio_pvalue()
            result_df.loc[f'{sig}_{group}', 'Group'] = sig

    result_df['Signature'] = result_df.index
    result_df.to_csv(f'{results_dir}/immotion151/umap/chi_square_test_response_immotion151_C.csv')

    fp.forestplot(result_df,
                  estimate="Log Odds Ratio",
                  ll="Lower Bound", hl="Upper Bound",
                  varlabel="Signature",
                  groupvar="Group",
                  pval="P value",
                  ylabel="Confidence interval",
                  xlabel="Log(Odds Ratio)")
    plt.title("IMmotion151: sunitinib")
    plt.savefig(f'{results_dir}/immotion151/umap/forestplot_immotion151_C.pdf', bbox_inches="tight")

    # Braunetal
    # treatment arm
    result_df = pd.DataFrame({'Log Odds Ratio': [], 'Lower Bound': [], 'Upper Bound': [], 'Group': []})
    for sig in signature:
        for group, df in {'all': survival_df_T, 'primary': survival_df_T_P, 'meta': survival_df_T_M}.items():
            df = df.loc[df['ORR'] != 'NE']
            df['response'] = np.nan
            df.loc[(df['ORR'] == 'PR') | (df['ORR'] == 'CR'), 'response'] = 1
            df.loc[(df['ORR'] == 'SD') | (df['ORR'] == 'PD'), 'response'] = 0
            median_score = df[sig].median()
            df['signature'] = df[sig].apply(lambda x: 1 if x > median_score else 0)
            df.loc[(df['response'] == 1) & (df['signature'] == 1)].shape
            contingency_table = pd.crosstab(df['response'], df['signature'])
            # rows are response (0, 1), columns are signature (0, 1)
            result_df.loc[f'{sig}_{group}', 'Log Odds Ratio'] = sm.stats.Table2x2(
                contingency_table.values).log_oddsratio
            # log(odds ratio) > 0: higher score in responders
            lower_bound, upper_bound = sm.stats.Table2x2(contingency_table.values).log_oddsratio_confint()
            result_df.loc[f'{sig}_{group}', 'Lower Bound'] = lower_bound
            result_df.loc[f'{sig}_{group}', 'Upper Bound'] = upper_bound
            result_df.loc[f'{sig}_{group}', 'P value'] = sm.stats.Table2x2(contingency_table.values).oddsratio_pvalue()
            result_df.loc[f'{sig}_{group}', 'Group'] = sig

    result_df['Signature'] = result_df.index
    result_df.to_csv(f'{results_dir}/braunetal/umap/chi_square_test_response_braunetal_T.csv')
    fp.forestplot(result_df,
                  estimate="Log Odds Ratio",
                  ll="Lower Bound", hl="Upper Bound",
                  varlabel="Signature",
                  groupvar="Group",
                  pval="P value",
                  ylabel="Confidence interval",
                  xlabel="Log(Odds Ratio)")
    plt.title("braunetal: NIVOLUMAB")
    plt.savefig(f'{results_dir}/braunetal/umap/forestplot_braunetal_T.pdf', bbox_inches="tight")

    # control arm
    result_df = pd.DataFrame({'Log Odds Ratio': [], 'Lower Bound': [], 'Upper Bound': [], 'Group': []})
    for sig in signature:
        for group, df in {'all': survival_df_C, 'primary': survival_df_C_P, 'meta': survival_df_C_M}.items():
            df = df.loc[df['ORR'] != 'NE']
            df['response'] = np.nan
            df.loc[(df['ORR'] == 'PR') | (df['ORR'] == 'CR') | (df['ORR'] == 'CRPR'), 'response'] = 1
            df.loc[(df['ORR'] == 'SD') | (df['ORR'] == 'PD'), 'response'] = 0
            median_score = df[sig].median()
            df['signature'] = df[sig].apply(lambda x: 1 if x > median_score else 0)
            df.loc[(df['response'] == 1) & (df['signature'] == 1)].shape
            contingency_table = pd.crosstab(df['response'], df['signature'])
            # rows are response (0, 1), columns are signature (0, 1)
            result_df.loc[f'{sig}_{group}', 'Log Odds Ratio'] = sm.stats.Table2x2(
                contingency_table.values).log_oddsratio
            # log(odds ratio) > 0: higher score in responders
            lower_bound, upper_bound = sm.stats.Table2x2(contingency_table.values).log_oddsratio_confint()
            result_df.loc[f'{sig}_{group}', 'Lower Bound'] = lower_bound
            result_df.loc[f'{sig}_{group}', 'Upper Bound'] = upper_bound
            result_df.loc[f'{sig}_{group}', 'P value'] = sm.stats.Table2x2(contingency_table.values).oddsratio_pvalue()
            result_df.loc[f'{sig}_{group}', 'Group'] = sig

    result_df['Signature'] = result_df.index
    result_df.to_csv(f'{results_dir}/braunetal/umap/chi_square_test_response_braunetal_C.csv')

    fp.forestplot(result_df,
                  estimate="Log Odds Ratio",
                  ll="Lower Bound", hl="Upper Bound",
                  varlabel="Signature",
                  groupvar="Group",
                  pval="P value",
                  ylabel="Confidence interval",
                  xlabel="Log(Odds Ratio)")
    plt.title("braunetal: EVEROLIMUS")
    plt.savefig(f'{results_dir}/braunetal/umap/forestplot_braunetal_C.pdf', bbox_inches="tight")

    # (3) Generate a forest plot for Lennert (Sep. 7, 2023)
    df = pd.read_excel('/juno/work/reznik/xiea1/MIRTH/revision/Table 2a.xlsx',
                       header=0)
    df['Log HR'] = np.log(df['HR'])
    df['Lower Bound'] = np.log(df['Lower Bound'])
    df['Upper Bound'] = np.log(df['Upper Bound'])
    fp.forestplot(df,
                  estimate="Log HR",
                  ll="Lower Bound", hl="Upper Bound",
                  varlabel="Characteristic",
                  groupvar="Group",
                  pval="p-value",
                  ylabel="Confidence interval",
                  xlabel="Log(Hazard Ratio)")
    plt.title("IMmotion151: Atezo_Bev v.s. Sunitinib")
    plt.savefig('revision/results_lennert_project/immotion151/umap/forestplot_immotion151_HR.pdf', bbox_inches="tight")
    plt.close()



