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


if __name__ == '__main__':
    signature = ['StromalScore', 'ImmuneScore', 'ESTIMATEScore', 'AdenoSig',
                 'JAVELIN_26', 'Angio_Score', 'Teff_Score', 'Myeloid_Score']
    results_dir = 'revision/results_lennert_project'

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
    #survival_df_T = survival_df.loc[survival_df['ARM']=='atezo_bev']
    #survival_df_C = survival_df.loc[survival_df['ARM']=='sunitinib']
    survival_df_P = survival_df.loc[survival_df['PRIMARY_VS_METASTATIC']=='PRIMARY']
    survival_df_M = survival_df.loc[survival_df['PRIMARY_VS_METASTATIC']=='METASTATIC']



    df = survival_df_C_P
    sub_dir = 'immotion151'
    plot_dir = f'{results_dir}/{sub_dir}/significant'
    # Univariate Cox’s proportional hazard model
    # PFS
    met_cox_df = pd.DataFrame(
        {'coef': [], 'exp(coef)': [], 'se(coef)': [], 'coef lower 95%': [], 'coef upper 95%': [],
         'exp(coef) lower 95%': [], 'exp(coef) upper 95%': [], 'cmp to': [], 'z': [], 'p': [], '-log2(p)': []})
    for i in signature:
        cph = CoxPHFitter()
        cph.fit(df, duration_col='PFS_MO', event_col='PFS_EVENT',
                formula=f'{i}')
        #plt.rcParams['figure.figsize'] = [10, 5]
        #cph.plot()
        #plt.savefig(f'{plot_dir}/{sub_dir}_{i}_uni_PFS_C_M.pdf')
        #plt.close()
        met_cox_df.loc[i]=cph.summary.loc[i]
    met_cox_df.index = signature
    met_cox_df['p_adj'] = multipletests(pvals=met_cox_df['p'], method="fdr_bh", alpha=0.05)[1]
    met_cox_df.to_csv(f'{results_dir}/{sub_dir}/{sub_dir}_coxdf_uni_PFS_T_M.csv')


    # OS
    met_cox_df = pd.DataFrame(
        {'coef': [], 'exp(coef)': [], 'se(coef)': [], 'coef lower 95%': [], 'coef upper 95%': [],
         'exp(coef) lower 95%': [], 'exp(coef) upper 95%': [], 'cmp to': [], 'z': [], 'p': [], '-log2(p)': []})
    for i in signature:
        cph = CoxPHFitter()
        cph.fit(df[df['OS_MO'].notna()], duration_col='OS_MO', event_col='OS_EVENT',
                formula=f'{i}')
        #plt.rcParams['figure.figsize'] = [10, 5]
        #cph.plot()
        #plt.savefig(f'{plot_dir}/{sub_dir}_{i}_uni_OS_C_P.pdf')
        #plt.close()
        met_cox_df.loc[i]=cph.summary.loc[i]
    met_cox_df.index = signature
    met_cox_df['p_adj'] = multipletests(pvals=met_cox_df['p'], method="fdr_bh", alpha=0.05)[1]
    met_cox_df.to_csv(f'{results_dir}/{sub_dir}/{sub_dir}_coxdf_uni_OS_C_M.csv')


    # Multivariate Cox’s proportional hazard model
    # PFS
    met_cox_df = pd.DataFrame(
        {'coef': [], 'exp(coef)': [], 'se(coef)': [], 'coef lower 95%': [], 'coef upper 95%': [],
         'exp(coef) lower 95%': [], 'exp(coef) upper 95%': [], 'cmp to': [], 'z': [], 'p': [], '-log2(p)': []})
    for i in signature:
        cph = CoxPHFitter()
        cph.fit(df[df['AGE'].notna()], duration_col='PFS_MO', event_col='PFS_EVENT',
                formula=f'{i}+AGE+SEX')
        #plt.rcParams['figure.figsize'] = [10, 5]
        #cph.plot()
        #plt.savefig(f'{plot_dir}/{sub_dir}_{i}_multi_PFS_C_P.pdf')
        #plt.close()
        met_cox_df.loc[i]=cph.summary.loc[i]
    met_cox_df.index = signature
    met_cox_df['p_adj'] = multipletests(pvals=met_cox_df['p'], method="fdr_bh", alpha=0.05)[1]
    met_cox_df.to_csv(f'{results_dir}/{sub_dir}/{sub_dir}_coxdf_multi_PFS_T_M.csv')


    # OS
    met_cox_df = pd.DataFrame(
        {'coef': [], 'exp(coef)': [], 'se(coef)': [], 'coef lower 95%': [], 'coef upper 95%': [],
         'exp(coef) lower 95%': [], 'exp(coef) upper 95%': [], 'cmp to': [], 'z': [], 'p': [], '-log2(p)': []})
    for i in signature:
        cph = CoxPHFitter()
        cph.fit(df[df['OS_MO'].notna()][df['AGE'].notna()], duration_col='OS_MO', event_col='OS_EVENT',
                formula=f'{i}+AGE+SEX')
        #plt.rcParams['figure.figsize'] = [10, 5]
        #cph.plot()
        #plt.savefig(f'{plot_dir}/{sub_dir}_{i}_multi_OS_C_P.pdf')
        #plt.close()
        met_cox_df.loc[i]=cph.summary.loc[i]
    met_cox_df.index = signature
    met_cox_df['p_adj'] = multipletests(pvals=met_cox_df['p'], method="fdr_bh", alpha=0.05)[1]
    met_cox_df.to_csv(f'{results_dir}/{sub_dir}/{sub_dir}_coxdf_multi_OS_C_M.csv')


    #  Check significant signatures
    'braunetal'
    'ImmuneScore' and "ESTIMATEScore" in multi_os_c_p and uni_os_c_p
    'StromalScore' in multi_os_c_m
    'Myeloid_Score' in uni_OS_c_m

    'immotion151'
    'Angio_Score' in multi_pfs_c_p
    'StromalScore' in uni_pfs_c_m and multi_pfs_c_m

    ########### K-M plot
    df = survival_df_C_P
    sig = 'Angio_Score'
    type = 'PFS'
    if type == 'OS':
        full_type = 'Overall Survival'
        test = df[df['OS_MO'].notna()]
    elif type == 'PFS':
        full_type = 'Progression Free Survival'
        test = df

    plt.hist(test[sig])
    plt.xlabel(f"{sig} scores")
    plt.ylabel("Counts of samples")
    plt.title(f'Histogram of {sig}')
    plt.show()
    print(test[sig].median())  # 3571.42526722469
    score = test[sig].median()

    T1 = test[test[sig] < score][f'{type}_MO']
    T2 = test[test[sig] >= score][f'{type}_MO']
    E1 = test[test[sig] < score][f'{type}_EVENT']
    E2 = test[test[sig] >= score][f'{type}_EVENT']
    kmf = KaplanMeierFitter(label=f"High {sig} score")
    kmf.fit(T2, E2)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)
    kmf = KaplanMeierFitter(label=f"Low {sig} score")
    kmf.fit(T1, E1)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)

    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    print(results.p_value)
    plt.xlabel(f"{full_type} (months)")
    plt.ylabel("Survival Probability")
    plt.show()

    plt.savefig(f'{plot_dir}/KM_{sig}_{sub_dir}_OS_C_M.pdf')

    #  KM plots for all signatures
    # split into 4 groups
    sub_dir = 'immotion151'
    df = survival_df_C
    arm = 'C'
    type = 'PFS'
    if sub_dir == 'braunetal':
        metastasis = 'Tumor_Sample_Primary_or_Metastas'
    elif sub_dir == 'immotion151':
        metastasis = 'PRIMARY_VS_METASTATIC'
    if type == 'OS':
        full_type = 'Overall Survival'
        df = df[df['OS_MO'].notna()]
    elif type == 'PFS':
        full_type = 'Progression Free Survival'
    km_dir = f'{results_dir}/{sub_dir}/KM_plots'

    for sig in signature:
        median_score = df[sig].median()
        df[f'{sig}_group'] = f'Low {sig}'
        df.loc[df[sig] >= median_score, f'{sig}_group'] = f'High {sig}'
        df['group'] = (df[metastasis] + '_' + df[f'{sig}_group'])
        df['group'] = df['group'].astype('category')

        # fit and plot
        ax = plt.subplot(111)
        plt.rcParams['figure.figsize'] = [10, 10]
        kmf = KaplanMeierFitter()
        for name, group in df.groupby('group'):
            kmf.fit(group[f"{type}_MO"], group[f"{type}_EVENT"], label=name)
            kmf.plot_survival_function(ax=ax, ci_alpha = 0.1,ci_show = False, show_censors= True)
        results = multivariate_logrank_test(df[f"{type}_MO"], df['group'], df[f"{type}_EVENT"])
        print(results.p_value)
        plt.xlabel(f"{full_type} (months)")
        plt.ylabel("Survival Probability")
        plt.title(f'{type} ~ {sig}')
        plt.text(0.95, 0.02, "p = %.3f" % (results.p_value), horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax.transAxes)
        # add_at_risk_counts(kmf, rows_to_show=['At risk'])
        plt.tight_layout()
        plt.savefig(f'{km_dir}/KM_{sub_dir}_{sig}_{type}_{arm}.pdf')
        #plt.show()
        plt.close()

    # Differential tests of signatures
    #-------------------------------------------------------------------------------------------------------------------
    sub_dir = 'braunetal'
    plot_dir = f'{results_dir}/{sub_dir}/differential_plots'
    for sig in signature:
        statistic, p_value = scipy.stats.ranksums(survival_df_P[sig], survival_df_M[sig])
        # Create a violin plot of the two groups
        ax = plt.subplot(111)
        sns.violinplot(data=[survival_df_P[sig], survival_df_M[sig]])
        plt.xticks([0, 1], ['Primary', 'Metastasis'])
        plt.ylabel(f'BraunEtAl {sig}')
        plt.text(0.5, 0.8, "p = %.3f" % (p_value), horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax.transAxes)
        #plt.show()
        plt.savefig(f'{plot_dir}/Violin_{sub_dir}_{sig}_P_M.pdf')
        plt.close()

        # Print the Wilcoxon rank-sum test statistic and p-value
        print('Wilcoxon rank-sum test statistic:', statistic)
        print('Wilcoxon rank-sum test p-value:', p_value)
    # Results showed only Angio_Score have significant differences between primary and metastasis.
    df = survival_df_T_M
    type = 'OS'
    arm = 'T'
    if type == 'OS':
        full_type = 'Overall Survival'
        df = df[df['OS_MO'].notna()]
    elif type == 'PFS':
        full_type = 'Progression Free Survival'
    for sig in signature:
        median_score = df[sig].median()
        df['group'] = f'Low {sig}'
        df.loc[df[sig] >= median_score, 'group'] = f'High {sig}'

        # fit and plot
        ax = plt.subplot(111)
        plt.rcParams['figure.figsize'] = [10, 10]
        kmf = KaplanMeierFitter()
        for name, group in df.groupby('group'):
            kmf.fit(group[f"{type}_MO"], group[f"{type}_EVENT"], label=name)
            kmf.plot_survival_function(ax=ax, ci_alpha = 0.1,ci_show = False, show_censors= True)
        results = multivariate_logrank_test(df[f"{type}_MO"], df['group'], df[f"{type}_EVENT"])
        print(results.p_value)
        plt.xlabel(f"{full_type} (months)")
        plt.ylabel("Survival Probability")
        plt.title(f'{type} ~ {sig} in Metastasis_Treatment arm')
        plt.text(0.92, 0.02, "p = %.3f" % (results.p_value), horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax.transAxes)
        # add_at_risk_counts(kmf, rows_to_show=['At risk'])
        plt.tight_layout()
        plt.savefig(f'{plot_dir}/KM_{sub_dir}_{sig}_{type}_{arm}_M.pdf')
        #plt.show()
        plt.close()

    # None of the signatures have significant differences in the Primary_treatment arm.
    # None of the signatures have significant differences in the Metastasis_treatment arm.

    # cluster analysis in IMmotion151
    #-------------------------------------------------------------------------------------------------------------------
    sub_dir = 'immotion151'
    plot_dir = f'{results_dir}/{sub_dir}/cluster_analysis'
    df = survival_df_M
    df.loc[df['ARM'] == 'sunitinib', 'ARM'] = 'control'
    df.loc[df['ARM'] == 'atezo_bev', 'ARM'] = 'treatment'
    met_cox_df = pd.DataFrame(
        {'coef': [], 'exp(coef)': [], 'se(coef)': [], 'coef lower 95%': [], 'coef upper 95%': [],
         'exp(coef) lower 95%': [], 'exp(coef) upper 95%': [], 'cmp to': [], 'z': [], 'p': [], '-log2(p)': []})
    for i in range(1,8):
        test = df.loc[df['NMF_GROUP'] == i]
        cph = CoxPHFitter()
        try:
            cph.fit(test[test['PFS_MO'].notna()][test['AGE'].notna()], duration_col='PFS_MO', event_col='PFS_EVENT',
                formula='ARM')
        #plt.rcParams['figure.figsize'] = [10, 5]
        #cph.plot()
        #plt.savefig(f'{plot_dir}/{sub_dir}_{i}_multi_OS_C_P.pdf')
        #plt.close()
            met_cox_df.loc[i]=cph.summary.iloc[0]
        except:
            pass
    met_cox_df.index = ['Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7']
    met_cox_df['A/B mPFS'] = np.nan
    met_cox_df['Sunitinib mPFS'] = np.nan

    # combine groups
    #-------------------------------------------------------------------------------------------------------------------
    test = df.loc[(df['NMF_GROUP'] == 2) | (df['NMF_GROUP'] == 1) | (df['NMF_GROUP'] == 7)]
    cph = CoxPHFitter()
    cph.fit(test[test['PFS_MO'].notna()][test['AGE'].notna()], duration_col='PFS_MO', event_col='PFS_EVENT',
            formula='ARM')
    met_cox_df.loc['Cluster 1+2+7'] = cph.summary.iloc[0]

    #K-M ------------------------------------------
    type = 'PFS'
    arm = 'C'
    if arm == 'T':
        df = survival_df_T_M
    elif arm == 'C':
        df = survival_df_C_M
    # fit and plot
    ax = plt.subplot(111)
    plt.rcParams['figure.figsize'] = [10, 10]
    kmf = KaplanMeierFitter()
    i=0
    for name, group in df.groupby('NMF_GROUP'):
        kmf.fit(group[f"{type}_MO"], group[f"{type}_EVENT"], label=name)
        kmf.plot_survival_function(ax=ax, ci_alpha=0.1, ci_show=False, show_censors=True)
        print(kmf.median_survival_time_)
        print(name)
        met_cox_df['Sunitinib mPFS'].iloc[i] = kmf.median_survival_time_
        i+=1
    results = multivariate_logrank_test(df[f"{type}_MO"], df['NMF_GROUP'], df[f"{type}_EVENT"])
    print(results.p_value)
    plt.xlabel("Progression Free Survival (months)")
    plt.ylabel("Survival Probability")
    plt.title(f'Metastasis_control arm')
    plt.text(0.92, 0.02, "p = %.3f" % (results.p_value), horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax.transAxes)
    # add_at_risk_counts(kmf, rows_to_show=['At risk'])
    plt.tight_layout()
    plt.savefig(f'{plot_dir}/KM_{sub_dir}_cluster_{arm}_M.pdf')
    #plt.show()
    plt.close()

    # combine groups
    #------------------------------------------
    group = df.loc[(df['NMF_GROUP'] == 1) | (df['NMF_GROUP'] == 2)]
    kmf = KaplanMeierFitter()
    kmf.fit(group[f"{type}_MO"], group[f"{type}_EVENT"])
    met_cox_df['Sunitinib mPFS'].iloc[7] = kmf.median_survival_time_

    met_cox_df.rename(columns={'exp(coef)': 'HR', 'exp(coef) lower 95%': 'HR lower 95%', 'exp(coef) upper 95%': 'HR upper 95%'}, inplace=True)
    met_cox_df.drop(['coef', 'se(coef)', 'coef lower 95%', 'coef upper 95%', 'cmp to', 'z', '-log2(p)'], axis=1, inplace=True)
    met_cox_df.to_csv(f'{plot_dir}/cluster_analysis_Metastasis.csv')

    # make forest plot
    fig, ax = plt.subplots(figsize=(10, 5))
    for i in range(met_cox_df.shape[0]):
        ax.plot([met_cox_df['HR lower 95%'].iloc[i], met_cox_df['HR upper 95%'].iloc[i]], [i, i], color='k', linewidth=1.5)
        ax.plot(met_cox_df['HR'].iloc[i], i, 'o', color='k', markersize=8)
        ax.axvline(x=1, ymin=i - 0.1, ymax=i + 0.1, color='k', linewidth=0.5)

    # Add labels to plot
    ax.set_yticks(range(met_cox_df.shape[0]))
    ax.set_yticklabels(met_cox_df.index)
    ax.set_xlabel('HR (95% CI)')
    ax.axvline(x=1, color='k', linestyle='--')
    # plt.show()
    plt.savefig(f'{plot_dir}/forest_plot_Metastasis.pdf')
    plt.close()

    # ------------------------May 19----------------------------------
    # IMmotion151
    # (1)
    rna = pd.read_csv('/juno/work/reznik/xiea1/MIRTH/RNA_raw_imputation/IMmotion151.csv', header=0, index_col=0)
    test = pd.concat([rna, survival_df], axis=1)
    met_wilcox_df = pd.DataFrame({'log2fc': [],  'statistic': [], 'p': []})
    for i in list(test.iloc[:,:rna.shape[1]].columns):
        metastasis = test.loc[test['PRIMARY_VS_METASTATIC'] == 'METASTATIC', i]
        primary = test.loc[(test['PRIMARY_VS_METASTATIC'] == 'PRIMARY', i)]
        try:
            met_wilcox_df.loc[i, 'log2fc'] = np.log2(metastasis.mean()/(primary.mean()))
            met_wilcox_df.loc[i, 'statistic'] = scipy.stats.ranksums(metastasis, primary).statistic
            met_wilcox_df.loc[i, 'p'] = scipy.stats.ranksums(metastasis, primary).pvalue
        except:
            pass
    met_wilcox_df['p_adj'] = multipletests(pvals=met_wilcox_df['p'], method="fdr_bh", alpha=0.1)[1]
    met_wilcox_df.to_csv(f'{results_dir}/immotion151/primary_met/differential_genes_statistic.csv')
    print(met_wilcox_df.loc[met_wilcox_df['p_adj'] < 0.1].shape[0]/met_wilcox_df.shape[0])

    # (2)
    sub_dir = 'immotion151'
    plot_dir = f'{results_dir}/{sub_dir}/primary_met'
    for sig in signature:
        statistic, p_value = scipy.stats.ranksums(survival_df_P[sig], survival_df_M[sig])
        # Create a violin plot of the two groups
        ax = plt.subplot(111)
        sns.boxplot(data=[survival_df_P[sig], survival_df_M[sig]])
        plt.xticks([0, 1], ['Primary', 'Metastasis'])
        plt.ylabel(f'IMmotion151 {sig}')
        plt.text(0.5, 0.8, "p = %.3f" % (p_value), horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax.transAxes)
        #plt.show()
        plt.savefig(f'{plot_dir}/Violin_{sub_dir}_{sig}_P_M.pdf')
        plt.close()

        # Print the Wilcoxon rank-sum test statistic and p-value
        print('Wilcoxon rank-sum test statistic:', statistic)
        print('Wilcoxon rank-sum test p-value:', p_value)

    # (3)
    data = {
        'category1': ['Cluster 1', 'Cluster 1', 'Cluster 2', 'Cluster 2', 'Cluster 3', 'Cluster 3',
                      'Cluster 4', 'Cluster 4', 'Cluster 5', 'Cluster 5', 'Cluster 6', 'Cluster 6',
                      'Cluster 7', 'Cluster 7'],
        'category2': ['Primary', 'Metastasis', 'Primary', 'Metastasis', 'Primary', 'Metastasis',
                      'Primary', 'Metastasis', 'Primary', 'Metastasis', 'Primary', 'Metastasis',
                      'Primary', 'Metastasis'],
        'value': [survival_df[survival_df['NMF_GROUP'] == 1][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0],
                  survival_df[survival_df['NMF_GROUP'] == 1][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0],
                  survival_df[survival_df['NMF_GROUP'] == 2][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0],
                  survival_df[survival_df['NMF_GROUP'] == 2][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 3][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 3][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 4][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 4][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 5][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 5][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 6][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 6][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 7][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0],
                    survival_df[survival_df['NMF_GROUP'] == 7][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  ]
    }
    df = pd.DataFrame(data)
    sns.barplot(x='category1', y='value', hue='category2', data=df)
    plt.xlabel('Cluster')
    plt.ylabel('# of Patients')
    plt.title('Frequency of Primary and Metastatic Samples in Clusters')
    plt.savefig(f'{plot_dir}/Cluster_frequency.pdf')
    plt.close()

    # percentage of metastatic samples in each cluster
    data = {
        'category1': ['Cluster 1', 'Cluster 2',  'Cluster 3', 'Cluster 4',  'Cluster 5',  'Cluster 6', 'Cluster 7'],
        'value': [(survival_df[survival_df['NMF_GROUP'] == 1][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  )/(survival_df[survival_df['NMF_GROUP'] == 1][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0]),
                  (survival_df[survival_df['NMF_GROUP'] == 2][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  )/(survival_df[survival_df['NMF_GROUP'] == 2][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0]),
                  (survival_df[survival_df['NMF_GROUP'] == 3][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  )/(survival_df[survival_df['NMF_GROUP'] == 3][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0]),
                  (survival_df[survival_df['NMF_GROUP'] == 4][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  )/(survival_df[survival_df['NMF_GROUP'] == 4][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0]),
                  (survival_df[survival_df['NMF_GROUP'] == 5][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  )/(survival_df[survival_df['NMF_GROUP'] == 5][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0]),
                  (survival_df[survival_df['NMF_GROUP'] == 6][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  )/(survival_df[survival_df['NMF_GROUP'] == 6][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0]),
                  (survival_df[survival_df['NMF_GROUP'] == 7][survival_df['PRIMARY_VS_METASTATIC'] == 'METASTATIC'].shape[0]
                  )/(survival_df[survival_df['NMF_GROUP'] == 7][survival_df['PRIMARY_VS_METASTATIC'] == 'PRIMARY'].shape[0])
                  ]
    }
    df = pd.DataFrame(data).sort_values('value', ascending=False)
    sns.barplot(x='category1', y='value', data=df, color='red')
    plt.xlabel('Cluster')
    plt.ylabel('Percentage of Metastatic Samples')
    plt.title('Percentage of Metastatic Samples in Clusters')
    plt.savefig(f'{plot_dir}/Cluster_metastatic_frequency.pdf')
    plt.close()


    # (4)
    df = survival_df_C
    type = 'PFS'
    if type == 'OS':
        full_type = 'Overall Survival'
        test = df[df['OS_MO'].notna()]
    elif type == 'PFS':
        full_type = 'Progression Free Survival'
        test = df

    T1 = test[test['PRIMARY_VS_METASTATIC'] == 'PRIMARY'][f'{type}_MO']
    T2 = test[test['PRIMARY_VS_METASTATIC'] == 'METASTATIC'][f'{type}_MO']
    E1 = test[test['PRIMARY_VS_METASTATIC'] == 'PRIMARY'][f'{type}_EVENT']
    E2 = test[test['PRIMARY_VS_METASTATIC'] == 'METASTATIC'][f'{type}_EVENT']
    kmf = KaplanMeierFitter(label="Metastasis")
    kmf.fit(T2, E2)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)
    kmf = KaplanMeierFitter(label="Primary")
    kmf.fit(T1, E1)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)

    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    print(results.p_value)
    plt.xlabel(f"{full_type} (months)")
    plt.ylabel("Survival Probability")
    plt.text(0.92, 0.02, "p = %.3f" % (results.p_value), horizontalalignment='center',
                 verticalalignment='center',
                 transform=ax.transAxes)
    plt.show()

    plt.savefig(f'{plot_dir}/KM_{sub_dir}_PFS_control_arm.pdf')
    plt.close()


    # Braunetal
    # (1)
    rna = pd.read_csv('/juno/work/reznik/xiea1/MIRTH/RNA_raw_imputation/BraunEtAl.csv', header=0, index_col=0)
    test = pd.concat([rna, survival_df], axis=1)
    met_wilcox_df = pd.DataFrame({'log2fc': [],  'statistic':[], 'p': []})
    for i in list(test.iloc[:, :rna.shape[1]].columns):
        metastasis = test.loc[test['Tumor_Sample_Primary_or_Metastas']=='METASTASIS', i]
        primary = test.loc[(test['Tumor_Sample_Primary_or_Metastas']=='PRIMARY', i)]
        try:
            met_wilcox_df.loc[i, 'log2fc'] = np.log2(metastasis.mean() / (primary.mean()))
            met_wilcox_df.loc[i, 'statistic'] = scipy.stats.ranksums(metastasis, primary).statistic
            met_wilcox_df.loc[i, 'p'] = scipy.stats.ranksums(metastasis, primary).pvalue
        except:
            pass
    met_wilcox_df['p_adj'] = multipletests(pvals=met_wilcox_df['p'], method="fdr_bh", alpha=0.1)[1]
    met_wilcox_df.to_csv(f'{results_dir}/braunetal/primary_met/differential_genes_statistic.csv')
    print(met_wilcox_df.loc[met_wilcox_df['p_adj'] < 0.1].shape[0] / met_wilcox_df.shape[0])

    # (4)
    sub_dir = 'braunetal'
    df = survival_df_C
    type = 'OS'
    if type == 'OS':
        full_type = 'Overall Survival'
        test = df[df['OS_MO'].notna()]
    elif type == 'PFS':
        full_type = 'Progression Free Survival'
        test = df

    T1 = test[test['Tumor_Sample_Primary_or_Metastas']=='PRIMARY'][f'{type}_MO']
    T2 = test[test['Tumor_Sample_Primary_or_Metastas']=='METASTASIS'][f'{type}_MO']
    E1 = test[test['Tumor_Sample_Primary_or_Metastas']=='PRIMARY'][f'{type}_EVENT']
    E2 = test[test['Tumor_Sample_Primary_or_Metastas']=='METASTASIS'][f'{type}_EVENT']
    kmf = KaplanMeierFitter(label="Metastasis")
    kmf.fit(T2, E2)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)
    kmf = KaplanMeierFitter(label="Primary")
    kmf.fit(T1, E1)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)

    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    print(results.p_value)
    plt.xlabel(f"{full_type} (months)")
    plt.ylabel("Survival Probability")
    plt.text(0.92, 0.02, "p = %.3f" % (results.p_value), horizontalalignment='center',
             verticalalignment='center',
             transform=ax.transAxes)
    plt.show()

    plt.savefig(f'{results_dir}/braunetal/primary_met/KM_{sub_dir}_OS_control_arm.pdf')
    plt.close()

    # --------------------------------------------- June 1 ---------------------------------------------
    # Immotion151
    rna = pd.read_csv('/juno/work/reznik/xiea1/MIRTH/RNA_raw_imputation/IMmotion151.csv', header=0, index_col=0)
    test = pd.concat([rna, survival_df_C_P], axis=1, join='inner')

    met_wilcox_df = pd.DataFrame({'log2fc': [],  'p': []})
    for i in list(test.iloc[:,:rna.shape[1]].columns):
        R = test.loc[(test['OBJECTIVE_RESPONSE'] == 'PR') | (test['OBJECTIVE_RESPONSE'] == 'CR'), i]
        NR = test.loc[(test['OBJECTIVE_RESPONSE'] == 'PD') | (test['OBJECTIVE_RESPONSE'] == 'SD'), i]
        #scipy.stats.ttest_ind(response, non_response)
        try:
            met_wilcox_df.loc[i, 'log2fc'] = np.log2(R.mean()/(NR.mean()))
            met_wilcox_df.loc[i, 'p'] = scipy.stats.ranksums(R, NR).pvalue
        except:
            pass
    met_wilcox_df['p_adj'] = multipletests(pvals=met_wilcox_df['p'], method="fdr_bh", alpha=0.1)[1]
    met_wilcox_df.to_csv(f'{results_dir}/immotion151/primary_met/differential_genes_response_C_P.csv')
    print(met_wilcox_df.loc[met_wilcox_df['p_adj'] < 0.1].shape[0]/met_wilcox_df.shape[0])

    # Braunetal
    rna = pd.read_csv('/juno/work/reznik/xiea1/MIRTH/RNA_raw_imputation/BraunEtAl.csv', header=0, index_col=0)
    test = pd.concat([rna, survival_df_T_P], axis=1, join='inner')

    met_wilcox_df = pd.DataFrame({'log2fc': [],  'p': []})
    for i in list(test.iloc[:,:rna.shape[1]].columns):
        R = test.loc[(test['ORR'] == 'PR') | (test['ORR'] == 'CR') | (test['ORR'] == 'CRPR'), i]
        NR = test.loc[(test['ORR'] == 'PD') | (test['ORR'] == 'SD'), i]
        #scipy.stats.ttest_ind(response, non_response)
        try:
            met_wilcox_df.loc[i, 'log2fc'] = np.log2(R.mean()/(NR.mean()))
            met_wilcox_df.loc[i, 'p'] = scipy.stats.ranksums(R, NR).pvalue
        except:
            pass
    met_wilcox_df['p_adj'] = multipletests(pvals=met_wilcox_df['p'], method="fdr_bh", alpha=0.1)[1]
    met_wilcox_df.to_csv(f'{results_dir}/braunetal/primary_met/differential_genes_response_T_P.csv')
    print(met_wilcox_df.loc[met_wilcox_df['p_adj'] < 0.1].shape[0]/met_wilcox_df.shape[0])
