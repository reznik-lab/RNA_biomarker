import pandas as pd
import scipy.stats
from statistics import mean
import numpy as np
from matplotlib import pyplot as plt
from lifelines import CoxPHFitter
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests

if __name__ == '__main__':
    signature = 'WP_FERROPTOSIS'
    signature = 'H.OXIDATIVE_PHOSPHORYLATION'

    # BraunEtAl
    sub_dir = 'BraunEtAl'
    len_samples = 311
    data_dir = 'revision'
    results_dir = 'revision/results'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'

    # Load the signature score
    score_df = pd.read_csv(f'{data_dir}/{sub_dir}.WP_FERROPTOSIS.ssGSEA.csv', header=0, index_col=0)
    score_df = pd.read_csv(f'{data_dir}/MergedDeconvolution.{sub_dir}.csv', header=0, index_col=0)[signature]

    # Load the clinical data
    clinical = pd.read_excel(f'{clinical_path}.xlsx', sheet_name="Braun_et_al.FOLH1", index_col="RNA_ID")
    clinical = clinical.iloc[:311, :]

    # Merge the data
    data = pd.concat([score_df, clinical], axis=1)
    data = data.iloc[:len_samples, :]

    # Merge datasets together
    cols = []
    cols.extend([signature])
    cols.extend(['Age', 'Sex', 'OS', 'OS_CNSR', 'PFS', 'PFS_CNSR','Arm'])
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
    survival_df = survival_df.drop(columns= ['PFS_CNSR','OS_CNSR'])
    survival_df['Dataset'] = 'BraunEtAl'
    survival_df_T = survival_df.loc[survival_df['Arm']=='NIVOLUMAB']
    survival_df_C = survival_df.loc[survival_df['Arm']=='EVEROLIMUS']

    # javelin_101
    sub_dir = 'javelin_101'
    data_dir = 'revision'
    results_dir = 'revision/results'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'
    len_samples = 726

    # Load the signature score
    score_df = pd.read_csv(f'{data_dir}/{sub_dir}.WP_FERROPTOSIS.ssGSEA.csv', header=0, index_col=0)
    score_df = pd.read_csv(f'{data_dir}/MergedDeconvolution.{sub_dir}.csv', header=0, index_col=0)[signature]

    # Load the clinical data
    clinical = pd.read_excel(f'{clinical_path}.xlsx', sheet_name="all_clinical_data",index_col='ID')

    # Merge the data
    data = pd.concat([score_df, clinical], axis=1)
    data = data.iloc[:len_samples, :]

    data = data.rename({'PFS_P': 'PFS_MO'},axis=1)
    data.loc[data['PFS_P_CNSR'] == 0, 'PFS_EVENT'] = 1
    data.loc[data['PFS_P_CNSR'] == 1, 'PFS_EVENT'] = 0
    cols = []
    cols.extend([signature])
    cols.extend(['AGE','SEX','PFS_MO','PFS_EVENT','TRT01P'])
    survival_df_add = data[cols]
    survival_df_add['Dataset'] = 'javelin_101'
    #survival_df = pd.concat([survival_df, survival_df_add], axis=0)
    survival_df_add_T = survival_df_add.loc[survival_df_add['TRT01P']=='Avelumab+Axitinib']
    survival_df_add_C = survival_df_add.loc[survival_df_add['TRT01P']=='Sunitinib']

    # IMmotion151
    sub_dir = 'IMmotion151'
    data_dir = 'revision'
    results_dir = 'revision/results'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'
    len_samples = 823

    # Load the signature score
    score_df = pd.read_csv(f'{data_dir}/{sub_dir}.WP_FERROPTOSIS.ssGSEA.csv', header=0, index_col=0)
    score_df = pd.read_csv(f'{data_dir}/MergedDeconvolution.{sub_dir}.csv', header=0, index_col=0)[signature]


    # Load the clinical data
    clinical = pd.read_csv(f'{clinical_path}.csv', header=0, index_col='RNASEQ_SAMPLE_ID')
    clinical = clinical.iloc[:823, :]

    # Merge the data
    data = pd.concat([score_df, clinical], axis=1)
    data = data.iloc[:len_samples, :]

    data = data.rename({'PFS_MONTHS': 'PFS_MO'}, axis=1)
    data.loc[data['PFS_CENSOR'] == 0, 'PFS_EVENT'] = 1
    data.loc[data['PFS_CENSOR'] == 1, 'PFS_EVENT'] = 0
    cols = []
    cols.extend([signature])
    cols.extend(['AGE', 'SEX', 'PFS_MO', 'PFS_EVENT','ARM'])
    survival_df_add = data[cols]
    survival_df_add['Dataset'] = 'IMmotion151'
    #survival_df = pd.concat([survival_df, survival_df_add], axis=0)
    survival_df_add_T = survival_df_add.loc[survival_df_add['ARM']=='atezo_bev']
    survival_df_add_C = survival_df_add.loc[survival_df_add['ARM']=='sunitinib']

    # Comparz
    sub_dir = 'Comparz'
    len_samples = 412
    data_dir = 'revision'
    results_dir = 'revision/results'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'


    # Load the signature score
    score_df = pd.read_csv(f'{data_dir}/COMPARZ.WP_FERROPTOSIS.ssGSEA.csv', header=0, index_col=0)
    score_df = pd.read_csv(f'{data_dir}/MergedDeconvolution.COMPARZ.csv', header=0, index_col=0)[signature]

    # Load the clinical data
    clinical = pd.read_excel(f'{clinical_path}.xlsx', sheet_name="Comparz.FOLH1",index_col="RNASampleID")
    clinical = clinical.iloc[25:, :]

    # Merge the data
    data = pd.concat([score_df, clinical], axis=1)
    data = data.iloc[:len_samples, :]
    data = data[data['SRVCFLCD'].notna()]
    data = data[data['SRVMO'].notna()]

    data = data.rename({'SRVCFLCD': 'OS_EVENT'}, axis=1)
    data = data.rename({'SRVMO': 'OS_MO'}, axis=1)
    data = data.rename({'PFSCFLCD': 'PFS_EVENT'}, axis=1)
    data = data.rename({'PFSMO': 'PFS_MO'}, axis=1)

    cols = []
    cols.extend([signature])
    cols.extend(['AGE', 'SEX', 'PFS_MO', 'PFS_EVENT','OS_EVENT','OS_MO','TRTGRP'])
    survival_df_add = data[cols]
    survival_df_add['Dataset'] = 'Comparz'
    #survival_df = pd.concat([survival_df, survival_df_add], axis=0)
    survival_df_add_T = survival_df_add.loc[survival_df_add['TRTGRP']=='pazopanib']
    survival_df_add_C = survival_df_add.loc[survival_df_add['TRTGRP']=='sunitinib']

    # CheckMate214
    sub_dir = 'CheckMate214'
    len_samples = 167
    data_dir = 'revision'
    results_dir = 'revision/results'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'

    # Load the signature score
    score_df = pd.read_csv(f'{data_dir}/{sub_dir}.WP_FERROPTOSIS.ssGSEA.csv', header=0, index_col=0)
    score_df = pd.read_csv(f'{data_dir}/MergedDeconvolution.{sub_dir}.csv', header=0, index_col=0)[signature]

    # Load the clinical data
    clinical = pd.read_csv(f'{clinical_path}.csv', header=0, index_col=0)
    # rename index for checkmate 214
    new_index = []
    for i in range(clinical.shape[0]):
        new_index.append(clinical.index[i].replace('-','.'))
    clinical['new_index'] = new_index
    clinical= clinical.set_index('new_index')

    # Merge the data
    data = pd.concat([score_df, clinical], axis=1)
    data = data.iloc[:len_samples, :]

    data = data.rename({'Age':'AGE'},axis=1)
    data = data.rename({'Sex':'SEX'},axis=1)
    #'Progression Free Survival per Investigator Primary Definition (months) (PFSINV)'
    data = data.rename({'Progression Free Survival per IRRC Primary Definition (months) (PFSIRC)':'PFS_MO'},axis=1)
    #'PFSINV Censor'
    data.loc[data['PFSIRC Censor'] == 0, 'PFS_EVENT'] = 1
    data.loc[data['PFSIRC Censor'] == 1, 'PFS_EVENT'] = 0

    data = data.rename({'Overall Survival (months) (OS)':'OS_MO'},axis=1)
    data.loc[data['OS Censor'] == 0, 'OS_EVENT'] = 1
    data.loc[data['OS Censor'] == 1, 'OS_EVENT'] = 0

    cols = []
    cols.extend([signature])
    cols.extend(['AGE', 'SEX', 'PFS_MO', 'PFS_EVENT','OS_EVENT','OS_MO','Arm'])
    survival_df_add = data[cols]
    survival_df_add['Dataset'] = 'CheckMate214'
    #survival_df = pd.concat([survival_df, survival_df_add], axis=0)
    survival_df_add_T = survival_df_add.loc[survival_df_add['Arm']=='Nivo+Ipi']
    survival_df_add_C = survival_df_add.loc[survival_df_add['Arm']=='Sunitinib']

    # IMmotion150
    sub_dir = 'IMmotion150'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    clinical_path = f'Other_RNA/clinical data/{sub_dir}'
    len_samples = 263

    # Load the predicted MET_RNA matrix
    ranked_predictions = pd.read_csv(f'{results_dir}/predicted_data_ave.csv', header=0, index_col=0)
    met_predicted = ranked_predictions.iloc[:len_samples, :]
    met_predicted = met_predicted.iloc[:, :1573]
    # met_predicted.to_csv(f'{results_dir}/predicted_metabolite_ave.csv')

    # save the overall file
    survival_df.to_csv(f'results_MET_RNA_imputation/all_OS_PFS_262met.csv')

    # Univariate Cox’s proportional hazard model
    # PFS
    cph = CoxPHFitter()
    cph.fit(survival_df_add_T[survival_df_add_T['AGE'].notna()], duration_col='PFS_MO', event_col='PFS_EVENT',
                formula=f'{signature}+AGE+SEX')
    cph.plot()
    plt.savefig(f'{results_dir}/T_{signature}_{sub_dir}_agesex_cox_PFS.pdf')
    met_cox_df = cph.summary
    # met_cox_df.to_csv(f'results_MET_RNA_imputation/all_PFS_CPH_262met.csv')
    met_cox_df.to_csv(f'{results_dir}/T_{signature}_{sub_dir}_agesex_cox_PFS.csv')

    # OS
    met_cox_df = pd.DataFrame(
        {'coef': [], 'exp(coef)': [], 'se(coef)': [], 'coef lower 95%': [], 'coef upper 95%': [],
         'exp(coef) lower 95%': [], 'exp(coef) upper 95%': [], 'cmp to': [], 'z': [], 'p': [], '-log2(p)': []})
    cph = CoxPHFitter()
    cph.fit(survival_df_add_C[survival_df_add_C['OS_MO'].notna()], duration_col='OS_MO', event_col='OS_EVENT',
                formula=f'{i}+AGE+SEX')
    met_cox_df.loc[i] = cph.summary.loc[i]
    met_cox_df.index = met_map.keys()
    met_cox_df['p_adj'] = multipletests(pvals=met_cox_df['p'], method="fdr_bh", alpha=0.1)[1]
    # met_cox_df.to_csv(f'results_MET_RNA_imputation/all_OS_CPH_262met.csv')
    met_cox_df.to_csv(f'{results_dir}/agesex_padj_cox_df_OS_T_262.csv')

    ### Fisher's Method to combine p values
    sub_dir = 'CheckMate214'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    checkmate214 = pd.read_csv(f'{results_dir}/multi_cox_df_PFS_C_35.csv', header=0, index_col=0)
    sub_dir = 'Comparz'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    comparz = pd.read_csv(f'{results_dir}/multi_cox_df_PFS_C_35.csv', header=0, index_col=0)
    sub_dir = 'IMmotion151'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    immotion151 = pd.read_csv(f'{results_dir}/multi_cox_df_PFS_C_35.csv', header=0, index_col=0)
    sub_dir = 'javelin_101'
    results_dir = f'results_MET_RNA_imputation/{sub_dir}'
    javelin = pd.read_csv(f'{results_dir}/multi_cox_df_PFS_C_35.csv', header=0, index_col=0)
    combined_cox_df = pd.DataFrame({'coef': [], 'p': []})
    for i in javelin.index:
        p_list = [javelin.loc[i, 'p'], immotion151.loc[i, 'p'], comparz.loc[i, 'p'], checkmate214.loc[i, 'p']]
        test_statistic, combined_p_value = scipy.stats.combine_pvalues(p_list, method='fisher', weights=None)
        ave_coef = mean(
            [javelin.loc[i, 'coef'], immotion151.loc[i, 'coef'], comparz.loc[i, 'coef'], checkmate214.loc[i, 'coef']])
        combined_cox_df.loc[i, 'coef'] = ave_coef
        combined_cox_df.loc[i, 'p'] = combined_p_value
    combined_cox_df.to_csv(f'results_MET_RNA_imputation/multi_combined_PFS_CPH_C_35met.csv')

    p_list = [0.008857636, 0.243205604, 0.002519368, 0.962589987]
    test_statistic, combined_p_value = scipy.stats.combine_pvalues(p_list, method='fisher', weights=None)

    ########### K-M plot
    # for 1-methylimidazoleacetate
    plt.hist(survival_df_add_T[signature])
    plt.xlabel(f"{signature} scores")
    plt.ylabel("Counts of samples")
    plt.title(f'Histogram of {signature}')
    plt.show()
    print(survival_df_add_T[signature].median())  # 3571.42526722469

    test = survival_df_add_T
    T1 = test[test[signature] < 3571.42526722469]['PFS_MO']
    T2 = test[test[signature] >= 3571.42526722469]['PFS_MO']
    E1 = test[test[signature] < 3571.42526722469]['PFS_EVENT']
    E2 = test[test[signature] >= 3571.42526722469]['PFS_EVENT']
    kmf = KaplanMeierFitter(label=f"High {signature} score")
    kmf.fit(T2, E2)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)
    kmf = KaplanMeierFitter(label=f"Low {signature} score")
    kmf.fit(T1, E1)
    kmf.plot(show_censors=True, ci_show=False)
    print(kmf.median_survival_time_)

    results = logrank_test(T1, T2, event_observed_A=E1, event_observed_B=E2)
    print(results.p_value)
    plt.xlabel("Progression Free Survival(months)")
    plt.ylabel("Survival Probability")
    #plt.show()

    plt.savefig(f'{results_dir}/KM_{signature}_{sub_dir}.pdf')

