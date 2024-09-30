import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import itertools
import os

#######################
#### Density plots ####
#######################

def add_vlines(df, col_name, surp_id, c_palette):

    '''
    Adding vertical mean lines to density plots. 
    '''
    
    col_vals = list(df[col_name].unique())
    for val in col_vals:
        df_val = df[df[col_name] == val].copy()
        val_surprisals = df_val[surp_id]
        mean = np.mean(val_surprisals)
        val_idx = col_vals.index(val)
        color = c_palette[val_idx]
        plt.axvline(mean, c=color, linestyle='--')


def kde_plot_conditions(df, model_id, make_title=True):

    '''
    Generating density plots (kernel density estimation) per study,llm and condition.

    df (DataFrame) : study-specific dataframe containing surprisal values
    model_id (str) : specifies which model's surprisal values will be plotted
    make_title (bool) : whether plot should include study as title
    '''

    study_id = str(df.study_id.unique()[0])

    if study_id == 'dbc19' or study_id == 'dbc19_corrected':
        r_dict= { 'a':'A: Baseline','b':'B: Event-related violation','c':'C: Event-unrelated violation'}
        df_new = df.replace({"Condition":r_dict})
        c_palette = ['black','red','blue']
        title = 'Delogu et al. (2019)'
        x_lim = 44
        y_lim = 0.04

    if study_id == 'adsbc21':
        r_dict ={'A':'A: A+E+','B':'B: A-E+','C':'C: A+E-', 'D':'D: A-E-'} 
        df_new = df.replace({"Condition":r_dict})
        c_palette = ["#000000", "#BB5566", "#004488", "#DDAA33"]
        title = 'Aurnhammer et al. (2021)'
        x_lim = 28
        y_lim = 0.06

    if study_id == 'adbc23':
        r_dict ={'A':'A: Plausible','B':'B: Less plausible, attraction','C':'C: Implausible, no attraction'} 
        df_new = df.replace({"Condition":r_dict})
        c_palette = ['black','red','blue']
        title = 'Aurnhammer et al. (2023)'
        x_lim = 21
        y_lim = 0.2

    if model_id == 'leo13b':
        surp_id = 'leo13b_surp'
        x_lab_name = 'Leo-13b surprisal'
    
    if model_id == 'secretgpt2':
        surp_id = 'secretgpt2_surp'
        x_lab_name = 'secret GPT-2 surprisal'

    if model_id == 'gerpt2':
        surp_id = 'gerpt2_surp'
        x_lab_name = 'GerPT-2 surprisal'

    if model_id == 'gerpt2large':
        surp_id = 'gerpt2large_surp'
        x_lab_name = 'GerPT-2 large surprisal'

    plt.figure(figsize=(4,4))
    sns.set(style='darkgrid')
    plot = sns.kdeplot(data=df_new,
                       x=surp_id,
                       hue='Condition',
                       palette=c_palette,
                       clip=(0,x_lim),
                       fill=True,
                       )
    
                        #clip=(0,df_new[surp_id].max()),
    
    #plot.set(xlabel=x_lab_name)
    plot.set_xlabel(x_lab_name, fontsize = 11)
    # plot.set_xlim(0,df_new[surp_id].max())
    plot.set_xlim(0,x_lim) # 43 was the max surprisal value overall
    plot.set_ylabel("Density", fontsize = 11)
    plot.set_ylim(0,y_lim) 

    if make_title:
        #plot.set_title(title,fontsize=13,pad=5,x=0.12)
        plot.set_title(title)

    add_vlines(df, 'Condition', surp_id, c_palette)
    plt.legend('', frameon=False)
    
    plt.tight_layout()
    plt.savefig(f'../results/{study_id}/plots/{study_id}_{model_id}_conditions.pdf')
    plt.clf()


#########################
#### BPE split check ####
#########################

def check_bpe_splits(df,model_id,bpe_dict):
    
    '''
    '''

    if bpe_dict:
        key = max(bpe_dict.keys())+1
    else:
        key = 0

    study_id = str(df.study_id.unique()[0])
    val_counts = df[f'{model_id}_bpe_split'].value_counts()
    bpe_prop = round(val_counts[1]/val_counts.sum(),3)
    bpe_dict[key] = [study_id, model_id, bpe_prop]
    return bpe_dict

######################
#### Correlations ####
######################

def correlations(df, model_ids):

    method = 'kendall'
    study_id = str(df.study_id.unique()[0])
    surp_ids = [i+'_surp' for i in model_ids]

    if study_id == 'dbc19' or study_id == 'dbc19_corrected':
        df_sub = df[['Cloze', 'Association', 'Plausibility',*surp_ids]]
    
    if study_id == 'adsbc21':
        df_sub = df[['Cloze', 'Association_Noun',*surp_ids]]

    if study_id == 'adbc23':
        df_sub = df[['Cloze', 'Plausibility',*surp_ids]]

    df_sub.corr(method = method).to_csv(f'../results/{study_id}/{study_id}_corr.csv', sep = ';')

###########################
#### Prepare rERP data ####
###########################

def prep_rERP_data(df, model_ids):

    study_id = str(df.study_id.unique()[0])
    print(study_id)
    surp_ids = [i+'_surp' for i in model_ids]

    erp_df = pd.read_csv(f'../data/{study_id}/{study_id}_erp.csv') # load ERP data
    if 'ItemNum' in erp_df.columns:
        erp_df = erp_df.rename(columns={'ItemNum':'Item'})
    print(f'ERP data shape {erp_df.shape}')
    erp_df.set_index(['Item','Condition'],inplace = True) # remove Item & Condition as columns, use as join index
    
    if study_id == 'dbc19' or study_id == 'dbc19_corrected':
        surp_df = df[['Item', 'Condition', 'Cloze', *surp_ids]] # the dbc19 erp data is missing Cloze
    else:
        surp_df = df[['Item','Condition',*surp_ids]]
    surp_df.set_index(['Item','Condition'], inplace=True) # remove Item & Condition as columns, use as join index
    print(f'Surp data shape {surp_df.shape}')
    
    merged_df = erp_df.join(surp_df, how='left')
    merged_df.reset_index(inplace=True) # get back Item & Condition as columns
    print(f'Merged data shape after reset {merged_df.shape}')
    merged_df.to_csv(f'../data/{study_id}/{study_id}_surp_erp.csv', index = False) # new df containing additional surprisal columns 


###################################################################################
###################################################################################


if __name__ == '__main__':

    adsbc21_df = pd.read_csv('../data/adsbc21/adsbc21.csv',sep=';')
    dbc19_df = pd.read_csv('../data/dbc19/dbc19.csv',sep=';')
    adbc23_df = pd.read_csv('../data/adbc23/adbc23.csv',sep=';')
    dbc19_corrected_df = dbc19_df.copy()
    dbc19_corrected_df["study_id"] = "dbc19_corrected"

    study_dfs = [adsbc21_df, dbc19_df, adbc23_df, dbc19_corrected_df]
    study_ids = ['adsbc21', 'dbc19', 'adbc23', 'dbc19_corrected']
    [os.makedirs(f'../results/{study}/plots/', exist_ok = True) for study in study_ids]
    model_ids = ['leo13b', 'secretgpt2','gerpt2', 'gerpt2large']
    make_title = False
    bpe_dict = dict()

    print("Density plots...")
    for args in itertools.product(study_dfs,model_ids):
        kde_plot_conditions(*args,make_title)
        check_bpe_splits(*args,bpe_dict)

    bpe_df = pd.DataFrame.from_dict(bpe_dict,
                                    orient='index',
                                    columns= ['study_id','model_id','bpe_prop'])
    bpe_df.to_csv('../results/bpe_splits.csv', sep = ';', index = False)

    print("Correlations...")
    [correlations(d, model_ids) for d in study_dfs]
    #print("Preparing rERPs...")
    #[prep_rERP_data(d, model_ids) for d in study_dfs]