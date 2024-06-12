import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def add_vlines(df, col_name, surp_id, c_palette):
    
    col_vals = list(df[col_name].unique())
    for val in col_vals:
        df_val = df[df[col_name] == val].copy()
        val_surprisals = df_val[surp_id]
        mean = np.mean(val_surprisals)
        val_idx = col_vals.index(val)
        color = c_palette[val_idx]
        plt.axvline(mean, c=color, linestyle='--')


def kde_plot_conditions(df, study_id, model_id, make_title=True):
    
    if study_id == 'dbc19':
        r_dict= { 'a':'A: Baseline','b':'B: Event-related violation','c':'C: Event-unrelated violation'}
        df_new = df.replace({"Condition":r_dict})
        c_palette = ['black','red','blue']
        title = 'Delogu et al. (2019)'

    if study_id == 'adsbc21':
        r_dict ={'A':'A: A+E+','B':'B: A-E+','C':'C: A+E-', 'D':'D: A-E-'} 
        df_new = df.replace({"Condition":r_dict})
        c_palette = ["#000000", "#BB5566", "#004488", "#DDAA33"]
        title = 'Aurnhammer et al. (2021)'

    if study_id == 'adbc23':
        r_dict ={'A':'A: Plausible','B':'B: Less plausible, attraction','C':'C: Implausible, no attraction'} 
        df_new = df.replace({"Condition":r_dict})
        c_palette = ['black','red','blue']
        title = 'Aurnhammer et al. (2023)'

    if model_id == 'leo13b':
        surp_id = 'leo13b_surp'
        x_lab_name = 'Leo-13b surprisal'
    
    if model_id == 'secret-gpt-2':
        surp_id = 'secretgpt2_surp'
        x_lab_name = 'GPT-2 surprisal'

    sns.set(style='darkgrid')
    plot = sns.kdeplot(data=df_new,
                       x=surp_id,
                       hue='Condition',
                       palette=c_palette,
                       clip=(0,df_new[surp_id].max()),
                       fill=True)
    
    plot.set(xlabel=x_lab_name)

    plot.set(xlabel=None,ylabel=None)
    plot.set_xlim(0,df_new[surp_id].max())

    if make_title:
        plot.set_title(title,fontsize=13,pad=5,x=0.12)

    add_vlines(df, 'Condition', surp_id, c_palette)

    plt.tight_layout()
    plt.savefig(f'../plots/{study_id}/{study_id}_{model_id}_conditions.pdf')
    plt.clf()



if __name__ == '__main__':
    adsbc21_df = pd.read_csv('../data/adsbc21/adsbc21.csv',sep=';')
    dbc19_df = pd.read_csv('../data/dbc19/dbc19.csv',sep=';')
    adbc23_df = pd.read_csv('../data/adbc23/adbc23.csv',sep=';')

    
    kde_plot_conditions(adsbc21_df,'adsbc21','leo13b')