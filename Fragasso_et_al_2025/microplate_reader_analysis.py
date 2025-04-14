# -*- coding: utf-8 -*-
"""
Created on Mon Apr 14 11:19:56 2025

@author: fragasso
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


'''Function to plot maximum growth rate vs AMP concentration. Uses the computed max_gr_df (see get_max_gr_params function)'''
def plot_max_gr_rate_vs_conc(max_gr_df, time_gr_df, conc, col, cond_name , cond = '', markers = markers, save_path = '', fluor = False):
    fig, ax = plt.subplots(figsize=(7,5.5))
    x_const_spacing = np.flip(np.arange(len(conc)))
    time_gr_arr = np.array(time_gr_df)
    max_gr_arr = np.array(max_gr_df)
    idx_gr_high_global = np.where(max_gr_arr>0.2)
    norm = plt.Normalize(vmin=0, vmax=20)
    for k in range(len(col)):
        marker = markers[k % len(markers)]  # Cycle through markers if there are more plots than markers
        if fluor:
            edge_colors = plt.cm.inferno(norm(time_gr_df[col[k]]))
        else:
            edge_colors = plt.cm.viridis(norm(time_gr_df[col[k]]))
        idx_gr_low = np.where(max_gr_df[col[k]]<0.2)
        edge_colors[idx_gr_low] = np.array([0, 0, 0, 1.0])
        face_colors='none'
        sc = ax.scatter(
            x_const_spacing, 
            max_gr_df[col[k]], 
            edgecolors=edge_colors, 
            facecolors=face_colors,  
            marker=marker, 
            label=col[k], 
            alpha=1,
            s=500/(k+1),
            linewidth=2  
        )
    if fluor:
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='inferno'), ax=ax)
    else: 
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap='viridis'), ax=ax)
    cbar.set_label('Time max gr (hrs)', rotation=90, labelpad=15, fontsize=12)
    ax.set_xticks(x_const_spacing)
    ax.set_xticklabels(np.round(conc,3))
    ax.set_xlabel('Concentration ($\mu$M)', fontsize=12)
    ax.set_ylabel('Max growth rate (log(OD)/hr)', fontsize=12)
    ax.set_title(cond_name)
    ax.set_ylim(-0.05,ax.get_ylim()[1]+0.05)
    legend_handles = []
    for i in range(len(col)):
        if cond == '':
            label=f'rep_{i+1}'
        else: label = cond[i]
        legend_handles.append(plt.Line2D([0], [0], marker=markers[i], color='black',  linestyle='None',markerfacecolor='white', markeredgewidth=2, markersize=10, label=label))
    
    ax.legend(handles=legend_handles)
    
    if len(save_path)>0:
        plt.savefig(save_path + '/'+cond_name+'_max_gr_rate_vs_conc.pdf', format = 'pdf')
    
    plt.show()

def plot_OD_vs_time(df,conc,col,cond,cond_name, row_labels = row_labels, save_path = ''):  
    fig = plt.figure(figsize=(5*len(col), 5))
    fig.suptitle(cond_name, fontsize=12)
    for k in range(len(col)):
        fig.add_subplot(1, len(col), k+1)
        for i in range(len(row_labels)):
            plt.semilogy(df['Time (hrs)'], df[row_labels[i]+str(col[k])],
                     linewidth=1.5, label=str(conc[i]))
        plt.xlabel('Time (hrs)', fontsize=12)
        plt.ylabel('OD (600 nm)', fontsize=12)
        plt.title('Col = '+str(col[k])+', ' + cond[k])
        plt.legend(title='Conc ($\mu$M)')
    plt.tight_layout()
    if len(save_path)>0:
        plt.savefig(save_path + '/'+cond_name+'OD_vs_time.pdf', format = 'pdf')
    plt.show()

def plot_growth_rate_vs_time(df,conc,col,cond,cond_name, row_labels = row_labels, save_path = ''):  
    fig = plt.figure(figsize=(5*len(col), 5))
    fig.suptitle('Growth rate, ' +cond_name, fontsize=12)
    for k in range(len(col)):
        ax = fig.add_subplot(1, len(col), k+1)
        for i in range(len(row_labels)):
            growth_rate=get_growth_rate(np.log(df[row_labels[i]+str(col[k])]),df['Time (hrs)'])
            plt.plot(df['Time (hrs)'], growth_rate,
                     linewidth=1.5, label=str(conc[i]))
        plt.xlabel('Time (hrs)', fontsize=12)
        plt.ylabel('Growth rate (1/hrs)', fontsize=12)
        plt.title('Col = '+str(col[k])+', ' + cond[k])
        plt.legend(title='Conc ($\mu$M)')
    plt.tight_layout()
    if len(save_path)>0:
        plt.savefig(save_path + '/'+cond_name+'growth_rate_vs_time.pdf', format = 'pdf')
    plt.show()

def plot_growth_rate_vs_OD(df,conc,col,cond,cond_name, row_labels = row_labels):      
    fig = plt.figure(figsize=(5*len(col), 5))
    fig.suptitle('Growth rate vs OD ' +cond_name, fontsize=12)
    for k in range(len(col)):
        ax = fig.add_subplot(1, len(col), k+1)
        for i in range(len(row_labels)):
            growth_rate=get_growth_rate(np.log(df[row_labels[i]+str(col[k])]),df['Time (hrs)'])
            plt.plot(df[row_labels[i]+str(col[k])], growth_rate,
                     linewidth=1.5, label=str(conc[i]))
        plt.xlabel('OD', fontsize=12)
        plt.ylabel('Growth rate (1/hrs)', fontsize=12)
        plt.title('Col = '+str(col[k])+', ' + cond[k])
        plt.legend(title='Conc ($\mu$M)')
    plt.tight_layout()
    plt.show()

'''Given input df with stored microplate reader data, row_labels ('A', 'B', etc') and evaluation window t_max (here is 20 hours)
    Compute dataframes with maximum growth rates (max_gr_df) and times at which maximum growth rate was found (time_gr_df)'''
def get_max_gr_params(df, row_labels = row_labels, t_max = 20):
    max_gr_df = pd.DataFrame()
    time_gr_df = pd.DataFrame()
    df = df[df['Time (hrs)']<t_max]
    for k in range(1,13):
        max_gr_rows_list = []
        max_gr_time_list = []
        for i in range(len(row_labels)):
            growth_rate=get_growth_rate(np.log(df[row_labels[i]+str(k)]),df['Time (hrs)'])
            idx = np.argmax(growth_rate)
            max_gr_rows_list.append(np.max(growth_rate))
            max_gr_time_list.append(df['Time (hrs)'].iloc[idx])
        max_gr_df[k] = max_gr_rows_list
        time_gr_df[k] =  max_gr_time_list
    return max_gr_df, time_gr_df

'''Growth rate definition:a 12-point moving window for computing linear regression is applied to the curve (OD,t)'''
def get_growth_rate(signal,time,window_size=12):
    derivatives = np.zeros_like(signal)
    for i in range(window_size, len(signal) - window_size):
        p = np.polyfit(time[i-window_size:i+window_size], signal[i-window_size:i+window_size], 1)
        derivatives[i] = p[0]
    return derivatives

def convert_to_minutes(timestamp):
    hours, minutes, seconds = map(int, timestamp.split(':'))
    return hours * 60 + minutes + seconds / 60


def convert_to_hours(timestamp):
    hours, minutes, seconds = map(int, timestamp.split(':'))
    return hours + minutes / 60 + seconds / 3600


'''Example of how to run the code. Input: text file containing raw microplate reader data'''
row_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
markers = ['s', 'o','^','D' ,  'v']

file_path = 'path_to_raw_data'
df = pd.read_csv(file_path, delimiter=r'\s+', encoding='ISO-8859-1')
time_in_minutes = df['Time'].apply(convert_to_minutes)
time_in_hours = df['Time'].apply(convert_to_hours)
time_col_index = df.columns.get_loc('Time')
df.insert(loc=time_col_index + 1, column='Time (hrs)', value=time_in_hours)

'''
Plot full matrix of microplate reader
'''
m = 8     # number of rows (letter)
n = 12    # number of columns (number)
max_per_column = df.max()
max_v = np.max(max_per_column[3:])
list_wells = list(df)

fig = plt.figure(figsize=(24, 16))  
fig.suptitle(experiment, fontsize=12)
for i in range(m*n):
    ax = fig.add_subplot(m, n, i+1)
    plt.plot(df['Time (hrs)'], df[list_wells[i+3]], linewidth=2)
    plt.ylim([0, max_v + max_v/5])
    plt.xlabel('Time (hrs)', fontsize=12)
    plt.ylabel('OD (600 nm)', fontsize=12)
plt.tight_layout()
plt.show()

'''Get growth rate parameters, used for MIC estimation'''
max_gr_df, time_gr_df = get_max_gr_params(df)

'''
AMP/condition
'''
conc = [64,32,16,8,4,2,1,0] # list of of concentrations, e.g. conc = [64,32,16,8,4,2,1,0]
col = [1,2,3]    # example of columns of microplate reader to analyse and compare (usually these are technical replicates)
cond = ['rep1','rep2','rep3']
cond_name = 'Name of the condition'

plot_OD_vs_time(df,conc,col,cond,cond_name)
plot_growth_rate_vs_time(df,conc,col,cond,cond_name)
plot_growth_rate_vs_OD(df,conc,col,cond,cond_name)
plot_max_gr_rate_vs_conc(max_gr_df, time_gr_df, conc, col, cond_name, cond = cond)


