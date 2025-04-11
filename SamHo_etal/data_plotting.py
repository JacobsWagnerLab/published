# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 14:00:32 2025

@author: Alexandros Papagianakis, Christine Jacobs-Wagner lab
"""


import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib

plt.rcParams.update(
    {
        "mathtext.fontset": "stix",
        "font.family": "STIXGeneral",
        "legend.fontsize": 12,  # this is the font size in legends
        "xtick.labelsize": 14,  # this and next are the font of ticks
        "ytick.labelsize": 14,
        "axes.titlesize": 20,
        "axes.labelsize": 20,  # this is the fonx of axes labels
        "savefig.format": "eps",  # how figures should be saved
        "legend.edgecolor": "0.0",
        "legend.framealpha": 0.5,
    }
)



def plot_individual_msd_trajectories(long_msd_df, plot_axis, title, yaxis_limit):
    for traj in long_msd_df.particle_complex.unique():
        traj_df = long_msd_df[long_msd_df.particle_complex==traj]
        plot_axis.plot(traj_df.lag_time_sec, traj_df.msd, color='gray', alpha=0.1)
    plot_axis.set_title(title, fontweight='bold')
    plot_axis.set_xscale('log')
    plot_axis.set_yscale('log')
    plot_axis.set_xlabel(r"$\tau$ (sec)")
    plot_axis.set_ylabel(r"MSD ($\mu$m$^2$)")
    plot_axis.set_ylim(*yaxis_limit)
    


def plot_all_individual_msds(msd_df_list, title_list, yaxis_limit=(10**-6, 5)):
    
    plt.figure(figsize=(len(msd_df_list)*6.7,6))
    gs = GridSpec(1,len(msd_df_list), height_ratios=[1], width_ratios=len(msd_df_list)*[1])
    plot_axis_dict = {}
    for i in range(len(msd_df_list)):
        plot_axis_dict[i] = plt.subplot(gs[0,i])
        plot_individual_msd_trajectories(msd_df_list[i], plot_axis_dict[i], title_list[i], yaxis_limit)
    plt.show()
    
    

def plot_gaussian_classification(msd_intercept_df, title):
    good_df = msd_intercept_df[msd_intercept_df.mixed_gaussian_class_verb=='included']
    bad_df = msd_intercept_df[msd_intercept_df.mixed_gaussian_class_verb=='excluded']
    plt.figure(figsize=(5,5))
    plt.title(title, fontweight='bold')
    # sns.kdeplot(good_df.msd, good_df.alpha2, 'o', label='Gaussian trajectories', color='royalblue')
    # sns.kdeplot(bad_df.msd, bad_df.alpha2, 'o', label='non-Gaussian trajectories', color='tomato')
    plt.plot(good_df.msd, good_df.alpha2, 'o', label='Gaussian trajectories', color='royalblue')
    plt.plot(bad_df.msd, bad_df.alpha2, 'o', label='non-Gaussian trajectories', color='tomato')
    plt.xlabel(r"MSD for $\tau$=50msec ($\mu$m$^2$)")
    plt.ylabel(r"$\alpha$$_{2}$")
    plt.legend(loc='upper left')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()
    
    

def plot_anomalous_exponent(gaussian_df, plot_axis, title, xaxis_limit, yaxis_limit):
    
    kde_fit = gaussian_kde(gaussian_df.alpha)
    plot_axis.hist(gaussian_df.alpha, density=True, histtype='step', 
                   bins=np.arange(xaxis_limit[0], xaxis_limit[1], 0.075), color='black', linewidth=1)
    y_hat = kde_fit(np.arange(xaxis_limit[0], xaxis_limit[1], 0.001))
    plot_axis.plot(np.arange(xaxis_limit[0], xaxis_limit[1], 0.001), y_hat, color='red', linewidth=2)
    plot_axis.set_title(title, fontweight='bold')
    plot_axis.set_xlabel(r"$\alpha$")
    plot_axis.set_ylabel(r"Probability density")
    plot_axis.set_xlim(*xaxis_limit)
    plot_axis.set_ylim(*yaxis_limit)
    
    return y_hat



def plot_all_anomalous_exponents(gaus_df_list, title_list, xaxis_limit, yaxis_limit):
    
    norm = matplotlib.colors.Normalize(vmin=0, vmax=len(gaus_df_list), clip=True)
    mapper = matplotlib.cm.ScalarMappable(norm=norm, cmap='viridis')
    starv_color = np.array([(mapper.to_rgba(v)) for v in list(range(len(gaus_df_list)))])
    
    plt.figure(figsize=(6, len(gaus_df_list)*6.7+6.7))
    gs = GridSpec(len(gaus_df_list)+1, 1, width_ratios=[1], height_ratios=len(gaus_df_list)*[1]+[1])
    plot_axis_dict = {}
    fitted_density_dict = {}
    for i in range(len(gaus_df_list)):
        plot_axis_dict[i] = plt.subplot(gs[i,0])
        fitted_density_dict[i] = plot_anomalous_exponent(gaus_df_list[i], plot_axis_dict[i], title_list[i], 
                                                         xaxis_limit, yaxis_limit)
    plot_axis_dict[i+1] = plt.subplot(gs[i+1,0])
    for cond in fitted_density_dict:
        plot_axis_dict[i+1].plot(np.arange(xaxis_limit[0], xaxis_limit[1], 0.001), 
                                 fitted_density_dict[cond], 
                                 label=title_list[cond], 
                                 color=starv_color[cond],
                                 linewidth=3)
    plot_axis_dict[i+1].legend()
    plot_axis_dict[i+1].set_title('Data comparison', fontweight='bold')
    plot_axis_dict[i+1].set_xlabel(r"$\alpha$")
    plot_axis_dict[i+1].set_ylabel(r"Probability density")
    plt.show()
    
    
    

    
    