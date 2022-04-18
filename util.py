import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np

def plot(x, y, plot_type='single', plot_mode='line', **setting):
    '''General basic plot setting for plot and subplots
    x, data for X Axis, 1-D array
    y, data for Y Axis, for single plot, 1-D array; for N subplots, N-D array
    
    '''
    default_key = ['title', 'xlabel', 'ylabel', 'plot_name']
    advanced_key = ['xlim', 'ylim']
    default_dict = {'title':'unname', 'xlabel':'X Axis', 'ylabel':'Y Axis', 'plot_name':'unname'}

    if plot_type  == 'single':
        fig, ax = plt.subplots()
        ax.plot(x,y)
        ax.set_xlabel(setting[['xlabel']])
        ax.set_ylabel(setting[['ylabel']])
        ax.set_title(setting['title'])
        fig.savefig(setting['plot_name'])

    elif plot_type == 'sub':
        n_plot = len(y)
        n_col_max = 3
        if n_plot <= 3:
            n_row = 1
            n_col = n_plot
        elif n_plot > 3:
            n_row = n_plot // n_col_max + 1
            n_col = n_col_max
        
        fig, axs = plt.subplots(n_row, n_col, sharex='all', sharey='all')
        axs.set_title(setting['title'])
        for i in range(0, n_row):
            for j in range(0, n_col):
                if plot_mode == 'line':
                    axs[i,j].plot(x,y[i*n_col_max+j])
                    axs[i,j].set_xlabel(setting[['xlabel']])
                    axs[i,j].set_ylabel(setting[['ylabel']])
                elif plot_mode == 'scatter':
                    axs[i,j].scatter(x,y[i*n_col_max+j])
                    axs[i,j].set_xlabel(setting[['xlabel']])
                    axs[i,j].set_ylabel(setting[['ylabel']])                    
                else:
                    raise NotImplemented
        
        fig.savefig(setting['plot_name'])
    else:
        raise NotImplemented

def printm(data):
    '''Print a matrix in a neat  way'''
    return NotImplemented

def printcm(data):
    '''Print a complex matrix in a neat way'''
    return NotImplemented