import os
import pickle
import numpy as np
from tqdm import tqdm

import matplotlib.pyplot as plt

from collections import defaultdict

from fun.selection_of_ROIs import load_data_mm, average_angles, subtract_angle

###########################################################################################################################
#########################################         1. Extract the parameters        ########################################
###########################################################################################################################

def parameters_histograms(MM, xs, ys, path_res):
    """
    generate the histogram for the four parameters

    Parameters
    ----------
    MM : dict
        the dictionnary containing the computed Mueller Matrices
    xs : list of length 2
        the values of the x axis border for the ROI
    ys : list of length 2
        the values of the y axis border for the ROI
    path_res : str
        the path in which to save the histograms
    
    Returns
    ----------
    data_all : dict
        the data for the ROI for each parameter
    """    
    parameters_map = {'linR': ['Linear retardance (°)', False, False, False, False],
    'totP': ['Depolarization', False, True, False, False],
     'totD': ['Diattenuation', False, False, True, False],
     'azimuth': ['Azimuth of optical axis (°)', True, False, False, False]}
    
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,11))
    
    data_all = defaultdict(list)
    keys = list(parameters_map.keys())
    dummy_img = np.zeros(MM[keys[0]].shape)
    
    data_complete = {'linR': MM['linR'],
                    'totP': MM['totP'],
                    'totD': MM['totD'],
                    'azimuth': MM['azimuth']}

    mask_MM = MM['Msk']
    
    for idx, x in enumerate(dummy_img):
        for idy, y in enumerate(x):
            if mask_MM[idx,idy] and xs[0] < idx < xs[1] and ys[0] < idy < ys[1]:
                for key in keys:
                    data_all[key].append(data_complete[key][idx, idy])


    for i, (key, param) in zip(range(0,4), parameters_map.items()):

        row = i%2
        col = i//2
        ax = axes[row, col]

        # change the range of the histograms
        if param[2]:
            range_hist = (0, 1)
        elif param[1]:
            range_hist = (0, 180)
        elif param[3]:
            range_hist = (0, 0.20)
        else:
            range_hist = (0, 60)
        
        data_MM = data_all[key]
        
        y, x = np.histogram(
            data_MM,
            bins=50,
            density=False,
            range = range_hist)
        
        x_plot = []
        for idx, x_ in enumerate(x):
            try: 
                x_plot.append((x[idx] + x[idx + 1]) / 2)
            except:
                pass
                # assert len(x_plot) == 30
        
        # get the mean, max and std
        max_ = x[np.argmax(y)]
        mean = np.nanmean(data_MM)
        std = np.nanstd(data_MM)
        
        y = y / np.max(y)
        
        # plot the histogram
        ax.plot(x_plot,y, c = 'black', linewidth=3)
        ax.axis(ymin=0,ymax=1.5)
        ax.locator_params(axis='y', nbins=4)
        ax.locator_params(axis='x', nbins=5)
    
        if max_:
            ax.text(0.8, 0.9, '$\mu$ = {:.3f}\n$\sigma$ = {:.3f}'.format(mean, std, max_), 
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, 
                    fontsize=16, fontweight = 'bold')
        else:
            ax.text(0.8, 0.9, '$\mu$ = {:.3f}\n$\sigma$ = {:.3f}'.format(mean, std), 
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes, 
                    fontsize=16, fontweight = 'bold')
        
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
            tick.label1.set_fontweight('bold')
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(20)
            tick.label1.set_fontweight('bold')
            
        # ax.set_title(param[0], fontdict = {'fontsize': 30, 'fontweight': 'bold'})
        ax.set_ylabel('Normalized pixel number', fontdict = {'fontsize': 20, 'fontweight': 'bold'})
        
    # save the figures
    plt.tight_layout()
    plt.savefig(os.path.join(path_res, 'parameters_histogram.png'))
    plt.savefig(os.path.join(path_res, 'parameters_histogram.pdf'))
    
    plt.close()
    
    return data_all







