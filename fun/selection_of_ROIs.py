from fun.border_zone import get_masks

import os
from os.path import exists
import traceback
import pickle
from tqdm import tqdm
import shutil
import imp

import matplotlib.pyplot as plt
import numpy as np
np.random.seed(seed=42)
import pandas as pd
import math

from PIL import Image
from scipy.ndimage import rotate

from datetime import datetime
import seaborn as sns

import matlab.engine


###########################################################################################################################
###########################         1. Set up the parameters and load the Mueller Matrix        ###########################
###########################################################################################################################

def get_mask_matter_and_grid(path_folder, tissue_type):
    _ = get_masks(path_folder)
    path_annotation = os.path.join(path_folder, 'annotation')
    WM_mask = plt.imread(os.path.join(path_annotation, 'WM_merged.png'))
    GM_mask = plt.imread(os.path.join(path_annotation, 'GM_merged.png'))
    grid = np.zeros(WM_mask.shape)
        
    if tissue_type == 'WM':
        mask_matter = WM_mask
    else:
        mask_matter = GM_mask
        
    return mask_matter, grid

def load_data_mm(path_folder, wavelength):
    with open(path_folder + '/polarimetry/' + str(wavelength) + 'nm/MM.pickle', 'rb') as handle:
        mat = pickle.load(handle)
    linear_retardance = load_and_verify_parameters(mat, 'totR')
    diattenuation = load_and_verify_parameters(mat, 'totD')
    azimuth = load_and_verify_parameters(mat, 'azimuth')
    depolarization = load_and_verify_parameters(mat, 'totP')
    return linear_retardance, diattenuation, azimuth, depolarization, mat

def load_and_verify_parameters(mat, name):
    out = mat[name]
    assert len(out) == 388
    assert len(out[0]) == 516
    return out




###########################################################################################################################
#######################################         2. Automatic square selection        ######################################
###########################################################################################################################


def select_region(shape, mask, idx, idy, square_size, grid, border = 1.5, offset = 15):
    """
    select the region that will be used to curate a pixel of the azimuth

    Parameters
    ----------
    shape : tuple
        the shape of the array
    azimuth : array of shape (388, 516)
        azimuth array
    idx, idy : int, int
        the index values of the pixel of interest
        
    Returns
    -------
    azimuth[min_x: max_x, min_y:max_y] : array
        the neighboring pixels
    """
    max_x, min_x = None, None
    max_y, min_y = None, None
    
    # special cases - borders of the image
    if idx - border * (square_size // 2 + 1) - offset < 0:
        min_x = 0 + offset
        min_x_reg = 0 + offset - border * square_size // 2
        max_x = square_size + offset
        max_x_reg = border * square_size // 2 + offset
        
    if idy - border * (square_size // 2 + 1) - offset < 0:
        min_y = 0 + offset
        min_y_reg = 0 + offset - border * square_size  // 2
        max_y = square_size + offset
        max_y_reg = border * square_size // 2 + offset
        
    if idx + border * (square_size // 2 + 1) + offset > shape[0]:
        min_x = shape[0] - square_size - offset
        min_x_reg = shape[0] - offset - border * (square_size  // 2) 
        max_x = shape[0] - offset
        max_x_reg = shape[0] - offset + border * (square_size  // 2)
        
    if idy + border * (square_size // 2 + 1) + offset > shape[1]:
        min_y = shape[1] - square_size - offset
        min_y_reg = shape[1] - border * (square_size  // 2) - offset
        max_y = shape[1] - offset
        max_y_reg = shape[1] - offset + border * (square_size  // 2)
        
    # middle of the image
    if max_x == None and min_x == None:
        min_x = idx - (square_size//2)
        min_x_reg = idx - border * (square_size//2)
        max_x = idx + (square_size//2)
        max_x_reg = idx + border * (square_size//2)
        
    if max_y == None and min_y == None:
        min_y = idy - (square_size//2)
        min_y_reg = idy - border * (square_size//2)
        max_y = idy + (square_size//2)
        max_y_reg = idy + border * (square_size//2)

    return mask[int(min_x_reg): int(max_x_reg), int(min_y_reg):int(max_y_reg)], grid[int(min_x):int(max_x), 
                                                        int(min_y):int(max_y)], [min_y, max_y, min_x, max_x]

def get_random_pixel(mask):
    randomRow = np.random.randint(mask.shape[0], size=1)
    randomCol = np.random.randint(mask.shape[1], size=1)
    return randomRow[0], randomCol[0]

def search_for_validity(mask, idx, mat = None, coordinates = None):
    positive = True
    for row in mask:
        for y in row:
            if y == idx:
                positive = False
    
    if positive and idx == 1:
        positive = sum(sum(mat['Msk'][coordinates[2]:coordinates[3], coordinates[0]:coordinates[1]])) > 0.8 * mask.shape[0]*mask.shape[1]
        if not positive:
            pass
    return positive

def update_grid(grided, coordinates):
    for idx, x in enumerate(grided):
        for idy, y in enumerate(x):
            if coordinates[0] <= idy <= coordinates[1] and coordinates[2] <= idx <= coordinates[3]:
                grided[idx, idy] = 1
    return grided

def get_square_coordinates(mask, square_size, grid, mat = None, coordinates = None):
    
    sum(sum(mat['Msk'][20:40, 40:60]))
    
    found = False
    counter = 0
    coord_ = None
    
    while not found and counter < 1000:
        randomRow, randomCol = get_random_pixel(mask)
        if mask[randomRow, randomCol] == 0:
            counter += 1
        else:
            # print(select_region(mask.shape, mask, randomRow, randomCol))
            region, grided, coordinates = select_region(mask.shape, mask, randomRow, randomCol, square_size, grid)
            
            positive_mask = search_for_validity(region, 0, mat = mat, coordinates = coordinates)
            positive_grid = search_for_validity(grided, 1, mat = mat, coordinates = coordinates)

            positive = positive_mask and positive_grid
            
            if positive:
                found = True
                coord_ = coordinates
                grid = update_grid(grid, coordinates)

            counter += 1
            
    if counter == 1000:
        pass
    
    return coordinates, grid




###########################################################################################################################
#######################################         3. Automatic image processing        ######################################
###########################################################################################################################

def square_selection(path_folder_temp, path_folder, path_folder_50x50, folder_name, wavelength, number_of_random_squares, 
                     square_size, type_rec_sq, WM, mask_matter, grid, MM_maps, mat, square_number = 1, propagated = False, 
                     pattern = '_FX_', fixation = False):
    
    continuation = True
    new_parameters = False
    new_folder = False
    propagation_list = []
    valerr = False
    
    [linear_retardance, diattenuation, azimuth, depolarization] = MM_maps
    
    # counter of number of squares (for automatic mode)
    counter = 0
    
    # get the name of the annotated image
    if WM:
        path_img_annotated = 'img_annotated_WM.png'
    else:
        path_img_annotated = 'img_annotated_GM.png'
        
    if exists(os.path.join(path_folder_50x50, path_img_annotated)):
        pass
    else:
        path_image = path_folder + '/polarimetry/' + str(wavelength) + 'nm/' + folder_name + '_' + str(wavelength) + 'nm_realsize.png'
        shutil.copyfile(path_image, os.path.join(path_folder_50x50, path_img_annotated))
        
    while continuation and counter < number_of_random_squares and not valerr:
        
        # get the path of the image and of the new output folder
        path_image = path_folder + '/polarimetry/' + str(wavelength) + 'nm/' + folder_name + '_' + str(wavelength) + 'nm_realsize.png'
        path_output, new_folder_name = get_new_output(path_folder_temp, path_folder_50x50, WM)
        with open(os.path.join(path_folder_temp, 'image_path.txt'), 'w') as file:
            file.write(path_image)
            
        # if automatic mode, get automatic square using get_square_coordinates
        try:
            coordinates_long, grid = get_square_coordinates(mask_matter, square_size, grid, mat = mat)
            [square_number, square_size_horizontal, square_size_vertical, orientation] = [1, square_size, square_size, 'None']
            all_coordinates = [coordinates_long]
        except:
            traceback.print_exc()
            valerr = True     
        
        if not valerr:
            # write the coordinates to a txt file to reuse later
            textfile = open(os.path.join(path_output, 'coordinates.txt'), 'w')
            for element in coordinates_long:
                textfile.write(str(element) + "\n")
            textfile.close()

            # get the values for the square and generate the histograms
            try:
                imnp_mask_single, imnp_mask = histogram_analysis(all_coordinates, linear_retardance, square_number, 
                                            square_size_horizontal, square_size_vertical, square_size, orientation, 
                                            diattenuation, azimuth, depolarization, coordinates_long, path_image, 
                                            path_output,
                                            path_image_ori = os.path.join('/'.join(path_output.split('/')[:-1]), 
                                            path_img_annotated), WM = WM, mat = mat)
            except RuntimeWarning:
                print(linear_retardance, diattenuation, azimuth, depolarization)
            except:
                imnp_mask = histogram_analysis(all_coordinates, linear_retardance, square_number, square_size_horizontal, 
                                               square_size_vertical, square_size, orientation, diattenuation, 
                                               azimuth, depolarization, coordinates_long, path_image, path_output,
                                               WM = WM, mat = mat)

            try:
                im = Image.fromarray(imnp_mask_single.T)
                im.save(os.path.join('/'.join(path_output.split('/')[:-1]), path_img_annotated))

            except:
                traceback.print_exc()
                im = Image.fromarray(imnp_mask.T)
                im.save(os.path.join('/'.join(path_output.split('/')[:-1]), path_img_annotated))

            new_folder = False
            new_parameters = False
            continuation = True
            propagation = True

            counter += 1

            # if yes, perform the propagation
            if propagation:
                propagated = True
                path_alignment, path_folders, all_folders = add_all_folders(path_folder, wavelength, pattern)
                img = create_and_save_mask(imnp_mask.T)
                new_folder_name = path_output.split('/')[-1]
                img.save(path_alignment + '/mask/' + path_folder.split('\\')[-1] + '_' + new_folder_name + '_selection.png')
                propagation_list.append([new_folder_name, all_folders, path_folders, wavelength, path_alignment, square_size])
                new_folder = False
                new_parameters = False
                continuation = True
                propagation = True
            
    return propagation_list


def get_new_output(path_folder_temp, path_folder_50x50, WM):
    
    if WM:
        folder_number = len([name for name in os.listdir(path_folder_50x50) if os.path.isdir(path_folder_50x50 + name) and name.replace('WM_', '').isdecimal()])
        path_output = path_folder_50x50 + 'WM_' + str(folder_number + 1)
    else:
        folder_number = len([name for name in os.listdir(path_folder_50x50) if os.path.isdir(path_folder_50x50 + name) and name.replace('GM_', '').isdecimal()])
        path_output = path_folder_50x50 + 'GM_' + str(folder_number + 1)
    
    try:
        os.mkdir(path_output)
    except:
        folder_number = len([name for name in os.listdir(path_folder_50x50) if os.path.isdir(path_folder_50x50 + name) and name.isdecimal()])
        path_output = path_folder_50x50 + str(folder_number)
    return path_output, 'GM_' + str(folder_number + 1)


def generate_summary_file_series(path_output, square_number):
    summaries_generated_fn = []
    for idx in range(square_number):
        try:
            summaries_generated_fn.append([path_output + str(idx) + '_summary.csv', idx])
        except:
            pass
    summaries_generated = []
    for file in summaries_generated_fn:
        df = pd.read_csv(file[0])
        fn = []
        for i in range(len(df)):
            fn.append(file[1])
        df['square_number'] = fn
        summaries_generated.append(df)

    if summaries_generated:
        result = pd.concat(summaries_generated)
        result = result.sort_values(['square_number', 'parameter']).set_index(['square_number', 'square_size', 'parameter'])
        result.to_csv(path_output + 'summaries.csv')
        result.to_excel(path_output + 'summaries.xlsx')
    
def histogram_analysis(all_coordinates, linear_retardance, square_number, square_size_horizontal, 
                       square_size_vertical, square_size, orientation, diattenuation, azimuth, depolarization, 
                       coordinates_long, path_image, path_output, path_image_ori = None, WM = True,
                       mat = None):
    data = []
    for idx, coordinates in enumerate(all_coordinates):
        data.append(analyze_and_get_histograms(linear_retardance, diattenuation, azimuth, depolarization, mat, 
                                               coordinates)[1:])
        retardance, diattenua, azi, depol = data[-1]
        fig = generate_histogram(retardance, diattenua, azi, depol, path_output + '/', idx)
        save_parameters(retardance, diattenua, azi, depol, path_output + '/', square_size, idx)
        
    imnp_mask = generate_pixel_image(coordinates_long, path_image, path_output + '/', [square_number, [square_size_horizontal, square_size_vertical], orientation], WM = True)
    
    try:
        imnp_mask_single = generate_pixel_image(coordinates_long, path_image_ori, path_output + '/', 
                                [square_number, [square_size_horizontal, square_size_vertical], orientation], combined = True, WM = WM)
    except:
        traceback.print_exc()
        pass
    
    generate_summary_file_series(path_output + '/', square_number)
    
    try:
        return imnp_mask_single, imnp_mask
    except:
        return imnp_mask
    
    
def analyze_and_get_histograms(linear_retardance, diattenuation, azimuth, depolarization, mat, coordinates, 
                               imnp = None, angle = 0):
    path_folder_temp = './temp'
    
    with open(path_folder_temp + '/histogram_parameters.pickle', 'rb') as handle:
        parameters = pickle.load(handle)
    
    if type(imnp) == np.ndarray:
        masked = True
    else:
        masked = False
        
    if masked:
        coordinates = imnp
    else:
        pass
    
    retardance = get_area_of_interest(coordinates, linear_retardance, parameters['retardance'], masked, mat = mat)
    diattenua = get_area_of_interest(coordinates, diattenuation, parameters['diattenuation'], masked, mat = mat)
    if type(angle) == float:
        azi = get_area_of_interest(coordinates, azimuth, parameters['azimuth'], masked, angle, mat = mat)
    else:
        azi = get_area_of_interest(coordinates, azimuth, parameters['azimuth'], masked, angle, mat = mat)
    depol = get_area_of_interest(coordinates, depolarization, parameters['depolarization'], masked, mat = mat)
    return coordinates, retardance, diattenua, azi, depol

def get_area_of_interest(params, matrix, param, masked, angle = 0, mat = None):
    
    if masked:
        listed = []
        for idx_x, x in enumerate(params):
            for idx_y, y in enumerate(x):
                if y != 0 and mat['Msk'][idx_x, idx_y]:
                    # listed.append((matrix[idx_x][idx_y] - angle)%180)
                    listed.append(matrix[idx_x][idx_y])
    else:
        [y_min, y_max, x_min, x_max] = params
        listed = []
        for idx_x, x in enumerate(mat['Msk']):
            for idx_y, y in enumerate(x):
                if x_min <= idx_x < x_max and y_min <= idx_y < y_max and y:
                    # listed.append((matrix[idx_x][idx_y] - angle)%180)
                    listed.append(matrix[idx_x][idx_y])

    mean = np.mean(listed)
    stdev = np.std(listed)
    median = np.median(listed)
    bins = np.linspace(param['borders'][0], param['borders'][1], num = param['num_bins'])
    data = plt.hist(listed, bins = bins)
    arr = data[0]
    max_idx = np.where(arr == np.amax(arr))[0][0]
    maximum = data[1][max_idx]
    return mean, stdev, maximum, listed, median


def get_params_summary(params, name):
    return [name, params[0], params[1], params[2], params[4]]

def save_parameters(retardance, diattenua, azi, depol, path_output, square_size, idx = None):
    retardance_params = get_params_summary(retardance, 'retardance')
    diatten_params = get_params_summary(diattenua, 'diattenuation')
    azi_params = get_params_summary(azi, 'azimuth')
    depol_params = get_params_summary(depol, 'depolarization')
    params = [retardance_params, diatten_params, azi_params, depol_params]
    df = pd.DataFrame(params, columns = ['parameter', 'mean', 'stdev', 'max', 'median'])
    df = df.set_index('parameter')
    fn = []
    for i in range(len(df)):
        fn.append(str(None))
    df['square_size'] = fn
    if idx == None:
        df.to_csv(path_output + 'summary.csv')
        df.to_excel(path_output + 'summary.xlsx')
    else:
        df.to_csv(path_output + str(idx) + '_summary.csv')
        df.to_excel(path_output + str(idx) + '_summary.xlsx')
        
    return df

def generate_pixel_image(coordinates, path_image, path_save, params = None, mask = None, path_save_align = None, combined = False,
                        save = True, WM = True):
    
    if WM:
        val = 0
    else:
        val = 255
        
    if type(mask) == np.ndarray:
        masked = True
    else:
        [x_min, x_max, y_min, y_max] = coordinates
        masked = False
    
    im = Image.open(path_image)
    imnp = np.array(im)
    imnp_mask = np.array(im)
        
    if masked:
        for idx_x, x in enumerate(imnp):
            
            min_idx_x = max(idx_x - 1, 0)
            max_idx_x = min(idx_x + 1, len(mask) - 1)
            
            for idx_y, y in enumerate(x):
                if mask[idx_x][idx_y] == 0 :
                    pass
                
                else:
                    
                    min_idx_y = max(idx_y - 1, 0)
                    max_idx_y = min(idx_y + 1, len(mask[0]) - 1)
                    
                    mask_idx_y = mask[:, idx_y]
                    
                    if mask_idx_y[min_idx_x] == 0:
                        imnp[idx_x][idx_y] = val
                    elif mask_idx_y[max_idx_x] == 0:
                        imnp[idx_x][idx_y] = val
                        
                    else:
                        mask_idx_x = mask[idx_x, :]
                        if mask_idx_x[min_idx_y] == 0:
                            imnp[idx_x][idx_y] = val
                        elif mask_idx_x[max_idx_y] == 0:
                            imnp[idx_x][idx_y] = val

    else:
        for y, row in enumerate(imnp):
            if y < y_min - 1 or y > y_max + 1:
                pass
            else:
                for x, column in enumerate(row):
                    if x < x_min - 1 or x > x_max + 1:
                        pass
                    else:
                        if x_min - 1 <= x <= x_min + 1 or x_max - 1 <= x <= x_max + 1:
                            imnp[y][x] = val
                        if y_min - 1 <= y <= y_min + 1 or y_max - 1 <= y <= y_max + 1:
                            imnp[y][x] = val
                        if x_min <= x <= x_max and y_min <= y <= y_max:
                            imnp_mask[y][x] = val

    if params:
        square_number, square_size, orientation = params
        all_coordinates = get_each_image_coordinates(coordinates, square_number, square_size, orientation)
        for coordinates in all_coordinates:
            [x_min, x_max, y_min, y_max] = coordinates
            for y, row in enumerate(imnp):
                if y < y_min or y > y_max:
                    pass
                else:
                    for x, column in enumerate(row):
                        if x < x_min or x > x_max:
                            pass
                        else:
                            if x_min <= x <= x_min or x_max <= x <= x_max:
                                imnp[y][x] = val
                            if y_min <= y <= y_min or y_max <= y <= y_max:
                                imnp[y][x] = val
                            else:
                                pass
        
    if type(path_save_align) == str and save:  
        Image.fromarray(imnp).save(path_save_align + path_image.split('/')[-1])
    if not combined:
        Image.fromarray(imnp).save(path_save + 'selection.png')
    
    return imnp_mask.T   

def get_each_image_coordinates(coordinates, square_number, square_size, orientation):
    if type(square_size) == int:
        square_size_horizontal = square_size
        square_size_vertical = square_size
    else:
        square_size_horizontal, square_size_vertical = square_size
        
    coordinates_all = []
    x_min, x_max, y_min, y_max = None, None, None, None
    if orientation == 'horizontal':
        for i in range(square_number):
            x_min = coordinates[0] + i*square_size_horizontal
            x_max = coordinates[0] + (i+1)*square_size_horizontal
            y_min = coordinates[2] 
            y_max = coordinates[3]
            coordinates_all.append([x_min, x_max, y_min, y_max])
    else:
        for i in range(square_number):
            x_min = coordinates[0]
            x_max = coordinates[1]
            y_min = coordinates[2] + i*square_size_vertical
            y_max = coordinates[2] + (i+1)*square_size_vertical
            coordinates_all.append([x_min, x_max, y_min, y_max])
    return coordinates_all


def generate_histogram(retardance, diattenua, azi, depol, path_folder, idx = None):
    fig, ((ax0, ax1), (ax2, ax3)) = plt.subplots(nrows=2, ncols=2, figsize=(15, 10))

    n_bins = 100

    sns.distplot(retardance[3], n_bins, ax= ax0)
    ax0.set_title('Linear retardance', fontsize = 16)
    ax0.get_yaxis().set_visible(False)
    sns.distplot(diattenua[3], n_bins, ax= ax1, color = 'salmon')
    ax1.set_title('Diattenuation', fontsize = 16)
    ax1.get_yaxis().set_visible(False)
    sns.distplot(azi[3], n_bins, ax= ax2, color = 'palegreen')
    ax2.set_title('Azimuth', fontsize = 16)
    ax2.set_xlim([0, 180])
    ax2.get_yaxis().set_visible(False)
    sns.distplot(depol[3], n_bins, ax= ax3, color = 'bisque')
    ax3.set_title('Depolarization', fontsize = 16)
    ax3.get_yaxis().set_visible(False)

    add_text(ax0, retardance)
    add_text(ax1, diattenua)
    add_text(ax2, azi)
    add_text(ax3, depol)

    fig.tight_layout()
    if idx == None:
        plt.savefig(path_folder + 'histogram.png')
        plt.savefig(path_folder + 'histogram.pdf')
        plt.close(fig)
    else:
        plt.savefig(path_folder + str(idx) + '_histogram.png')
        plt.savefig(path_folder + str(idx) + '_histogram.pdf')
        plt.close(fig)
        
def add_text(ax, res):
    ax.text(0.7, 0.8, '$\mu$=' + "%.2f" % res[0] + '$,\ median$=' + "%.2f" % res[4] +'\n$\sigma=$' + "%.2f" % res[1] + '$,\ \max = $' + "%.2f" % res[2], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
    
    
###########################################################################################################################
####################################         4. Automatic propagation of the ROIs       ###################################
###########################################################################################################################


def add_all_folders(path_folder, wavelength, pattern = '_FX_', automatic = True):
    
    end_pattern = path_folder.split('P-T0_FX_')[-1]
    all_folders = find_other_measurements(path_folder, pattern, automatic, end_pattern = end_pattern)
    
    folder_name = path_folder.split('/')[-1]
    path_alignment = 'alignment/to_align/'
    aligned = False
    path_aligned = None
    for folder in os.listdir(path_alignment):
        if folder_name in folder:
            aligned = True
            path_aligned = os.path.join(path_alignment, folder)
    
    path_folders = '/'.join(path_folder.split('\\')[:-1])
    
    if aligned:
        path_alignment = path_aligned
    else:
        now = datetime.now()
        dt_string = path_folder.split('\\')[-1] +'__' + now.strftime("%d/%m/%Y %H:%M:%S").replace(' ', '_').replace('/', '_').replace(':', '_')
        path_alignment = './alignment/to_align/' + dt_string
        os.mkdir(path_alignment)
        os.mkdir(path_alignment + '/mask')

        target = path_alignment
        for folder in all_folders:
            for file in os.listdir(path_folders + '/' + folder + '/polarimetry/' + str(wavelength) + 'nm/'):
                if file.endswith('realsize.png'):
                    shutil.copy(path_folders + '/' + folder + '/polarimetry/' + str(wavelength) + 'nm/' + file, target)


        for file in os.listdir(path_alignment):
            old_name = path_alignment + '/' + file
            new_name = path_alignment + '/' + file.split('.')[0] + '_ref_align.png'
            if path_folder.split('/')[-1] in file:
                os.rename(old_name, new_name)

    return path_alignment, path_folders, all_folders


def find_other_measurements(path_folder, pattern = '_FX_', automatic = False, end_pattern = ''):
    CUSA = 'CUSA' in path_folder
        
    if CUSA:
        measurement = path_folder.split(pattern)[0].split('_')[-1]
    else:
        measurement = path_folder.split(pattern)[1]

    all_folders = []
    for folder in os.listdir('/'.join(path_folder.split('\\')[:-1])):
        if 'P-T0' and measurement in folder and end_pattern in folder:
            all_folders.append(folder)

    return all_folders


def create_and_save_mask(imnp_mask):
    imnp_masked = []
    for line in imnp_mask:
        line_mask = []
        for pixel in line:
            if pixel[-1] > 0:
                line_mask.append(0)
            else:
                line_mask.append(1)
        imnp_masked.append(line_mask)
        
    img = Image.new('1', (len(imnp_masked), len(imnp_masked[0])))
    pixels = img.load()

    for i in range(img.size[0]):
        for j in range(img.size[1]):
            pixels[i, j] = imnp_masked[i][j]
    img = img.transpose(Image.FLIP_LEFT_RIGHT).rotate(90, expand=True)
    return img




###########################################################################################################################
##########################################         5. Automatic propagation        ########################################
###########################################################################################################################


def generate_combined_mask(propagation_list):
    masks = []
    imgs = {}
    for file in os.listdir('alignment/to_align/'):
        for img_path in os.listdir(os.path.join('alignment/to_align/', file, 'mask')):
            im = plt.imread(os.path.join('alignment/to_align/', file, 'mask', img_path))
            imgs[int(img_path.split('WM_')[-1].split('GM_')[-1].split('_')[0])] = im
            
    base = np.zeros(imgs[1].shape)
    for val, img in imgs.items():
        assert val < 255
        for idx, x in enumerate(img):
            for idy, y in enumerate(x):
                if y != 0:
                    base[idx, idy] = val

    for img_path in os.listdir(os.path.join('alignment/to_align/', file, 'mask')):
        os.remove(os.path.join('alignment/to_align/', file, 'mask', img_path))

    mask = Image.fromarray(base)
    mask = mask.convert("L")
    mask.save(os.path.join('alignment/to_align/', file, 'mask', 'mask.png'))
    
    for folder in os.listdir('alignment/to_align/'):
        path_folder = os.path.join('alignment/to_align/', folder)
        if file in folder:
            path_folder_propagation = path_folder
            for fname in os.listdir(path_folder):
                if 'P-T0_' in fname:
                    os.rename(os.path.join(path_folder, fname), 
                              os.path.join(path_folder, fname.replace('realsize', 'realsize_ref_align')))
                    
        else:
            pass
            shutil.rmtree(path_folder)
            
    for prop in propagation_list:
        prop[-2] = path_folder_propagation
    return propagation_list

def do_alignment():
    path_alignment_batch = os.getcwd() + '\\alignment\\to_align'
    with open('RegistrationElastix/temp/path_alignment_batch.txt', 'w') as f:
        f.write(path_alignment_batch)
    f.close()

    FixPattern = '_ref_align'
    with open('RegistrationElastix/temp/FixPattern.txt', 'w') as f:
        f.write(FixPattern)
    f.close()

    Tag = 'AffineElastic'
    with open('RegistrationElastix/temp/Tag.txt', 'w') as f:
        f.write(Tag)
    f.close()

    eng = matlab.engine.start_matlab()
    path = os.getcwd() + '/RegistrationElastix/RegistrationScripts'
    eng.cd(path, nargout=0)
    s = eng.genpath('0_NIfTI_IO')
    eng.addpath(s, nargout=0)
    eng.python_call(nargout=0)
    
def move_computed_folders():
    
    folder_names = []
    log_name = []
    for fname in os.listdir('alignment/to_align/'):
        if not os.path.isfile('alignment/to_align/' + fname):
            folder_names.append('alignment/to_align/' + fname)
        else:
            log_name.append('alignment/to_align/' + fname)

    output_folders = {}
    for folder_name in folder_names:
        output_folder = 'alignment/aligned/'
        if folder_name.split('/')[-1] in os.listdir(output_folder):
            output_folder = output_folder + folder_name.split('/')[-1]
            output_folders[folder_name] = output_folder
        else:
            os.mkdir(output_folder + folder_name.split('/')[-1])
            output_folder = output_folder + folder_name.split('/')[-1]
            output_folders[folder_name] = output_folder
        
    for source, target in output_folders.items():
        folder_name_ = target.split('/')[-1]
        shutil.move(source, target)
        os.rename(target + '/' + folder_name_, target + '/results')
    
    output_folder = 'alignment/aligned/logbooks/'
    for log in log_name:
        shutil.move(log, output_folder)
    
    return output_folders

def create_shortcut(path_images, path_folder):
    path = os.path.join(path_folder, 'aligned_images.lnk')
    target = path_images
    icon = path_images
    shell = win32com.client.Dispatch("WScript.Shell")
    shortcut = shell.CreateShortCut(path)
    shortcut.Targetpath = target
    shortcut.IconLocation = icon
    shortcut.save()
    
    
    
    
###########################################################################################################################
#######################################         6. Propagated data collection        ######################################
###########################################################################################################################

def check_outliers_propagation(all_folders, path_alignment, new_folder_name, mask_matter_after, mask_matter_after_opposite, 
                               elastic = True, check_outliers_bool = False):
    
    elastic = True

    for directory in os.listdir('alignment/aligned/'):
        if path_alignment.split('/')[-1] in directory:
            path_aligned = 'alignment/aligned/' + directory + '/results'
            path_aligned_root = 'alignment/aligned/' + directory
    assert path_aligned != None and path_aligned_root != None

    to_remove = []
    
    base_folder = path_alignment.split('/')[-1].split('__')[0]
    
    for folder in all_folders:

        path_image = None
        img_of_interest = []
        for img in os.listdir(path_aligned + '/invReg'):
            if folder in img and '_PrpgTo_' in img and 'AffineElastic_TransformParameters' in img and '.png' in img:
                img_of_interest.append(img)
        assert len(img_of_interest) == 2

        for img in img_of_interest:
            if img.endswith('1.png') and elastic:
                path_image = path_aligned + '/invReg/' + img
            elif img.endswith('0.png') and not elastic:
                path_image = path_aligned + '/invReg/' + img
            
        mask = load_image(path_image, new_folder_name)
        
        if check_outliers(mask, mask_matter_after, mask_matter_after_opposite):
            to_remove.append(new_folder_name)
        
    return to_remove


def check_outliers(mask, mask_matter_after, mask_matter_after_opposite):
    corrupt = False
    
    opposite = 0
    same = 0
    total = 0
    
    for idx, x in enumerate(mask):
        for idy, y in enumerate(x):
            if y > 0:
                total += 1
                if idx == 0 or idy == 0:
                    corrupt = True
                elif idx == mask.shape[0] - 1 or idx == mask.shape[1] - 1:
                    corrupt = True
                elif mask_matter_after[idx, idy] == 1:
                    same += 1
                elif mask_matter_after_opposite[idx, idy] == 1:
                    opposite += 1
    
    try:
        if same/total > 0.50 or (same/total > 0.30 and same/opposite == 0):
            pass
        else:
            corrupt = True
    except:
        pass
    return corrupt


def collect_data_propagated(WM, new_folder_names, new_dates, old_folder_name, old_date, path_folder, propagation_list, output_folders):
    if WM:
        type_ = 'WM'
    else:
        type_ = 'GM'

    mask_matter_afters = []
    mask_matter_after_opposites = []
    
    for idx, (new_name, new_date) in enumerate(zip(new_folder_names,new_dates)):
        _ = get_masks(path_folder.replace(old_folder_name, new_name).replace(old_date, new_date))
        path_annotation = os.path.join(path_folder.replace(old_folder_name, new_name).replace(old_date, new_date),
                                           'annotation')
        WM_mask = plt.imread(os.path.join(path_annotation, 'WM_merged.png'))
        GM_mask = plt.imread(os.path.join(path_annotation, 'GM_merged.png'))
        if WM:
            mask_matter_afters.append(WM_mask)
            mask_matter_after_opposites.append(GM_mask)
        else:
            mask_matter_afters.append(GM_mask)
            mask_matter_after_opposites.append(WM_mask)

    remove = []

    for element in propagation_list:
        new_folder_name, all_folders, path_folders, wavelength, path_alignment, square_size = element

        for idx, (all_folder, mask_matter_after, mask_matter_after_opposite, new_name) in enumerate(zip(all_folders[1:],
                                                            mask_matter_afters, mask_matter_after_opposites, new_folder_names)):

            to_remove = check_outliers_propagation([all_folder], path_alignment, new_folder_name, mask_matter_after, 
                                                       mask_matter_after_opposite, elastic = True, check_outliers_bool = True)

            if len(to_remove) > 0:
                remove.append([to_remove[0], new_name])

            data, dfs, aligned_images_path = propagate_measurements(new_folder_name, [all_folder], 
                                                    path_folders, wavelength, output_folders[path_alignment], 
                                                    square_size, mask_matter_after, mask_matter_after_opposite,
                                                    check_outliers_bool = False, create_dir = True)

            with open(os.path.join(path_folders, all_folder, 'polarimetry', '550nm', '50x50_images', 
                                new_folder_name + '_align', 'data_raw' + '.pickle'), 'wb') as handle:
                pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    with open(os.path.join(path_folders, all_folders[0], 'polarimetry', '550nm', '50x50_images', 
                        'to_rm_' + type_ + '.pickle'), 'wb') as handle:
        pickle.dump(remove, handle, protocol=pickle.HIGHEST_PROTOCOL)

    generate_summary_file(propagation_list)
    
    
def propagate_measurements(new_folder_name, all_folders, path_folders, wavelength, path_alignment, 
                           square_size, mask_matter_after, mask_matter_after_opposite, check_outliers_bool = False, 
                           create_dir = False):
    output_directory = get_output_directory_name(new_folder_name, all_folders, path_folders, wavelength)
    
    if create_dir:
        create_output_dir(output_directory, all_folders, path_folders, wavelength)

    if check_outliers_bool:
        to_remove = check_outliers_propagation(output_directory, all_folders, path_folders, wavelength, path_alignment, 
                                               square_size, new_folder_name, mask_matter_after, 
                                               mask_matter_after_opposite, elastic = True, 
                                               check_outliers_bool = False, create_dir = False)
        return to_remove
    else:
        data, dfs, aligned_images_path = get_data_propagation(output_directory, all_folders, 
                                                                         path_folders, wavelength, 
                                                                         path_alignment, square_size, 
                                                                         new_folder_name, 
                                                                         mask_matter_after, 
                                                                         mask_matter_after_opposite,
                                                                         check_outliers_bool = check_outliers_bool,
                                                                         create_dir = create_dir)
        # if create_dir:
            # save_dfs(dfs, all_folders, path_folders, wavelength, output_directory, aligned_images_path)
        return data, dfs, aligned_images_path
    
    
def create_output_dir(output_directory, all_folders, path_folders, wavelength):
    for folder in all_folders:
        path_folder_50x50 = path_folders + '/' + folder + '/polarimetry/' + str(wavelength) + 'nm/50x50_images/'
        os.mkdir(path_folder_50x50 + output_directory)
    
    
def get_output_directory_name(new_folder_name, all_folders, path_folders, wavelength):
    new_folder_name = new_folder_name + '_align'
    new_folder_free = True
    if not new_folder_name:
        new_folder_free = False
        
    number = 0
    for folder in all_folders:
        path_folder_50x50 = path_folders + '/' + folder + '/polarimetry/' + str(wavelength) + 'nm/50x50_images/'
        folders_ = os.listdir(path_folder_50x50)
        if new_folder_name in folders_:
            new_folder_free = False
        folder_number = len(folders_)
        if folder_number > number:
            number = folder_number
            
    if new_folder_free:
        return new_folder_name
    else:
        return str(number) + '_align'
    
    
def get_data_propagation(output_directory, all_folders, path_folders, wavelength, path_alignment, square_size, 
                         new_folder_name, mask_matter_after, mask_matter_after_opposite, elastic = True, 
                         check_outliers_bool = False, create_dir = False):
    
    # get the path to the folder containing the aligned images
    path_aligned = None
    path_aligned_root = None
    elastic = True
    for directory in os.listdir('alignment/aligned/'):
        if path_alignment.split('/')[-1] in directory:
            path_aligned = 'alignment/aligned/' + directory + '/results'
            path_aligned_root = 'alignment/aligned/' + directory
    assert path_aligned != None and path_aligned_root != None
    try:
        aligned_images_path = path_aligned_root + '/aligned_images'
        os.mkdir(aligned_images_path)
    except FileExistsError:
        pass

    
    data = []
    dfs = []
    
    base_folder = path_alignment.split('/')[-1].split('__')[0]
    
    for folder in all_folders:

        path_folder = path_folders + '/' + folder
        path_output = path_folders + '/' + folder + '/polarimetry/' + str(wavelength) + 'nm/50x50_images/' + output_directory
        path_original_image = path_folder + '/polarimetry/' + str(wavelength) + 'nm/' + path_folder.split('/')[-1] + '_' + str(wavelength) + 'nm_realsize.png'
        
        if folder == base_folder:
            dfs.append([pd.read_csv(os.path.join(path_folder, 'polarimetry', str(wavelength) + 'nm', '50x50_images',
                               new_folder_name, '0_summary.csv')).set_index('parameter'), base_folder])
            selected_image = path_output.split('_align')[0] + '/selection.png'
            output_path = path_aligned_root + '/aligned_images/selection' + folder + '.png'
            shutil.copy(selected_image, output_path)
            
        else:
            
            # get the paths to the images
            path_image = None
            if folder in path_alignment:
                assert len(os.listdir(path_aligned + '/mask')) == 1
                path_image = path_aligned + '/mask/' + os.listdir(path_aligned + '/mask')[0]
            else:
                img_of_interest = []
                for img in os.listdir(path_aligned + '/invReg'):
                    if folder in img and '_PrpgTo_' in img and 'AffineElastic_TransformParameters' in img and '.png' in img:
                        img_of_interest.append(img)
                assert len(img_of_interest) == 2

                for img in img_of_interest:
                    if img.endswith('1.png') and elastic:
                        path_image = path_aligned + '/invReg/' + img
                    elif img.endswith('0.png') and not elastic:
                        path_image = path_aligned + '/invReg/' + img

            # load the propagated mask image
            angle = 0
            mask = load_image(path_image, new_folder_name)
            
            # load the data (polarimetric parameters)
            linear_retardance, diattenuation, azimuth, depolarization, mat = load_data_mm(path_folder, wavelength)
            data.append(analyze_and_get_histograms(linear_retardance, diattenuation, azimuth, depolarization, mat, None, 
                                                   imnp = mask, angle = angle)[1:])
            retardance, diattenua, azi, depol = data[-1]

            if create_dir:
                # fig = generate_histogram(retardance, diattenua, azi, depol, path_output + '/', None)
                dfs.append([save_parameters(retardance, diattenua, azi, depol, path_output + '/', square_size), folder])

            path_original_image = path_folder + '/polarimetry/' + str(wavelength) + 'nm/' + path_folder.split('/')[-1] + '_' + str(wavelength) + 'nm_realsize.png'

            if create_dir:
                try:
                    imnp_mask = generate_pixel_image(None, path_original_image, path_output + '/', mask = mask, 
                                                     path_save_align = path_aligned_root + '/aligned_images/', save = False)
                except:
                    pass

            if create_dir:
                selected_image = path_output + '/selection.png'
                output_path = path_aligned_root + '/aligned_images/selection' + folder + '.png'
                shutil.copy(selected_image, output_path)
        
    return data, dfs, aligned_images_path


def generate_summary_file(propagation_list):
    
    new_folder_names = []
    all_folders = None
    path_folders = None
    wavelength = None
    
    for element in propagation_list:
        new_folder_name, all_folders, path_folders, wavelength, _, _ = element
        new_folder_names.append(new_folder_name)
        
    for new_folder_name in new_folder_names:
        
        summaries = []
        for folder in all_folders:
            path_50x50_folder = get_path_50x50_folder(path_folders, folder, new_folder_name, wavelength)
            df = pd.read_csv(path_50x50_folder)
            df['fname'] = [folder] * len(df)
            summaries.append(df)
            
        result = pd.concat(summaries)
        result = result.sort_values(['fname', 'parameter']).set_index(['fname', 'parameter'])
        
        for folder in all_folders:
            path_50x50_folder = get_path_50x50_folder(path_folders, folder, new_folder_name, wavelength, aligned = True)
            result.to_csv(os.path.join(path_50x50_folder, 'summary_aligned.csv'))
            result.to_excel(os.path.join(path_50x50_folder, 'summary_aligned.xlsx'))
            
            
def get_path_50x50_folder(path_folders, folder, new_folder_name, wavelength, aligned = False):
    path_50x50 = os.path.join(path_folders, folder, 'polarimetry', str(wavelength) + 'nm', '50x50_images')
    if os.path.exists(os.path.join(path_50x50, new_folder_name)):
        if aligned:
            path_50x50_folder = os.path.join(path_50x50, new_folder_name)
        else:
            path_50x50_folder = os.path.join(path_50x50, new_folder_name, '0_summary.csv')
    else:
        if aligned:
            path_50x50_folder = os.path.join(path_50x50, new_folder_name + '_align')
        else:
            path_50x50_folder = os.path.join(path_50x50, new_folder_name + '_align', 'summary.csv')
        assert os.path.exists(path_50x50_folder)
    return path_50x50_folder


def load_image(path_image, new_folder_name = 'WM_1'):
    val = int(new_folder_name.split('_')[-1])
    
    im = Image.open(path_image)
    imnp = np.array(im)
    
    for idx, x in enumerate(imnp):
        for idy, y in enumerate(x):
            if y != val :
                imnp[idx, idy] = 0
    return imnp




###########################################################################################################################
###########################################         7. Show the histograms        #########################################
###########################################################################################################################

def parameters_histograms(MM, xs, ys, path_res):
    """
    generate the histogram for the four parameters

    Parameters
    ----------
    MuellerMatrices : dict
        the dictionnary containing the computed Mueller Matrices
    folder : path
        the name of the current processed folder
    parameters_map : dict
        the map linking the keys in the MM and the name of the parameter
    logbook : file
        the logbook file
    max_ : bool
        boolean indicating wether or not the max_ should be printed
    """    
    parameters_map = {'linR': ['Linear retardance (°)', False, False, False, False],
    'totP': ['Depolarization', False, True, False, False],
     'totD': ['Diattenuation', False, False, True, False],
     'M11': ['Intensity', False, False, False, True],
     'azimuth': ['Azimuth of optical axis (°)', True, False, False, False]}
    
    try:
        parameters_map.pop('M11')
    except:
        pass
    
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(15,11))
    
    data_all = {}
    
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
        
        data_MM = []
        for idx, x in enumerate(MM[key]):
            for idy, y in enumerate(x):
                if MM['Msk'][idx, idy] and xs[0] < idx < xs[1] and ys[0] < idy < ys[1]:
                    data_MM.append(y)
                    
        data_all[param[0]] = data_MM
        
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

def subtract_angle(targetA, sourceA):
    a = targetA - sourceA
    return abs((a + 90) % 180 - 90)

def average_angles(angles):
    """Average (mean) of angles

    Return the average of an input sequence of angles. The result is between
    ``0`` and ``2 * math.pi``.
    If the average is not defined (e.g. ``average_angles([0, math.pi]))``,
    a ``ValueError`` is raised.
    """

    x = sum(math.cos(a * 2*math.pi/180) for a in angles)
    y = sum(math.sin(a * 2*math.pi/180) for a in angles)

    if x == 0 and y == 0:
        raise ValueError(
            "The angle average of the inputs is undefined: %r" % angles)

    # To get outputs from -pi to +pi, delete everything but math.atan2() here.
    return 180 * math.fmod(math.atan2(y, x) + 2 * math.pi, 2 * math.pi) / (2*math.pi)