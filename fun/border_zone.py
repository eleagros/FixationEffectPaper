import os
import pickle
import copy
from collections import defaultdict
import numpy as np
from PIL import Image
import cv2
from matplotlib import pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from tqdm.notebook import trange, tqdm
from scipy.optimize import curve_fit
from scipy import stats
from kneed import KneeLocator


###########################################################################################################################
##########################         1. Load variables and MMs and create the combined masks        #########################
###########################################################################################################################

def load_variables_and_paths(path_fixation_folder, match_sequence = '', wavelength = '550nm', to_keep = ['M11', 'Msk', 'totD', 'linR', 'azimuth', 'totP']):
    """
    load_variables_and_paths is the function called to load the images and MMs for the fixation measurements

    Parameters
    ----------
    path_fixation_folder : str
        the path to the folder containing the fixation measurements
    match_sequence : str
        the sequence that is used to identify the folders of interest (default '')
    wavelength : str
        the wavelength of interest (default '550nm')
    to_keep : list of str
        the keys of the MMs to keep (default ['M11', 'Msk', 'totD', 'linR', 'azimuth', 'totP'])
        
    Returns
    -------
    paths : dict of str
        the paths to the annotation, MM, and image for each measurement folder
    variables : dict of var
        the annotated image, MM and grayscale image for each measurement folder
    """
    folders = os.listdir(path_fixation_folder)
    folder_of_interests = []
    for folder in folders:
        if match_sequence in folder:
            folder_of_interests.append(os.path.join(path_fixation_folder, folder))
            
    paths = {}
    variables = {}

    for folder_of_interest in tqdm(folder_of_interests):
        path_annotation = os.path.join(folder_of_interest, 'annotation', 'merged.jpeg')

        path_polarimetry = os.path.join(folder_of_interest, 'polarimetry', wavelength)
        path_MM = os.path.join(path_polarimetry, 'MM.npz')

        for f in os.listdir(path_polarimetry):
            if '_realsize.png' in f:
                path_image = os.path.join(path_polarimetry, f)

        paths[folder_of_interest] = {'annotation': path_annotation, 'MM': path_MM, 'grayscale': path_image}

        im = np.asarray(Image.open(path_annotation)).astype("uint8")
        im_annotated = np.zeros(im.shape)
        for idx, x in enumerate(im):
            for idy, y in enumerate(x):
                if y < 80:
                    im_annotated[idx, idy] = 0
                elif 80 <= y <= 176:
                    im_annotated[idx, idy] = 128
                else:
                    im_annotated[idx, idy] = 255

        MM = np.load(path_MM)
        
        MM_light = {}
        for key, val in MM.items():
            if key in to_keep:
                MM_light[key] = val
            else:
                pass

        im_gs = np.asarray(Image.open(path_image))
        variables[folder_of_interest] = {'annotation': im_annotated, 'MM': MM_light, 'grayscale': im_gs}
        
    return paths, variables


def create_the_masks(path_fixation_folder, match_sequence):
    """
    master function to create the combined match for all the folders located in path_fixation_folder

    Parameters
    ----------
    path_fixation_folder : str
        the path to the measurement folder
    match_sequence : str
        the sequence to match to be considered a measurement folder
    """
    folders = os.listdir(path_fixation_folder)

    folder_of_interests = []
    for folder in folders:
        if match_sequence in folder:
            folder_of_interests.append(os.path.join(path_fixation_folder, folder))

    for folder_of_interest in tqdm(folder_of_interests):
        _ = get_masks(folder_of_interest)
        

def combine_masks(masks):
    """
    combine previsously manually drawn masks

    Parameters
    ----------
    masks : list of array of shape (388,516)
        the manually drawn masks

    Returns
    -------
    base : array of shape (388,516)
        the combined mask
    """
    if masks:
        # use the first mask as a base
        base = masks[0]
        
        # for each of the mask, search for pixels equal to 255 and add them as positive values to the base mask
        for id_, mask in enumerate(masks[1:]):
            for idx, x in enumerate(mask):
                for idy, y in enumerate(x):
                    if base[idx, idy] == 255 or y == 255:
                        base[idx, idy] = 255
    
    # if no mask is found, everything is set to 0
    else:
        base = np.zeros((388, 516))
    return base


def get_masks(path, bg = True):
    """
    obtain the masks (white matter, grey matter and background) by combining previsously manually drawn masks

    Parameters
    ----------
    path : str
        the path to the folder containing the annotations considered

    Returns
    -------
    BG_merged : array of shape (388,516)
        the background mask
    GM_merged : array of shape (388,516)
        the grey matter mask
    WM_merged : array of shape (388,516)
        the white matter mask
    all_merged : array of shape (388,516)
        the three masks combined (one color = one class)
    """
    BG = []
    WM = []
    WM_found = None
    BG_found = None
    path = os.path.join(path, 'annotation')
    
    # add the masks to the different lists
    for file in os.listdir(path):
        im = Image.open(os.path.join(path, file))
        imarray = np.array(im)
        
        if 'BG' in file and not 'merged' in file:
            BG.append(imarray)
        elif 'WM' in file and not 'merged' in file:
            WM.append(imarray)
        elif 'merged' in file or 'mask-viz' in file or 'mask_total' in file:
            pass
        else:
            raise(NotImplementedError)
    
    
    # combine the WM and BG masks
    WM = combine_masks(WM)
    BG = combine_masks(BG)
    GM = np.zeros(WM.shape)
    
    # return the merged masks
    return merge_masks(BG, WM, GM, path, bg)


def combine_masks(masks):
    """
    combine previsously manually drawn masks

    Parameters
    ----------
    masks : list of array of shape (388,516)
        the manually drawn masks

    Returns
    -------
    base : array of shape (388,516)
        the combined mask
    """
    if masks:
        # use the first mask as a base
        base = masks[0]
        
        # for each of the mask, search for pixels equal to 255 and add them as positive values to the base mask
        for id_, mask in enumerate(masks[1:]):
            for idx, x in enumerate(mask):
                for idy, y in enumerate(x):
                    if base[idx, idy] == 255 or y == 255:
                        base[idx, idy] = 255
    
    # if no mask is found, everything is set to 0
    else:
        base = np.zeros((388, 516))
    return base


def merge_masks(BG, WM, GM, path, bg):
    """
    merge masks is used to merge the previsously combined manually drawn masks

    Parameters
    ----------
    BG : array of shape (388,516)
        the combined background mask
    GM : array of shape (388,516)
        the combined grey matter mask
    WM : array of shape (388,516)
        the combined white matter mask
    path : str
        the path to the folder containing the annotations considered

    Returns
    -------
    BG_merged : array of shape (388,516)
        the merged background mask
    GM_merged : array of shape (388,516)
        the merged grey matter mask
    WM_merged : array of shape (388,516)
        the merged white matter mask
    all_merged : array of shape (388,516)
        the three masks combined (one color = one class)
    """
    WM_merged = np.zeros(WM.shape)
    BG_merged = np.zeros(WM.shape)
    GM_merged = np.zeros(WM.shape)
    all_merged = np.zeros(WM.shape)

    for idx, x in enumerate(WM):
        for idy, y in enumerate(x):
            
            if bg:
                # 1. check if it is background
                if BG[idx, idy] == 255:
                    BG_merged[idx, idy] = 255
                    all_merged[idx, idy] = 128
                
                # 2. if not, check if the pixel is white matter
                elif WM[idx, idy] == 255:
                    WM_merged[idx, idy] = 255
                    all_merged[idx, idy] = 0
                    
                # 3. if not, it is grey matter
                else:
                    GM_merged[idx, idy] = 255
                    all_merged[idx, idy] = 255
            
            else:
                # 1. check if it is background
                if WM[idx, idy] == 255:
                    WM_merged[idx, idy] = 255
                    all_merged[idx, idy] = 0

                # 2. if not, check if the pixel is white matter
                elif BG[idx, idy] == 255:
                    BG_merged[idx, idy] = 255
                    all_merged[idx, idy] = 128
                    
                # 3. if not, it is grey matter
                else:
                    GM_merged[idx, idy] = 255
                    all_merged[idx, idy] = 255

    # save the masks
    save_image(path, WM_merged, 'WM_merged')
    save_image(path, BG_merged, 'BG_merged')
    save_image(path, GM_merged, 'GM_merged')
    new_p = Image.fromarray(all_merged)
    if new_p.mode != 'L':
        new_p = new_p.convert('L')
        new_p.save(os.path.join(path, 'merged.jpeg'))
        new_p.save(os.path.join(path, 'merged.png'))
        
    return BG_merged, WM_merged, GM_merged, all_merged


def save_image(path, img, name):
    """
    save_image is used to save an image as a .png file

    Parameters
    ----------
    path : str
        the path to the folder containing the annotations considered
    img : array of shape (388,516)
        the image to be saved
    name : str
        the name of the image
    """
    new_p = Image.fromarray(img)
    if new_p.mode != 'L':
        new_p = new_p.convert('L')
        new_p.save(os.path.join(path, name + '.png'))
        



###########################################################################################################################
#########################################         2. Create the border masks        #######################################
###########################################################################################################################


def get_depolarization_values(img, MM, normalize = False, iterations_dilation = 5, number_of_dilations = 9, parameter = 'depolarization'):
    """
    get_depolarization_values is the function used to dilatte the image and get the depolarization values in the 
    consecutive regions

    Parameters
    ----------
    img : np.array
        the annotation mask of shape (388,516)
    MM : dict
        the MM for the measurement of interest
    normalize : boolean
        boolean indicating if the values of depolarization should be normalized between 0 and 1
    iterations_dilation : int
        the number of pixels for one dilation step
    number_of_dilations : int
        the number of dilation steps
    
    Returns
    -------
    x : list of int
        the x values
    y : list of int
        the mean values for each dilated region
    y_vals_sorted : list of list
        the values for each dilated regions
    border_pixels : np.array
        the mask used to identify the pixels at the border of white and gray matter
    regions : np.array
        the dilated image
    """
    border_pixels = generate_border_pixel_img(img)
    
    dilated_images = {}
    img_dilation = border_pixels
    kernel = np.ones((3, 3), np.uint8)

    img_BG = np.zeros(img.shape)
    for idx, x in enumerate(img):
        for idy, y in enumerate(x):
            if y == 128:
                img_BG[idx, idy] = 1

    dilated_BG = cv2.dilate(img_BG, kernel, iterations=5)

    for i in range(number_of_dilations):
        if i == 0:
            img_dilation = cv2.dilate(img_dilation, kernel, iterations=iterations_dilation//2)
        else:
            img_dilation = cv2.dilate(img_dilation, kernel, iterations=iterations_dilation)
        dilated_images[i + 1] = img_dilation

    regions = np.zeros(border_pixels.shape)

    for id_img, dilatted_img in dilated_images.items():
        for idx, x in enumerate(dilatted_img):
            for idy, y in enumerate(x):
                if y == 1:
                    if id_img == 1:
                        regions[idx, idy] = id_img
                    else:
                        if regions[idx, idy] == 0:
                            regions[idx, idy] = id_img
                        else:
                            pass

    for idx, x in enumerate(regions):
        for idy, y in enumerate(x):
            if regions[idx, idy] == 1 and img[idx, idy] != 128:
                pass
            if img[idx, idy] == 128:
                regions[idx, idy] = 0
            elif img[idx, idy] == 255 and regions[idx, idy] != 0:
                regions[idx, idy] = number_of_dilations + regions[idx, idy] - 1
            elif img[idx, idy] == 0 and regions[idx, idy] != 0:
                regions[idx, idy] = abs(number_of_dilations - regions[idx, idy]) + 1
                
    depolarization = defaultdict(list)
    retardance = defaultdict(list)
    diattenuation = defaultdict(list)

    for idx, x in enumerate(regions):
        for idy, y in enumerate(x):
            if y != 0 and dilated_BG[idx, idy] == 0:
                depolarization[y].append(MM['totP'][idx, idy])
                retardance[y].append(MM['linR'][idx, idy])
                diattenuation[y].append(MM['totD'][idx, idy])

    x = []
    y = []
    y_vals = []
    if parameter == 'depolarization':
        listed = depolarization
    else:
        listed = retardance
        
    for key, val in listed.items():
        x.append(key)
        y.append(np.nanmean(val))
        y_vals.append(val)
        
    y = [y_s for _, y_s in sorted(zip(x, y))]
    if normalize:
        y = (y - min(y))/(max(y) - min(y))
    
    y_vals_sorted = [y_s for _, y_s in sorted(zip(x, y_vals))]

    x = sorted(x)
    
    kernel = np.ones((3, 3), np.uint8)
    dilated_border_pixels = cv2.dilate(border_pixels, kernel, iterations=2)

    return x, y, y_vals_sorted, dilated_border_pixels, regions


def generate_border_pixel_img(img):
    """
    generate_border_pixel_img creates a mask with the pixel values for the one the border of white 
    and gray matter values set at 1

    Parameters
    ----------
    img : np.array
        the annotation mask of shape (388,516)
    
    Returns
    -------
    border_pixels : np.array
        the mask used to identify the pixels at the border of white and gray matter
    """
    border_pixels = np.zeros(img.shape)

    for idx, x in enumerate(img):
        min_x, max_x, debug = get_idx(idx, 0, img.shape[0])
        for idy, y in enumerate(x):
            if y == 0:
                min_y, max_y, _ = get_idx(idy, 0, img.shape[1])
                
                # get the neighbors
                neighbors = img[min_x:max_x, min_y:max_y]
                
                if sum(sum(neighbors)) == 0:
                    pass
                else:
                    border_gray = True
                    for idx_s, xs in enumerate(neighbors):
                        for idy_s, ys in enumerate(xs):
                            if ys == 128:
                                border_gray = False
                            else:
                                pass
                    if border_gray:
                        border_pixels[idx, idy] = 1
                        
    return border_pixels


def get_idx(idx, distance, shape):
    """
    get_idx is used to find the window that should be used to determine the border pixels

    Parameters
    ----------
    idx : int
        the index of the pixel considered
    distance : int
        the distance to be used for the window
    shape : int
        the length of the image in the considered axis
    
    Returns
    -------
    min_x : int
        the minimal index for the window
    max_x : int
        the maximal index for the window
    debug : boolean
        boolean debug
    """
    debug = False
    
    # if the window would start before the first pixel...
    if idx - distance <= 0:
        min_x = 0
        max_x = idx + distance + 2
        debug = True
        
    # ... or after the last one
    elif idx + distance >= shape - 1:
        min_x = idx - distance - 2
        max_x = shape
    else:
        min_x = idx - distance - 1
        max_x = idx + distance + 2
    return min_x, max_x, debug


def dilate_images(variables, iterations_dilation, number_of_dilations, parameter, color_map):
    """
    perform the image dillation and returns the border pixel map, the values for a parameter of interest as well as the
    dilated images

    Parameters
    ----------
    variables : dict
        a dict containing the MMs
    iterations_dilation : int
        the number of iteration for each dilation step
    number_of_dilations : int
        the number of dilation steps
    parameter : str
        the parameter of interest
    color_map : dict
        the color map to plot the dilated regions
    
    Returns
    -------
    values_combined : dict
        the values of the parameter in the regions created during dilation
    border_pixels : dict
        a dictionnary containing the map of the pixels located at the border
    dilated_images, dilated_images_plot : dict
        the dilated images
    """
    # define the dictionnaries to collect the information
    values_combined = {}
    border_pixels = {}
    dilated_images = {}
    dilated_images_plot = {}
    
    # loop over all the measurements performed
    for var, imgs in tqdm(variables.items()):

        # 4.1. get the annotated image and the mueller matrix for the measurement
        img = imgs['annotation']
        MM = imgs['MM']

        # 4.2. get the parameter values for each dilated region and the corresponding distance values
        x, y, y_vals_sorted, border_pixel, dilated_image = get_depolarization_values(img, MM, 
                                         iterations_dilation = iterations_dilation, number_of_dilations = number_of_dilations, 
                                                                                     parameter = parameter)
        dilated_image_rescaled = copy.deepcopy(dilated_image)
        dilated_image_rescaled = (255/np.max(dilated_image_rescaled) * dilated_image_rescaled).astype(np.uint8)

        # 4.3. save the variables in the dictionnaries
        values_combined[var] = [x, y, y_vals_sorted]
        border_pixels[var] = border_pixel
        dilated_images[var] = dilated_image

        # 4.4. create an image with the different regions colored in different colors
        image_figure = np.zeros((388,516,3))
        for idx, x in enumerate(dilated_image):
            for idy, y in enumerate(x):
                if y == 0:
                    image_figure[idx, idy] = [0, 0, 0]
                else :
                    image_figure[idx, idy] = color_map[y]
        dilated_images_plot[var] = image_figure.astype(np.uint8)
        
    return values_combined, border_pixels, dilated_images, dilated_images_plot



###########################################################################################################################
#######################         2.1. Save the border pixel images and the dilated images        ###########################
###########################################################################################################################


def save_the_images(border_pixels, dilated_images_plot, folder_save):
    """
    save the border pixel images and the dilated images

    Parameters
    ----------
    border_pixels : dict of array of shape (388, 516)
        the border pixel arrays
    dilated_images_plot : dict of array of shape (388, 516)
        the dilatted images to be plotted
    folder_save : str
        the folder in which to save the images
    """
    for var, border_pixel in border_pixels.items():
    
        # create the folders that will host the results, if not existing
        try:
            os.mkdir(os.path.join(folder_save))
        except FileExistsError:
            pass
        try:
            os.mkdir(os.path.join(folder_save, 'border_regions'))
        except FileExistsError:
            pass
        try:
            os.mkdir(os.path.join(folder_save, 'dilatted_imgs'))
        except FileExistsError:
            pass

        # save the border region pixel map
        im = Image.fromarray(border_pixel * 255).convert("L")
        im.save(os.path.join(folder_save, 'border_regions', var.split('/')[-1] + '.png'))

        # and the so-called rainbow plot
        im = Image.fromarray(dilated_images_plot[var])
        im.save(os.path.join(folder_save, 'dilatted_imgs', var.split('/')[-1] + '_rainbow.png'))
        
        

###########################################################################################################################
######################         2.2. Plot the the parameter values according to the distance        ########################
###########################################################################################################################


def get_complete_plot_data(values_combined, folder_save, parameter, iterations_dilation, number_of_dilations):
    """
    get_complete_plot_data allows one to get the complete plot data

    Parameters
    ----------
    values_combined : dict
        the values of the parameter in the regions created during dilation
    parameter : str
        the parameter of interest
    folder_save : str
        the folder in which to save the images
    iterations_dilation : int
        the number of iteration for each dilation step
    number_of_dilations : int
        the number of dilation steps
        
    Returns
    -------
    fig : python figure
        the generated figure
    x_scaled_sorted : list
        the distances from the border zone
    organized_by_distance_means : list
        the means organized according to the distance from the border zone
    """
    # 1. organize the data that we have by distance to the border zone
    vals_distance = defaultdict(list)
    for folder, [x, y, y_vals] in values_combined.items():
        for x_, y_val_ in zip(x, y_vals):
            vals_distance[x_].append(y_val_) 
    xs = []
    for x, y_val in vals_distance.items():
        xs.append(x)
        vals_distance[x] = flatten(y_val)
        
    means, stds = get_data_plot_individual(vals_distance)
    
    # 2. rescale the distances to make them "physical distances"
    xs = list(vals_distance.keys())
    x_scaled = []
    for x in xs:
        x_scaled.append(round(iterations_dilation * (x - number_of_dilations) * 0.05, 2))
        
    # 3. re-organize the means and data to sort them according to growing distance   
    organized_by_distance_means = []
    organized_by_distance_stds = []

    for (x, val), (x, std) in zip(means.items(), stds.items()):
        organized_by_distance_means.append(val)
        organized_by_distance_stds.append(std)

    ys_avg_sorted = [y_s for _, y_s in sorted(zip(x_scaled, organized_by_distance_means))]
    ys_std_sorted = [y_s for _, y_s in sorted(zip(x_scaled, organized_by_distance_stds))]
    x_scaled_sorted = sorted(x_scaled)

    # 4. actually do the plot
    fig = plot_distances_and_std(x_scaled_sorted, organized_by_distance_means, organized_by_distance_stds, folder_save, parameter)
    
    return fig, x_scaled_sorted, organized_by_distance_means


def get_data_plot_individual(times_vals):
    """
    get_data_plot_individual is used to get the means and stds

    Parameters
    ----------
    times_vals : dict of dict
        the values to obtain the parameters for
    
    Returns
    -------
    means : dict of dict
        the means
    stds : dict of dict
        the standard deviations
    """
    means = {}
    stds = {}
 
    for x, vals in tqdm(times_vals.items()):
        means[x] = np.nanmean(vals)
        stds[x] = np.nanstd(vals)
            
    return means, stds


def flatten(l):
    """
    flatten a list
    Parameters
    ----------
    l : list
        the list to flatten
    
    Returns
    -------
    flattened_list : l
        the flattened list
    """
    return [item for sublist in l for item in sublist]


def plot_distances_and_std(x_scaled_sorted, organized_by_distance_means, organized_by_distance_stds, folder_save, parameter):
    """
    generate the plot with the distances and stds

    Parameters
    ----------
    x_scaled_sorted : list
        the distances from the border zone
    organized_by_distance_means : list
        the means organized according to the distance from the border zone
    organized_by_distance_stds : list
        the stds organized according to the distance from the border zone
    folder_save : str
        the folder in which to save the figure
    parameter : str
        the parameter of interest
    
    Returns
    -------
    fig : python figure
        the generated figure
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(x_scaled_sorted, organized_by_distance_means, yerr = organized_by_distance_stds)
    ax.set_xlabel('Distance to border (mm)', fontweight="bold", fontsize=16)
    if parameter == 'retardance':
        ax.set_ylabel('Linear retardance (°)', fontweight="bold", fontsize=16)
        # ax.set_ylim([0.7, 1])
    else:
        ax.set_ylabel('Depolarization', fontweight="bold", fontsize=16)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim([0.7, 1])
        ax.set_xlim([-2.5, 2.5])
  
    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontweight('bold') for label in labels]
    [label.set_fontsize(14) for label in labels]

    ax.text(-2.2, 0.82, 'White matter', fontsize=16, fontweight = 'bold')
    ax.text(0.75, 0.96, 'Grey matter', fontsize=16, fontweight = 'bold')

    ax.annotate("", xy=(-2, 0.8), xytext=(-0.5, 0.8),
                arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(2.2, 0.94), xytext=(1.0, 0.94),
                arrowprops=dict(arrowstyle="->"))

    plt.tight_layout()

    try:
        os.mkdir(os.path.join(folder_save, 'figures'))
    except FileExistsError:
        pass

    plt.savefig(os.path.join(folder_save, 'figures', 'combined_depolarization.png'))
    plt.savefig(os.path.join(folder_save, 'figures', 'combined_depolarization.pdf'))

    return fig



###########################################################################################################################
##########################################         2.3. Try to fit a function        ######################################
###########################################################################################################################

def fit_and_plot_sigmoid(x_scaled_sorted, organized_by_distance_means, folder_save, parameter, wavelength):
    """
    fit a sigmoid and generate the plot

    Parameters
    ----------
    x_scaled_sorted : list
        the distances from the border zone
    organized_by_distance_means : list
        the means organized according to the distance from the border zone
    folder_save : str
        the folder in which to save the figure
    parameter : str
        the parameter of interest
    wavelength : str
        the wavelength of interest
    
    Returns
    -------
    fig : python figure
        the generated figure
    kn_neg : double
        the elbow for white matter
    kn_pos : double
        the elbow for grey matter
    """
    # 1. initial guess for the fitting
    p0 = [max(organized_by_distance_means), np.median(x_scaled_sorted),1,min(organized_by_distance_means)] 

    # 2. fit the curve to th sigmoid
    popt, pcov = curve_fit(sigmoid, x_scaled_sorted, organized_by_distance_means, p0, method='dogbox')
    print(popt, pcov)

    # 3. find the elbows of the function
    x_long = np.arange(min(x_scaled_sorted) - 0.2, 0.5, 0.0001)
    y_fitted = sigmoid(x_long, *popt)
    kn_neg = KneeLocator(x_long, y_fitted, curve='concave', direction='decreasing').knee

    x_long = np.arange(0, max(x_scaled_sorted) + 0.2, 0.0001)
    y_fitted = sigmoid(x_long, *popt)
    kn_pos = KneeLocator(x_long, y_fitted, curve='convex', direction='decreasing').knee

    # 4. and prepare the data points to plot it
    x_long = np.arange(min(x_scaled_sorted) - 0.2, max(x_scaled_sorted) + 0.2, 0.1)
    y_fitted = sigmoid(x_long, *popt)

    # 5. actually plot
    fig = plot_with_sigmoid(x_scaled_sorted, organized_by_distance_means, x_long, y_fitted, folder_save, parameter, wavelength,
                           kn_neg, kn_pos)
    
    return fig, [kn_neg, kn_pos]

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)


def plot_with_sigmoid(x_scaled_sorted, organized_by_distance_means, x_long, y_fitted, folder_save, parameter, wavelength, kn_neg, kn_pos):
    """
    generate the plot with the distances and stds

    Parameters
    ----------
    x_scaled_sorted : list
        the distances from the border zone
    organized_by_distance_means : list
        the means organized according to the distance from the border zone
    x_long : list
        the xs used to plot the sigmoid function
    y_fitted : list
        the fitted sigmoid
    folder_save : str
        the folder in which to save the figure
    parameter : str
        the parameter of interest
    wavelength : str
        the wavelength of interest
    kn_neg : double
        the elbow for white matter
    kn_pos : double
        the elbow for grey matter
    
    Returns
    -------
    fig : python figure
        the generated figure
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_scaled_sorted, organized_by_distance_means, 'o', label='data')
    ax.plot(x_long, y_fitted, label='fitted')

    if parameter == 'depolarization':
        if wavelength == '550nm':
            neg_min = kn_neg 
            neg_max = kn_pos

    # negative part
    ax.axvspan(neg_min, neg_max, alpha = 0.15, color = '#ffb618')
    ax.axvspan(-10, neg_min, alpha = 0.15, color = '#52a736')
    ax.axvspan(neg_max, 10, alpha = 0.15, color = '#FF3030')


    ax.set_xlabel('Distance to border (mm)', fontweight="bold", fontsize=16)
    if parameter == 'retardance':
        ax.set_ylabel('Linear retardance (°)', fontweight="bold", fontsize=16)
        # ax.set_ylim([0.7, 1])
    else:
        ax.set_ylabel('Depolarization', fontweight="bold", fontsize=16)
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.set_ylim([0.7, 1])
        ax.set_xlim([-2.5, 2.5])

    labels = ax.get_xticklabels() + ax.get_yticklabels()
    [label.set_fontweight('bold') for label in labels]
    [label.set_fontsize(14) for label in labels]

    ax.text(-2.3, 0.82, 'White matter', fontsize=16, fontweight = 'bold')
    ax.text(1, 0.96, 'Grey matter', fontsize=16, fontweight = 'bold')

    ax.annotate("", xy=(-2.2, 0.8), xytext=(-0.8, 0.8),
                arrowprops=dict(arrowstyle="->"))
    ax.annotate("", xy=(2.3, 0.94), xytext=(1.3, 0.94),
                arrowprops=dict(arrowstyle="->"))

    plt.tight_layout()

    plt.savefig(os.path.join(folder_save, 'figures', 'fitted_depolarization.png'))
    plt.savefig(os.path.join(folder_save, 'figures', 'fitted_depolarization.pdf'))
    
    return fig




###########################################################################################################################
##############################         3. Generate the transition region on the image        ##############################
###########################################################################################################################

def generate_transition_region_on_the_image(border_pixel, dilated_image, variables, distance_pixel, elbows, var, folder_save):
    """
    function used to plot the transition region on an image

    Parameters
    ----------
    border_pixels : dict
        a dictionnary containing the map of the pixels located at the border
    dilated_images : dict
        the dilated images
    variables : dict of var
        the annotated image, MM and grayscale image for each measurement folder
    distance_pixel : double
        the distance for each pixel
    elbows : list
        the elbows for white and grey matter
    var : str
        the folder of interest
    folder_save : str
        the folder in which to save the image
    
    Returns
    -------
    new_img : image
        the created image
    
    """
    
    [kn_neg, kn_pos] = elbows
    border_pixels_ = generate_uncertainty(border_pixel, dilated_image, kn_neg, kn_pos, distance_pixel)
    
    uncertainty_region = np.zeros((388, 516, 3))

    for idx, x in tqdm(enumerate(border_pixels_), total = len(border_pixels_)):
        for idy, y in enumerate(x):
            if y != 0:
                uncertainty_region[idx, idy] = [255, 182, 24]
            else:
                if dilated_image[idx, idy] != 0:
                    if dilated_image[idx, idy] <= 6:
                        uncertainty_region[idx, idy] = [82,167,54]
                    elif dilated_image[idx, idy] > 6:
                        uncertainty_region[idx, idy] = [255,48,48]
                else:
                    im = variables[var]['grayscale'][idx, idy]
                    uncertainty_region[idx, idy] = [im[0],im[1],im[2]]
                    
    im = Image.fromarray(uncertainty_region.astype(np.uint8))
    
    background = im
    overlay = Image.fromarray(variables[var]['grayscale'])
    background = background.convert("RGBA")
    overlay = overlay.convert("RGBA")

    new_img = Image.blend(background, overlay, 0.8)
    new_img.save(os.path.join(folder_save, "uncertainty_rainbow.png"),"PNG")
    
    return new_img


def get_idx(idx, distance, shape, add = True):
    """
    get_idx is used to find the window that should be used to determine the border pixels

    Parameters
    ----------
    idx : int
        the index of the pixel considered
    distance : int
        the distance to be used for the window
    shape : int
        the length of the image in the considered axis
    
    Returns
    -------
    min_x : int
        the minimal index for the window
    max_x : int
        the maximal index for the window
    debug : boolean
        boolean debug
    """
    debug = False
    
    # if the window would start before the first pixel...
    if idx - distance <= 0:
        min_x = 0
        
        if add:
            max_x = idx + distance + 2
        else:
            max_x = idx + distance
            
        debug = True
        
    # ... or after the last one
    elif idx + distance >= shape - 1:
        if add:
            min_x = idx - distance - 2
        else:
            min_x = idx - distance
        max_x = shape
    else:
        if add:
            min_x = idx - distance - 1
            max_x = idx + distance + 2
        else:
            min_x = idx - distance
            max_x = idx + distance

    return min_x, max_x, debug


def generate_uncertainty(img, dilated_image, distance_white, distance_grey, distance_pixel):
    """
    generate_border_pixel_img creates a mask with the pixel values for the one the border of white 
    and gray matter values set at 1

    Parameters
    ----------
    img : np.array
        the annotation mask of shape 516x388
    
    Returns
    -------
    border_pixels : np.array
        the mask used to identify the pixels at the border of white and gray matter
    """
    border_pixels = np.zeros(img.shape)
    distance = int(abs(distance_white)/distance_pixel)

    # white matter
    for idx, x in tqdm(enumerate(img), total = len(img)):
        min_x, max_x, debug = get_idx(idx, distance, img.shape[0], add = False)
        
        for idy, y in enumerate(x):
            
            min_y, max_y, _ = get_idx(idy, distance, img.shape[1], add = False)
                
            # get the neighbors
            neighbors = img[min_x:max_x, min_y:max_y]
            
            if dilated_image[idx, idy] <= 6: 
                try:
                    if sum(sum(neighbors)) != 0:
                        border_pixels[idx, idy] = 150
                except TypeError:
                    pass
                    
               
            
    distance = int(abs(distance_grey)/distance_pixel)
    
    # grey matter
    for idx, x in tqdm(enumerate(img), total = len(img)):
        min_x, max_x, debug = get_idx(idx, distance, img.shape[0], add = False)
        
        for idy, y in enumerate(x):
            
            min_y, max_y, _ = get_idx(idy, distance, img.shape[1], add = False)
                
            # get the neighbors
            neighbors = img[min_x:max_x, min_y:max_y]
            
            if dilated_image[idx, idy] > 6: 
                try:
                    if sum(sum(neighbors)) != 0:
                        border_pixels[idx, idy] = 150
                except TypeError:
                    pass
    return border_pixels

def make_overlay(background, overlay, alpha = 0.1):
    """
    creates an overlay of background and overlay image

    Parameters
    ----------
    background : Pillow image
        the image to be used as background
    overlay : Pillow image
        the image to overlay on the background
    fname : str
        the path in which the image should be saved
    pol_maps : boolean
        indication of if we are working with polarimetric maps (in this case, crop the image - default = False)
    alpha : double
        the transparency level (default = 0.1)
    """
    overlay = overlay.convert('RGBA')
    new_img = Image.blend(background, overlay, alpha)
    return new_img