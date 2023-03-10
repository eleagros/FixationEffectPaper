import os
import copy
import pickle

import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict
import math
import numpy as np

from scipy import stats
from scipy.stats import mannwhitneyu, wilcoxon

import re



###########################################################################################################################
###########################         1. Obtain the parameters and data folders of interest        ##########################
###########################################################################################################################

def get_parameters():
    """
    returns the different parameters to be used to run the script depending on the experiment considered
    
    Returns
    -------
    path_folder : str
        the path to the folder containing the data
    wavelength : str
        the wavelength to be used for the analysis
    match_sequence : str
        the match sequence to be found in the folders to be analyzed (i.e. 'PIG-CUSA')
    measurement : str
        the name of the output folder in the data subfolder
    """
    path_folder = './data/fixation_over_time/'
    wavelength = '550nm'
    match_sequence = 'T_P-T0_'
    measurement = 'comparaison'
    max_ = 25
        
    return [path_folder, wavelength, match_sequence, measurement, max_]


def get_data_folders(path_folder, wavelength = '550nm', match_sequence = 'PIG-CUSA'):
    """
    get_data_folders is used to find all the folders containing the match_sequence given as an input

    Parameters
    ----------
    path_folder : str
        the path to the folder containing all the subfolders
    wavelength : str
        the wavelength to be considered (by default '550nm')
    match_sequence : str
        the sequence that needs to be matched (by default 'PIG_CUSA')

    Returns
    -------
    data_folders : list of str
        all the subsubfolders containing the xlsx files produced with the interactive square selector
    base_folders : list of str
        all the subfolders containing the match sequence
    to_be_removed : dict of list
        the subfolders that should be removed for further analyses (i.e. switching from grey to white or full
        of invalid pixels)
    """
    folders = os.listdir(path_folder)

    # check for the folders containing the match sequence
    base_folders = []
    for folder in folders:
        if re.findall(match_sequence, folder):
            if 'P-T0_FX' in folder:
                base_folders.append(os.path.join(path_folder, folder))


    
    # define the matching sequence for the first measurement folder
    matching = 'T0'
    
    # check for subfolders containing the xlsx files produced with the interactive square selector
    data_folders = []
    to_be_removed = {}
    
    for folder in base_folders:
        
        # if the folder is the first measurement
        if matching in folder:
            path50x50 = os.path.join(folder, 'polarimetry', wavelength, '50x50_images')

            # add the measurement that should be removed to the dict to_be_removed
            with open(os.path.join(path50x50, 'to_rm_GM.pickle'), 'rb') as handle:
                rm_gm = pickle.load(handle)
            with open(os.path.join(path50x50, 'to_rm_WM.pickle'), 'rb') as handle:
                rm_wm = pickle.load(handle)
            to_be_removed[folder] = [rm_gm, rm_wm]

            for folder_50x50 in os.listdir(path50x50):
                
                # check if the ROI of interest should be removed
                to_remove = False
                
                # if the folder should be removed - add it to the to_be_removed dictionnary
                if to_remove:
                    if fixation:
                        raise ValueError
                    else:
                        pass
                
                
                # if not, add it to the data_folders list and load the data in the next function
                else:
                    if matching in folder:
                        if os.path.isdir(os.path.join(path50x50, folder_50x50)):
                            data_folders.append(os.path.join(path50x50, folder_50x50))

    return data_folders, base_folders, to_be_removed




###########################################################################################################################
##############################         2. Get the data obtained during the measurement        #############################
###########################################################################################################################


def get_data(data_folders, data_types = ['GM', 'WM']):
    """
    get_data is used to retrieve all the data that was produced with the interactive square selector

    Parameters
    ----------
    data_folders : list of str
        all the subsubfolders containing the xlsx files produced with the interactive square selector
    data_types : list of str
        the different types of data existing (by default ['GM', 'WM']

    Returns
    -------
    data_types : list of str
        the different types of data existing given as an input
    data_fixation : defaultdict
        all the data produced with the interactive square selector in a dict format with type and folder names as keys
    """
    # create the structure that will receive the data
    data_fixation = defaultdict(list)
    for type_ in data_types:
        data_fixation[type_] = dict()

    match_sequence = '[A-Z]{2}'

    # iterate over each data folder
    for folder in data_folders:
        
        # find the type of the tissue
        type_ = re.findall(match_sequence, folder.split('\\')[-1])
        
        assert len(type_) == 1
        type_ = type_[0]
        
        # and add to the dictionnary the measured values
        data_fixation[type_][folder] = read_xlsx(os.path.join(folder, 'summary_aligned.xlsx'))
    return data_types, data_fixation

def read_xlsx(path):
    """
    read_xlsx is used to read excel files into pandas dataframes

    Parameters
    ----------
    path : str
        the path to the excel file
        
    Returns
    -------
    df : pandas dataframe
        the loaded dataframe
    """
    df = pd.read_excel(path)
    
    fname = None

    for idx, row in df.iterrows():
        if type(row['fname']) == str:
            fname = row['fname']
        else:
            
            # add the missing filenames
            df.at[idx, 'fname'] = fname
    return df

def find_match_seq(row, match_seq):
    """
    find_match_seq is used in an apply loop to find if match_seq can be found in the filename, and if yes put all
    the values to NaN

    Parameters
    ----------
    row : str
        the row of the measured dataframe
    match_seq : str
        the match_sequence to be found

    Returns
    -------
    boolean
        indicating wether or not match_seq could be found in the fname of the row
    """
    if match_seq in row['fname']:
        for col in row.items():
            if col[0] == 'fname' or col[0] == 'parameter' or col[0] == 'square_size':
                pass
            else:
                row[col[0]] = math.nan
        return row
    else:
        return row
    
    
def create_output_directories(measurement, data_types, 
                              parameters = ['retardance', 'depolarization', 'diattenuation', 'azimuth']):
    """
    create_output_directories generated the structure that will be used to store the output of the analysis

    Parameters
    ----------
    measurement : str
        the name of the output folder in the data subfolder
    data_types : list of str
        the different types of data existing
    parameters : list of str
        the list of the parameters to be studied (default = ['retardance', 'depolarization', 'diattenuation', 'azimuth'])
    
    Returns
    -------
    folder : str
        the path to the folder in which data should be saved
    """
    # create the main subfolder
    try:
        folder = os.path.join('results', measurement)
        os.mkdir(folder)
    except FileExistsError:
        pass

    for data_type in data_types:
        # create a folder for each data_type...
        try:
            folder_data_type = os.path.join(folder, data_type)
            os.mkdir(folder_data_type)
        except FileExistsError:
            pass
        
        # ... and for each parameter
        for parameter in parameters:
            try:
                folder_data_parameter = os.path.join(folder_data_type, parameter)
                os.mkdir(folder_data_parameter)
            except FileExistsError:
                pass
    return folder


def get_and_organize_data(data_fixation, metrics, data_types, times, path_folder,
                          parameters = ['retardance', 'depolarization', 'diattenuation', 'azimuth']):
    """
    get_and_organize_data is used to organize the data and save it into xlsx file

    Parameters
    ----------
    data_fixation : dict of dict of df
        the data from the different measurements
    metrics : list of str
        the list of metrics to be used (i.e. ['mean', 'median', 'max'])
    data_types : list of str
        the different types of data existing
    times : list of str
        the list of times of the measurements (i.e. ['T0', 'D12', 'D24', 'D36', 'D48', 'D7d'])
    path_folder : str
        the path to the folder in which data should be saved
    parameters : list of str
        the list of the parameters to be studied (default = ['retardance', 'depolarization', 'diattenuation', 'azimuth'])
        
    Returns
    -------
    to_combine : dict
        the data from different measurements combined with a tuple key containing the parameter, the tissue type and the metric
    """
    to_combine = {}

    for parameter in parameters:
        for metric in metrics:
            data_organized = {}
            data_organized_normalized = {}
            for data_type in data_types:
                data_organized[data_type], data_organized_normalized[data_type] = get_data_per_time(data_fixation[data_type], 
                                                                                    times, parameter, metric)
            
            for data_type, dfs in data_organized.items():
                df = create_df(dfs)
                df_normalized = create_df(data_organized_normalized[data_type])
                
                folder_data_type = os.path.join(path_folder, data_type)
                folder_data_parameter = os.path.join(folder_data_type, parameter)

                path_df = os.path.join(folder_data_parameter, metric + '_raw_data.xlsx')
                df.to_excel(path_df)

                path_df_normalized = os.path.join(folder_data_parameter, metric + '_normalized.xlsx')
                df_normalized.to_excel(path_df_normalized)
                
                if (parameter == 'retardance' or 
                                         parameter == 'depolarization' or parameter == 'azimuth'):
                    while(df.shape[1] < 100):
                        df.insert(df.shape[1], df.shape[1], [math.nan, math.nan, math.nan, math.nan, math.nan, math.nan])
                    
                    to_combine[(parameter, data_type, metric)] = df
                
    return to_combine


def get_data_per_time(test_df, times, parameter, metric):
    """
    get_data_per_time is used to reorganize the data obtained by time

    Parameters
    ----------
    test_df : dict
        all the data produced with the interactive square selector in a dict format with type and folder names as keys
    times : list of str
        the list of times of the measurements (i.e. ['T0', 'D12', 'D24', 'D36', 'D48', 'D7d'])
    parameter : str
        the polarimetric parameter to be studied (i.e. 'retardance', 'depolarization', 'diattenuation', 'azimuth')
    parameter : str
        the metric of interest

    Returns
    -------
    data_time : dict
        all the data produced with the interactive square selector organized by times
    data_time_normalized : dict
        data_time but with normalized data
    """
    # create the stuctures that will receive the data
    data_time = dict()
    for time in times:
        data_time[time] = list()
        
    data_time_normalized = dict()
    for time in times:
        data_time_normalized[time] = list()

    # iterate over each df
    for fname, df in test_df.items():
        df['fname'] = find_names(list(df['fname']))
        
        # search for the data relative to the parameters of interest
        parameter_df = df[df['parameter'] == parameter]
        
        reference = 0
        for idx, row in parameter_df.iterrows():
            match = []
            for time in times:
                if re.findall(time, row['fname']):
                    match.append(re.findall(time, row['fname']))
                    
            try:
                assert len(match) == 1 or len(match) == 2
            
                # if necessary, add as reference
                if match[-1][0] == times[0]:
                    reference = row[metric]

                if math.isnan(row[metric]):
                    data_time[match[-1][0]].append(row[metric])
                    data_time_normalized[match[-1][0]].append(row[metric])
                else:
                    data_time[match[-1][0]].append(row[metric])
                    try:
                        data_time_normalized[match[-1][0]].append(row[metric]/reference)
                    except:
                        data_time_normalized[match[-1][0]].append(math.nan)
            except:
                pass # print(row['fname'])
                
    return data_time, data_time_normalized


def find_names(fnames):
    """
    find_names is used to match the actual file names with the convention 

    Parameters
    ----------
    fnames : list of str
        the list containing all the file names

    Returns
    -------
    fnames : list of str
        the modified file names
    """
    numbers = []
    for file in fnames:
        try:
            numbers.append(int(file.split('_')[-1]))
        except:
            pass
        
    if min(numbers) > 3 and max(numbers) > 3:
        to_modify = {'M_1': 'M_' + str(min(numbers)), 'M_3': 'M_' + str(max(numbers))}

        modified = []
        for file in fnames:
            new_f = file
            for source, new in to_modify.items():
                new_f = new_f.replace(new, source)
            modified.append(new_f)

        return modified
    else:
        return fnames
    
def create_df(dfs):
    """
    create_df is used to create a dataframe from a dictionnary

    Parameters
    ----------
    dfs : dict
        the dictionnary to be converted to a dataframe

    Returns
    -------
    df : pd.DataFrame
        the converted dataframe
    """
    values = []
    for time_point, val in dfs.items():
        values.append(val)
    df = pd.DataFrame(values)
    return df

def average_angles(angles):
    """
    computes the average (mean) of angles
    
    Parameters
    ----------
    angles : list
        the list of angles: range [0, 180]

    Returns
    -------
    the average angle
    """
    x = sum(math.cos(a * 2*math.pi/180) for a in angles)
    y = sum(math.sin(a * 2*math.pi/180) for a in angles)

    if x == 0 and y == 0:
        raise ValueError(
            "The angle average of the inputs is undefined: %r" % angles)

    # To get outputs from -pi to +pi, delete everything but math.atan2() here.
    return 180 * math.fmod(math.atan2(y, x) + 2 * math.pi, 2 * math.pi) / (2*math.pi)


def subtract_angles(column):
    """
    subtract_angles is used to get the difference of the angles in a list with the first one
    
    Parameters
    ----------
    column : list
        the list of angles: range [0, 180]

    Returns
    -------
    the differences of the angles in a list with the first one
    """
    original_angle = column[0]
    res = []
    for c in column:
        res.append(subtract_angle(original_angle, c))
    return res

def subtract_angle(targetA, sourceA):
    """
    subtract_angle substract two angles
    
    Parameters
    ----------
    targetA : int
        the first angle: range [0, 180]
    sourceA : int
        the second angle: range [0, 180]

    Returns
    -------
    the differences of the angles
    """
    a = targetA - sourceA
    return abs((a + 90) % 180 - 90)

def determine_outlier_thresholds_std(values):
    """
    determine_outlier_thresholds_std allows to determine outlier values (located at +- 2 stds from the mean) - not used in the analysis
    
    Parameters
    ----------
    values : list
        the series of values for which to find the outliers

    Returns
    -------
        the lower and upper boundaries for the outliers
    """
    upper_boundary = values.mean() + 2 * values.std()
    lower_boundary = values.mean() - 2 * values.std()
    return lower_boundary, upper_boundary
    
def remove_outliers(df_cv):
    """
    remove_outliers allows to check column to remove to remove outlier values (located at +- 2 stds from the mean) - not used in the analysis
    
    Parameters
    ----------
    df_cv : pandas dataframe
        the dataframe for which to check the presnce of outliers

    Returns
    -------
    df_cv : pandas dataframe
        the dataframe without the outliers
    ~cv_big : boolean list
        the columns to remove
    """
    lower_boundary, upper_boundary = determine_outlier_thresholds_std(df_cv.loc[1])
    cv_big = np.logical_or(df_cv.loc[1] > upper_boundary, df_cv.loc[1] < lower_boundary)

    lower_boundary, upper_boundary = determine_outlier_thresholds_std(df_cv.loc[2])
    cv_big = np.logical_or(cv_big, np.logical_or(df_cv.loc[2] > upper_boundary, df_cv.loc[2] < lower_boundary))

    lower_boundary, upper_boundary = determine_outlier_thresholds_std(df_cv.loc[3])
    cv_big = np.logical_or(cv_big, np.logical_or(df_cv.loc[3] > upper_boundary, df_cv.loc[3] < lower_boundary))

    lower_boundary, upper_boundary = determine_outlier_thresholds_std(df_cv.loc[4])
    cv_big = np.logical_or(cv_big, np.logical_or(df_cv.loc[4] > upper_boundary, df_cv.loc[4] < lower_boundary))

    lower_boundary, upper_boundary = determine_outlier_thresholds_std(df_cv.loc[5])
    cv_big = np.logical_or(cv_big, np.logical_or(df_cv.loc[5] > upper_boundary, df_cv.loc[5] < lower_boundary))
    
    return df_cv, ~cv_big

def combine_data_cv(df, azimuth = False):
    """
    combine_data_cv allows to remove outlier values (not used in the analysis) and prepare the data for fold change analysis
    
    Parameters
    ----------
    df : pandas dataframe
        the dataframe to convert
    azimuth : boolean
        boolan indicating if we are working with azimuth values (default: False)

    Returns
    -------
    df_cv : pandas dataframe
        the dataframe ready for fold change analysis
    """
    if azimuth:
        df_substracted = df.apply(lambda x: subtract_angles(x), axis=0)
    else:
        df_substracted = df/df.loc[0]
    df_cv, cv_big = remove_outliers(df_substracted)
    for idx, cv in enumerate(cv_big):
        if cv:
            pass
        else:
            pass
            # df_cv[idx] = [1.0,math.nan,math.nan,math.nan,math.nan,math.nan]
        
    return df_cv



###########################################################################################################################
##############################         3. Get the data to be used for the comparisons        ##############################
###########################################################################################################################

def recombine_data_tests(to_combine):
    """
    recombine_data_tests is used to recombine the data in a different format that is suited for the function t_test

    Parameters
    ----------
    to_combine : dict
        the data to be recombined

    Returns
    -------
    data_all_recombined : dict of dict of dict
        the recombined data in a format suited for the function t_test
    """
    data_all_recombined = {'GM':{'retardance': {}, 'depolarization': {}, 'azimuth': {}},
                      'WM':{'retardance': {}, 'depolarization': {}, 'azimuth': {}}}
    for idx, val in to_combine.items():
        val_T = copy.deepcopy(val.T)
        val_T.columns = ['T0', 'D12', 'D24', 'D36', 'D48', 'D7d']
        data_all_recombined[idx[1]][idx[0]][idx[2]] = val_T
    return data_all_recombined


def t_test(data_all, times, param_interest):
    """
    perform the t-test to determine if the parameters have changed after CUSA or not

    Parameters
    ----------
    data_all : dict
        a dictionnary containing all the data previously computed for each data type and each parameters
    times : list of str
        the list of times to be considered (i.e. M_1, M_3)
    param_interest : str
        the parameter of interest

    Returns
    -------
    paired_t_test : dict
        a dictionnary containing the results of the statistical tests
    """
    paired_t_test = {}
    
    # iterate over all the data types...
    for data_type, values in data_all.items():

        paired_t_test_param = {}
        
        # ...over all the parameters...
        for parameter, value in values.items():
            
            # do nothing for the azimuth
            if parameter == 'azimuth':
                
                paired_t_test_metric = {}

                # iterate over each metric
                for metric, val in value.items():
                    
                    if metric == 'mean':
                        paired_t_test_time = {}

                        # and perform the t_test, comparing either the two series or the CV with 0
                        for time in times:
                            dof = len(val.dropna()[[times[1], time]]) - 1
                            before = val[times[1]].dropna()
                            after = val[time].dropna()
                            n_1 = len(before)
                            n_2 = len(after)

                            stat, p_val = mannwhitneyu(before, after)
                            paired_t_test_time[time] = [stat, p_val, n_1, n_2, np.mean(before), 
                                                            np.std(before), np.mean(after), np.std(after)]
                        paired_t_test_metric[metric] = paired_t_test_time
                    else:
                        pass
                    
                paired_t_test_param[parameter] =  paired_t_test_metric
                
                
            else:
                paired_t_test_metric = {}

                # iterate over each metric
                for metric, val in value.items():
                    
                    if metric == param_interest:
                        paired_t_test_time = {}

                        # and perform the t_test, comparing either the two series or the CV with 0
                        for time in times:
                            dof = len(val.dropna()[[times[1], time]]) - 1
                            before = val[times[1]].dropna()
                            after = val[time].dropna()
                            n_1 = len(before)
                            n_2 = len(after)

                            stat, p_val = mannwhitneyu(before, after)
                            
                            before = val[times[1]].dropna()
                            after = val[time].dropna()
                            paired_t_test_time[time] = [stat, p_val, n_1, n_2, np.mean(before), 
                                                            np.std(before), np.mean(after), np.std(after)]
                                                            # stats.ttest_rel(before, after).pvalue, dof,
                                                            # np.mean(before), 
                                                                   # np.std(before), np.mean(after), np.std(after)]
                            # paired_t_test_time[time] = [stats.ttest_rel(before, after).statistic, 
                                                            # stats.ttest_rel(before, after).pvalue , dof,
                                                            # np.mean(before), 
                                                                   # np.std(before), np.mean(after), np.std(after)]
                        paired_t_test_metric[metric] = paired_t_test_time
                    else:
                        pass
                    
                paired_t_test_param[parameter] =  paired_t_test_metric
                
            paired_t_test[data_type] = paired_t_test_param
            
    return paired_t_test


import traceback

def create_test_df(paired_t_test, parameter = 'median'):
    """
    create_test_df is used to summarize the results of the tests in a pandas dataframe

    Parameters
    ----------
    paired_t_test : dict
        a dictionnary containing the results of the statistical tests
    parameter : str
        the parameter of interest (default 'median')

    Returns
    -------
    df_grouped : pandas dataframe
        the results of the test
    """
    dfs = []
    
    try:
        # iterate over each data_type
        for data_type, values in paired_t_test.items():

            # and each parameters
            for param, value in values.items():

                # val = []
                # for _, va in value[parameter].items():
                    # val.append(va[1])
                try:
                    df = pd.DataFrame.from_dict(value[parameter]).T
                except:
                    assert param == 'azimuth'
                    df = pd.DataFrame.from_dict(value['mean']).T


                df.columns = ['Z', 'p_value', 'n_1', 'n_2', 'before/GM', 'before/GM std', 'after/WM', 'after/WM std']
                df['parameter'] = param
                df['data_type'] = data_type
                df['time'] = df.index

                dfs.append(df)
    except:
        traceback.print_exc()
            
    df = pd.concat(dfs)
    df_grouped = df.set_index(['parameter', 'data_type', 'time']).sort_values(['parameter', 'data_type'])
    return df_grouped


def create_df_prism(data_dict, max_):
    """
    create_test_df is used to create the dataframe that was used to create the plots in prism

    Parameters
    ----------
    data_dict : dict
        the data in the form of a dictionnary
    max_ : int
        the number of ROI per image

    Returns
    -------
    df : pandas dataframe
        the dataframe that can be used to generate the plots in prism
    """
    data = np.hstack([data_dict['GM'], data_dict['WM']])
    df = pd.DataFrame(data)
    return df