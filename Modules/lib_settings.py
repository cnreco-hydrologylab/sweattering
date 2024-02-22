"""
Library Features:

Name:          lib_settings
Author(s):     Martina Natali
Date:          '20240212'
Version:       '1.0.0'
"""
# ----------------------------------------------------------------------------
# libraries
import os
import json
import logging
# ----------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to get data settings
def get_data_settings(file_name):
    if os.path.exists(file_name):
        with open(file_name) as file_handle:
            data_settings = json.load(file_handle)
    else:
        logging.error(' ===> Error in reading settings file "' + file_name + '"')
        raise IOError('File not found')
    return data_settings
# -------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------
# method to parse data settings
def parse_data_settings(data_settings):

    # Retrieve variables from the 'options' dictionary
    verbose = data_settings['options']['verbose']
    units = data_settings['options']['units']
    opt_veg = data_settings['options']['opt_veg']
    opt_calib = data_settings['options']['opt_calib']
    opt_cost = data_settings['options']['opt_cost']
    opt_field = data_settings['options']['opt_field']
    opt_year = data_settings['options']['opt_year']

    # Retrieve variables from the 'paths' dictionary
    file_shapefile = data_settings['paths']['data_input']['file_shapefile']
    file_meteo = data_settings['paths']['data_input']['file_meteo']
    file_temp = data_settings['paths']['data_input']['file_temp']
    file_pet = data_settings['paths']['data_input']['file_pet']
    file_ndvi = data_settings['paths']['data_input']['file_ndvi']
    file_sm = data_settings['paths']['data_input']['file_sm']
    file_s0 = data_settings['paths']['data_input']['file_s0']

    # Retrieve variables from the 'paths' -> 'data_output' dictionary
    root = data_settings['paths']['root']
    folder_plot = data_settings['paths']['data_output']['folder_plot']
    filename_template = data_settings['paths']['data_output'][
        'filename_template']
    filename_params = data_settings['paths']['data_output'][
        'filename_params']
    filename_table_machine = data_settings['paths']['data_output'][
        'filename_table_machine']
    filename_table_params = data_settings['paths']['data_output'][
        'filename_table_params']
    filename_triple = data_settings['paths']['data_output'][
        'filename_triple']
    filename_scatter = data_settings['paths']['data_output'][
        'filename_scatter']
    extension_plot = data_settings['paths']['data_output']['extension_plot']
    add_description = data_settings['paths']['data_output'][
        'add_description']

    # Retrieve variables from the 'field_params' dictionary
    sand = data_settings['field_params']['sand']
    clay = data_settings['field_params']['clay']
    total_organic_carbon = data_settings['field_params'][
        'total_organic_carbon']
    saturation = data_settings['field_params']['saturation']
    field_capacity = data_settings['field_params']['field_capacity']
    wilting_point = data_settings['field_params']['wilting_point']
    rho_st = data_settings['field_params']['rho_st']

    # Retrieve variables from the 'backscattering' dictionary
    freq = data_settings['backscattering']['freq']

    # Retrieve variables from the 'calibration' dictionary
    PAR = data_settings['calibration']['PAR']
    PAR_str = data_settings['calibration']['PAR_str']
    bounds_up = data_settings['calibration']['bounds']['up']
    bounds_low = data_settings['calibration']['bounds']['low']
    nrun = data_settings['calibration']['run_params']['nrun']
    n_particles = data_settings['calibration']['run_params']['n_particles']
    n_step = data_settings['calibration']['run_params']['n_step']
    optim = data_settings['calibration']['run_params']['optim']
    norma = data_settings['calibration']['run_params']['norma']
    verbose_calib = data_settings['calibration']['run_params'][
        'verbose_calib']
    automate = data_settings['automate']['flag']

    return data_settings

# -------------------------------------------------------------------------------------

def substitute_keywords(template, **kwargs):
    for key, value in kwargs.items():
        template = template.replace('{' + key + '}', str(value))
    return template

# -------------------------------------------------------------------------------------