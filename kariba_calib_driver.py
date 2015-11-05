import os
import json
import multiprocessing as mp
import threading as th
from Queue import Queue 
import time
import subprocess as sp

from dtk.utils.simulation.OutputParser import CompsDTKOutputParser as parser
from dtk.utils.simulation.COMPSJavaInterop import Experiment, QueryCriteria, Client, Configuration, Priority
from dtk.utils.core.DTKSetupParser import DTKSetupParser

from utils import copy_all 
from plot_utils import PlotUtils
from viz_config import VizConfig
from kariba_utils import combine_tags_reports

from kariba_settings import num_procs, sim_data_dir, best_fits_file, all_fits_file, fit_terms_file, calibration_data_file, traces_plots_dir, traces_base_file_name, cc_traces_plots_dir, cc_traces_base_file_name, err_surfaces_plots_dir, err_surfaces_base_file_name, cc_penalty_model, kariba_viz_dir, residuals_file

def multi_proc_run(sweep_name, sweep, command):
    
    root_sweep_dir = os.path.join(sim_data_dir, sweep_name)
    if not os.path.exists(root_sweep_dir):
        os.mkdir(root_sweep_dir)
        
    viz_root_sweep_dir = os.path.join(kariba_viz_dir, sweep_name)
    if not os.path.exists(viz_root_sweep_dir):
        os.mkdir(viz_root_sweep_dir)
    
    procs_available = num_procs
    
    categories = sweep.keys()
    print categories
    
    cat_files = []
    cat_out_files = []
    cat_idx = 0
    
    while cat_idx < len(categories):
    
        category = categories[cat_idx]
        sweep_dir = os.path.join(sim_data_dir,category)
        if not os.path.exists(sweep_dir):
            os.mkdir(sweep_dir)
        cat_f_name = os.path.join(sweep_dir,category+'.json')
        cat_out_f_name = os.path.join(sweep_dir,category+'_out.txt')

        with open(cat_f_name, 'w') as cat_f:
            cat_record = {}
            cat_record['categories'] = {}
            cat_record['categories'][category] = sweep[category]
            cat_record['sweep_dir'] = sweep_dir
            cat_record['root_sweep_dir'] = root_sweep_dir
            cat_record['viz_root_sweep_dir'] = viz_root_sweep_dir 
            json.dump(cat_record, cat_f)
            
        cat_files.append(cat_f_name)
        cat_out_files.append(cat_out_f_name)
        
        procs_available = procs_available - 1
        cat_idx = cat_idx + 1
        
        if procs_available == 0 or cat_idx == len(categories):
            processes = [] 
            for idx,cat_file in enumerate(cat_files):
                with open(cat_out_files[idx], 'w') as out_f:
                    processes.append( sp.Popen(['python', command, cat_file], stdout = out_f) )
                #os.system('python kariba_calib.py ' + cat_file)
            for p in processes:
                p.wait()
            procs_available = num_procs
            cat_files = []
            cat_out_files = []
    
    

if __name__ == '__main__':
    
    sweep_name = cc_penalty_model + '_climate_categories_pilot_categories'
    
    
    
    sweep =  {              
                
                "Gwembe_pilot":{
                              'Gwembe_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Gwembe_1_node_MSAT_0.7_w_pilot_calib_a84a10ba-ef48-e511-93f8-f0921c16849c.json'],
                              'Gwembe_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Gwembe_1_node_MSAT_0.55_w_pilot_calib_bc4d78c1-ee48-e511-93f8-f0921c16849c.json'],
                              'Gwembe_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Gwembe_1_node_MSAT_0.35_w_pilot_calib_330dee5a-ee48-e511-93f8-f0921c16849c.json']
                              },

                "Munumbwe_no_pilot":{
                                   'Munumbwe_1_node_MSAT_0.7_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.7_w_NO_pilot_calib_5706b2fd-ed48-e511-93f8-f0921c16849c.json'],
                                   'Munumbwe_1_node_MSAT_0.55_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.55_w_NO_pilot_calib_d16f43a0-ed48-e511-93f8-f0921c16849c.json'],
                                   'Munumbwe_1_node_MSAT_0.35_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.35_w_NO_pilot_calib_55816a27-ed48-e511-93f8-f0921c16849c.json']
                                   },
                
                "Munumbwe_pilot":{
                                'Munumbwe_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.7_w_pilot_calib_0e55abc7-ec48-e511-93f8-f0921c16849c.json'],
                                'Munumbwe_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.55_w_pilot_calib_b6195a52-ec48-e511-93f8-f0921c16849c.json'],
                                'Munumbwe_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.35_w_pilot_calib_af4a9ce1-eb48-e511-93f8-f0921c16849c.json']
                                },
               
                "Lukonde_pilot":{
                               'Lukonde_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Lukonde_1_node_MSAT_0.7_w_pilot_calib_b324f46a-eb48-e511-93f8-f0921c16849c.json'],
                               'Lukonde_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Lukonde_1_node_MSAT_0.55_w_pilot_calib_44b497de-ea48-e511-93f8-f0921c16849c.json'],
                               'Lukonde_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Lukonde_1_node_MSAT_0.35_w_pilot_calib_7072677b-ea48-e511-93f8-f0921c16849c.json']                           
                               },
              
               "Sinamalima_no_pilot":{
                                     'Sinamalima_1_node_MSAT_0.7_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.7_w_NO_pilot_calib_d1297ec8-e848-e511-93f8-f0921c16849c.json'],
                                     'Sinamalima_1_node_MSAT_0.55_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.55_w_NO_pilot_calib_71ff6f4c-e848-e511-93f8-f0921c16849c.json'],
                                     'Sinamalima_1_node_MSAT_0.35_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.35_w_NO_pilot_calib_98ee99d3-e748-e511-93f8-f0921c16849c.json']
                                     },
          
               'Sinamalima_pilot':{
                                  'Sinamalima_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.55_w_pilot_calib_90937da8-e648-e511-93f8-f0921c16849c.json'],
                                  'Sinamalima_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.7_w_pilot_calib_8beb772d-e748-e511-93f8-f0921c16849c.json'],
                                  'Sinamalima_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_nodeMSAT_0.35_w_pilot_calib_a433fc45-e648-e511-93f8-f0921c16849c.json']
                                  }
              }
    
    
    '''
    sweep =  {              
           "Sinamalima_no_pilot":{
                                     'Sinamalima_1_node_MSAT_0.7_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.7_w_NO_pilot_calib_d1297ec8-e848-e511-93f8-f0921c16849c.json'],
                                     'Sinamalima_1_node_MSAT_0.55_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.55_w_NO_pilot_calib_71ff6f4c-e848-e511-93f8-f0921c16849c.json'],
                                     'Sinamalima_1_node_MSAT_0.35_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.35_w_NO_pilot_calib_98ee99d3-e748-e511-93f8-f0921c16849c.json']
                                     }
              }
    '''
    
    
    root_sweep_dir = os.path.join(sim_data_dir, sweep_name)
    if not os.path.exists(root_sweep_dir):
        os.mkdir(root_sweep_dir)
        
    viz_root_sweep_dir = os.path.join(kariba_viz_dir, sweep_name)
    if not os.path.exists(viz_root_sweep_dir):
        os.mkdir(viz_root_sweep_dir)
    
    
    multi_proc_run(sweep_name, sweep, 'kariba_calib.py')
    
    
    best_fits = {}
    all_fits = {}
    fit_terms = {}
    residuals = {}
    # combine all categories best fits in a single file in the sweep dir folder for spatial sims
    # also, create all error surface, prevalence time_series plots if any in the root sweep dir
    print "Best fits per category found."

    print "Combining per category best fits, all_fits, residuals; plotting best fits and residuals, and preparing d3js visualization."
    min_residual = float('inf')
    max_residual = 0.0
        
    for category in sweep:
        
        sweep_dir = os.path.join(sim_data_dir,category)
        
        with open(os.path.join(sweep_dir, best_fits_file), 'r') as best_fits_f:
            best_fits_cat = json.load(best_fits_f)
            
        with open(os.path.join(sweep_dir, all_fits_file), 'r') as all_fits_f:
            all_fits_cat = json.load(all_fits_f)
            
        with open(os.path.join(sweep_dir, fit_terms_file), 'r') as fit_terms_f:
            fit_terms_cat = json.load(fit_terms_f)
            
        with open(os.path.join(sweep_dir, residuals_file), 'r') as res_f:
            res_cat = json.load(res_f)
            
        if res_cat['min_residual'] < min_residual:
            min_residual = res_cat['min_residual']
        if res_cat['max_residual'] > max_residual:
            max_residual = res_cat['max_residual'] 
            
        residuals['min'] = min_residual
        residuals['max'] = max_residual
        
        print category + ' '  + str(len(best_fits_cat))    
       
        best_fits.update(best_fits_cat)
        all_fits.update(all_fits_cat)
        fit_terms.update(fit_terms_cat)
        
        print "Best fits updated for category " + category
    
    print "Processed " + str(len(best_fits)) + " clusters!"

    with open(os.path.join(root_sweep_dir, best_fits_file), 'w') as best_fits_f:
        json.dump(best_fits, best_fits_f, indent = 4)
        
    print "Stored best fit parameters json file to " + best_fits_file + " in " + root_sweep_dir
    
    with open(os.path.join(root_sweep_dir, all_fits_file), 'w') as all_fits_f:
        json.dump(all_fits, all_fits_f, indent = 4)
        
    with open(os.path.join(sim_data_dir, fit_terms_file), 'w') as fit_terms_f:
        json.dump(fit_terms, fit_terms_f, indent = 5)
        
    with open(os.path.join(root_sweep_dir, residuals_file), 'w') as res_f:
        json.dump(residuals, res_f, indent = 2)
   

        
    multi_proc_run(sweep_name, sweep, 'kariba_plots.py')
    
    
    print "Generating gazetteer"
    
    # combining tags reports from sweep categories
    
    
    sweep_dirs = []
    for category in sweep:
        sweep_dir = os.path.join(sim_data_dir, category)
        sweep_dirs.append(sweep_dir)
        
    combine_tags_reports(sweep_dirs, root_sweep_dir)
    
        
    with open(os.path.join(root_sweep_dir, best_fits_file), 'r') as best_fits_f:
        best_fits = json.load(best_fits_f)
    
    viz_conf = VizConfig(best_fits, sweep_name)
    viz_conf.update_d3js()
    viz_conf.generate_gazetteer()
    viz_conf.generate_gazetteer_header()
    viz_conf.generate_gazetteer_map()
    
    
    print "Gazetteer generated. Index file stored in " + kariba_viz_dir
    