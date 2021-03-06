import os
import sys
import json
import gc

from comps_2_sim_data import get_sweep_results, combine_sweep_results
from sim_data_2_models import calib_data_2_models_list

from kariba_settings import  load_cc_penalty, load_prevalence_mse, load_reinf_penalty, sim_data_dir, calibration_data_file, tags_data_file, objectives_channel_codes, reports_channels, channels, best_fits_file, all_fits_file, residuals_file, fit_terms_file, cc_penalty_model, scale_fit_terms, load_scaled_fit_terms
from kariba_utils import load_fit_terms, normalize_fit_terms
from kariba_fit import KaribaFit


def prepare_calib_files(cat_record):

    sweep_dir = cat_record['sweep_dir']
      
    for category, experiments in cat_record.iteritems():
    
        for base_dir, sim_meta_config_files in experiments.iteritems():
            
            calib_file_path = os.path.join(sim_data_dir, sweep_dir, base_dir, calibration_data_file)
            tags_data_file_path = os.path.join(sim_data_dir, sweep_dir, base_dir, tags_data_file)
            
            get_sweep_results(sim_meta_config_files, calib_file_path, tags_data_file_path)
            
        comb_base_dirs = [os.path.join(sim_data_dir, base_dir) for base_dir in experiments.keys()] 
        comb_calib_dir = os.path.join(sim_data_dir, category) 
    
        combine_sweep_results(comb_base_dirs, sweep_dir)


def calibrate(category, sweep_dir):
    
    calib_file_path = os.path.join(sweep_dir, calibration_data_file)
    
    with open(calib_file_path, 'r') as calib_f:
        calib_data = json.load(calib_f)
        
    fit_terms_file_path = os.path.join(sweep_dir, fit_terms_file)
    
    if scale_fit_terms and not load_scaled_fit_terms:
        fit_terms = normalize_fit_terms(fit_terms_file_path)
    else:
        fit_terms = load_fit_terms(fit_terms_file_path)
        
    kariba_fit_cat = KaribaFit(category, calib_data, fit_terms = fit_terms)

    best_fits, all_fits = kariba_fit_cat.fit()
    
    del calib_data
    gc.collect()
    
    
    # record best fit parameters
    best_fits_f_path = os.path.join(sweep_dir, best_fits_file)
    with open(best_fits_f_path, 'w') as best_fits_f:
        json.dump(best_fits, best_fits_f, indent = 4)
        
        
    # go through all fit models, extract params and individual fit function terms for each
    fit_entries = {}
    for cluster_id, models in all_fits['models'].iteritems():
        fit_entries[cluster_id] = {}
        for model in models:
            fit_entries[cluster_id].update(model.fit_entry())
        
        # making sure we record optimal values from previous computation back in fit_terms 
        if 'min_terms' in fit_terms[cluster_id]:
            fit_entries[cluster_id]['min_terms'] = fit_terms[cluster_id]['min_terms']
        
        if 'max_terms' in fit_terms[cluster_id]:
            fit_entries[cluster_id]['max_terms'] = fit_terms[cluster_id]['max_terms'] 
    
    
    all_fits_f_path = os.path.join(sweep_dir, all_fits_file)
    with open(all_fits_f_path, 'w') as all_fits_f:
        json.dump(fit_entries, all_fits_f, indent = 5)
        
    with open(fit_terms_file_path, 'w') as fit_terms_f:
        json.dump(fit_entries, fit_terms_f, indent = 5)
        
    residuals = {'min_residual':all_fits['min_residual'], 'max_residual':all_fits['max_residual']}
    
    residuals_f_path = os.path.join(sweep_dir, residuals_file)
    with open(residuals_f_path, 'w') as res_f:
        json.dump(residuals, res_f, indent = 2)
    
    return best_fits, all_fits, residuals
    


if __name__ == '__main__':
        
    print "Processing " + sys.argv[1]

    with open(sys.argv[1], 'r') as cat_f:
        cat_record = json.load(cat_f)
    
    # only do once per calibration sweep
    #prepare_calib_files(cat_record)
    
    
    # once calib files are stored just load them
    sweep_dir =  cat_record['sweep_dir']   
     
    for category, experiments in cat_record['categories'].iteritems():
                            
        best_fits = None
        all_fits = None
        
        # only calibrate if needed otherwise comment out 
        best_fits, all_fits, residuals = calibrate(category, sweep_dir)
        
        
        # generate plots
        
        '''
        # if no calibration results are available try to load them
        if not (best_fits and all_fits):
            
            best_fits_f_path = os.path.join(sweep_dir, best_fits_file)
            with open(best_fits_f_path, 'r') as best_fits_f:
                best_fits = json.load(best_fits_f)
            

            all_fits_f_path = os.path.join(sweep_dir, all_fits_file)
            with open(all_fits_f_path, 'r') as all_fits_f:
                all_fits = json.load(all_fits_f)
        
        '''
        
        
            
            
    print "Config file for calibration category " + cat_record.keys()[0] + " has been successfully loaded!"
    print "All output files will be saved to " + sweep_dir