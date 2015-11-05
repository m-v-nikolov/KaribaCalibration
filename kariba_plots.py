import os
import sys
import json
import gc

from comps_2_sim_data import get_sweep_results, combine_sweep_results
from sim_data_2_models import calib_data_2_models_list

from kariba_settings import sim_data_dir, best_fits_file, all_fits_file, calibration_data_file, traces_plots_dir, traces_base_file_name, cc_traces_plots_dir, cc_traces_base_file_name, err_surfaces_plots_dir, err_surfaces_base_file_name, cc_penalty_model, kariba_viz_dir, residuals_file, err_surface_types
from kariba_fit import KaribaFit
from plot_utils import PlotUtils


if __name__ == '__main__':
        
    print "Processing " + sys.argv[1]

    with open(sys.argv[1], 'r') as cat_f:
        cat_record = json.load(cat_f)
    
    # only do once per calibration sweep
    #prepare_calib_files(cat_record)
    
    
    # once calib files are stored just load them
    sweep_dir =  cat_record['sweep_dir']
    root_sweep_dir = cat_record['root_sweep_dir']
    viz_root_sweep_dir = cat_record['viz_root_sweep_dir']
    
    with open(os.path.join(root_sweep_dir, residuals_file),'r') as res_f:
        residuals = json.load(res_f)   
     
     
    # slight chance of race conditions per current design; to avoid, sync processes in the future
     
    for category, experiments in cat_record['categories'].iteritems():

        with open(os.path.join(sweep_dir, best_fits_file), 'r') as best_fits_f:
            best_fits_cat = json.load(best_fits_f)
            
        with open(os.path.join(sweep_dir, all_fits_file), 'r') as all_fits_f:
            all_fits_cat = json.load(all_fits_f)
        
        with open(os.path.join(sweep_dir, calibration_data_file), 'r') as calib_data_f:
            calib_data_cat = json.load(calib_data_f)
         
        # plot clinical cases
        pp = PlotUtils(best_fits_cat, all_fits_cat, calib_data_cat, residuals, viz_root_sweep_dir, category)
        
        
        #create respective figure dirs if not existing;         
        
        
        # create clinical cases plots directory if it doesn't exist
        cc_plots_dir = os.path.join(viz_root_sweep_dir, cc_traces_plots_dir)
        if not os.path.exists(cc_plots_dir):
            print "Creating clinical incidence directory " + str(cc_plots_dir)
            os.mkdir(cc_plots_dir)
        
        print "Plotting per cluster clinical cases for category " + category
        print "Plots stored in "  + cc_plots_dir
        
        pp.plot_calib_cc_traces_clusters()
        
        print "Plotting clinical cases done for category " + category
        
        
        # create err surfaces plots directory if it doesn't exist
        err_plots_dir = os.path.join(viz_root_sweep_dir, err_surfaces_plots_dir)
        if not os.path.exists(err_plots_dir):
            print "Creating error surfaces directory " + str(err_plots_dir)
            os.mkdir(err_plots_dir)
        
        print "Plotting per cluster error surfaces for category " + category
        
        pp.plot_calib_err_surfaces(err_surface_types)
        
        print "Plotting error surfaces done for category " + category
        print "Plots stored in "  + err_plots_dir
        
        
        # create prevalence plots directory if it doesn't exist
        prev_plots_dir = os.path.join(viz_root_sweep_dir, traces_plots_dir)
        if not os.path.exists(prev_plots_dir):
            print "Creating prevalence directory " + str(prev_plots_dir)
            os.mkdir(prev_plots_dir)
        
        print "Plotting per cluster prevalence for category " + category
        
        pp.plot_calib_prev_traces()
        
        print "Plotting prevalence traces done for category " + category 
        print "Plots stored in "  + prev_plots_dir
            