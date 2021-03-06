import json
import math
import pandas as pd
from datetime import datetime,timedelta, date

from reference_data_scripts.nodenames import *

campaign_days_offset = [6*365 + 157, 6*365 + 157 + 60, 6*365 + 157 + 60 + 60, 7*365 + 157, 7*365 + 157 + 60, 7*365 + 157 + 60 + 60]
campaign_days = [6*365 + 160, 6*365 + 160 + 60, 6*365 + 160 + 60 + 60, 7*365 + 160, 7*365 + 160 + 60, 7*365 + 160 + 60 + 60]

pilot_round = 5*365 + 344
pilot_round_offset = 5*365 + 340

max_ref_reinfection_points = 4
max_ref_prevalence_points = 6 

channels = ['New Clinical Cases', 'New Diagnostic Prevalence']

objectives_channel_codes = ['prevalence']

channels_2_code = {
                   'New Clinical Cases': 'cc', 
                   'New Diagnostic Prevalence':'prevalence'
                   } 

channels_sample_points = {
                         'prevalence': [pilot_round_offset] + campaign_days_offset
                         } 

reports_channels = ['reinfections']

debug_flag = True
verbosity_flag = False
warnings_flag = True

calib_node_pop= 1000 

num_procs = 2

cc_weight = 0.0051
#cc_weight = 100001
reinf_weight = 0

#cc_penalty_model = 'trunc_ls_norm'
cc_penalty_model = 'ls_folded_norm_pop_cc_2011_2013_cc_w_' + str(cc_weight)

cc_agg_period = 6 # weeks
cc_agg_fold = True
cc_num_fold_bins = int(math.ceil(365/(cc_agg_period*7 + 0.0)))  
fold_start_date = pd.to_datetime('1/1/2011')
fold_end_date = pd.to_datetime('2/9/2014')



cc_correction_factor = 0.33

# rdt threshold above which the dtk sim cannot reach prevalence using current detection threshold and hrp2 model
rdt_max = 0.52

traces_plots_dir = 'prev_traces'
traces_base_file_name = 'prev_trace_'

cc_traces_plots_dir = 'cc_traces'
cc_traces_base_file_name = 'cc_trace_'

err_surfaces_plots_dir = 'err_surfaces'
err_surfaces_base_file_name = 'surf_'

weighted_cc_traces_plots_dir = 'weighted_cc_traces'
weighted_cc_traces_base_file_name = 'weighted_cc_trace_'



cc_sim_start_date = datetime(2005,1,1)
cc_ref_start_date = datetime(2010,5,24)
#cc_ref_end_date = datetime(2015,7,20)
#cc_ref_end_date = datetime(2013,12,30)
cc_ref_end_date = datetime(2014,2,9)



ref_data_dir = 'reference_data_scripts'
sim_data_dir = 'sim_data'
fit_terms_file = 'fit_terms.json'
residuals_file = 'residuals.json' 
calibration_data_file = 'calibration_cc.json'
tags_data_file = 'tags.json'
tags_report_data_file = 'tags_report.json'
best_fits_file = cc_penalty_model + '_best_fits.json'
all_fits_file =  cc_penalty_model + '_all_fits.json'


# fit terms settings
load_cc_penalty = True
load_prevalence_mse = True
load_reinf_penalty = False

scale_fit_terms = False # scale individual fit function components to 0-1 range
load_scaled_fit_terms = False # recompute min/max values for individual fit function components for 0-1 scaling; set to True if those have already been computed ('min_terms/max_terms' are present in fit_terms.json for each cluster)

load_prevalence_based_ccs = True

#use_scaled_fit_terms = False

# add fit terms here (e.g. penalty types) as objective function options change
# reflect the fit_terms.json schema
fit_terms_types = {
                   'corr_folded':['cc_penalty','corr_folded','penalty'],
                   'ls_norm':['cc_penalty', 'ls_norm'],
                   'mse':['mse']
                  }


# LIKELIHOOD SAMPLING SETTINGS
weighted_ccs_by_hfca_id_file = 'weighted_ccs_by_hfca.json'
sample_size_percentile = 5 # sample from the top 5% points
sample_size = 40
#sample_size_ccs_per_cluster = 2
min_num_combos_per_hfca = 1000


# VISUALIZATION SETTINGS
subopt_plots_threshold = 0.1 # only plot traces, surfaces that are suboptimal (when required) if the corresponding fit is within some fraction of the optimal (e.g. 0.1)

err_surface_types = {'fit':{'title':'Clinical cases + Prevalence', 'marker':'*'}, 'cc_penalty':{'title':'Clinical cases', 'marker':'s'}, 'mse':{'title':'Prevalence', 'marker':'o'}}

root_viz_dir = 'visualization'
kariba_viz_dir = os.path.join(root_viz_dir, 'kariba_viz')
d3js_src_dir = 'C:\\Users\\Mnikolov\\workspace\\KaribaCalibrationDashboard'

d3js_src_files = [
                    'topolakes.json',
                    'hfcas.json',
                    'hhs.json',
                    'hhs_shapes.json',
                    'cluster_pop_ts.tsv',
                    'clusters_hulls.json',
                    'navigation.js',
                    'require.js',
                    'buttons.js',
                    'colorbar.js',
                    'colorbrewer2.js',
                    'dtk_multinode.css',
                    'dtk_multinode.js',
                    'index.html',
                    'zambia_flag.png'
                ]



gazetteer_params_file_name = 'gazetteer.json'
gazetteer_params_header_file_name = 'gazetteer_header.json'
gazetteer_base_map_file_name = 'map.json'
gazetteer_sim_mn_base_map_file_name = cc_penalty_model + '_mn_map.json'
 
markers = ['o', 'v', '^', '<', '>', '1', '2', '3', '4', '8', 's', 'p', '*', 'h', 'H', '+', 'x', 'd', '|', '_']

opt_marker = 'D'
opt_marker_size = 6  


def get_fold_bins(): 
    
    fold_bins = {}
    for i in range(0, cc_num_fold_bins):
        #fold_bins[i] = {'dates':[], 'cases_model':'nan', 'cases_ref':'nan', 'cases_ref_entries':[]}
        fold_bins[i] = {'cases_model':'nan', 'cases_ref':'nan'}

    '''    
    with open(os.path.join(ref_data_dir,'fold_bins.json'), 'r') as folds_f:
        fold_bins = json.load(folds_f) 
    '''
    return fold_bins
        
def get_hfca_ids_true():
    return health_facility_names.keys()

def get_hfca_ids():
    #return ['80206']
    return health_facility_names.keys()


def hfca_id_2_cluster_ids(hfca_id):
    print(hfca_id)
    return hfca_id_cluster_ids[int(hfca_id)]

def hfca_id_2_facility(hfca_id):
    return health_facility_names[int(hfca_id)]

def get_cc_cluster_weight_factor(cluster_id):
    
    with open(os.path.join(ref_data_dir,'clusters_hs.json'),'r') as cl_f:
        clusters_hs_weights = json.load(cl_f)
        
    if cluster_id in clusters_hs_weights:
        return clusters_hs_weights[cluster_id]
    else:
        # return a generic average rate of health seeking if the cluster is not present in the data
        # issue a warning 
        print 'Cluster ' + cluster_id + ' was not found in cluster health seeking data. Using generic health seeking rate of ' + str(cc_correction_factor) 
        return cc_correction_factor

# get all clusters mapped to a given category
def category_2_clusters(category):
    
    # get category to subsets of clusters map; e.g. categories Sinamalima_pilot, Sinamalima_no_pilot, Gwembe_pilot, etc.
    with open(os.path.join(ref_data_dir,'cluster_categories.json'), 'r') as cats_f:
        categories_2_clusters = json.load(cats_f) 
    
    if category in categories_2_clusters:
        return categories_2_clusters[category]
    else:
        return None

def get_cluster_category(cluster_id):
    # get category to subsets of clusters map; e.g. categories Sinamalima_pilot, Sinamalima_no_pilot, Gwembe_pilot, etc.
    with open(os.path.join(ref_data_dir,'cluster_categories.json'), 'r') as cats_f:
        categories_2_clusters = json.load(cats_f)
        
    for category, cluster_ids in categories_2_clusters.iteritems():
        if cluster_id in cluster_ids:
            return category
    # if category is not found for the cluster return Sinamalima_no_pilot category and issue a warning
    print 'Cluster ' + cluster_id + ' category was not found. Assigning to category Sinamalima_no_pilot.'
    return 'Sinamalima_no_pilot'
    
# get prevalences and pops across different rounds for all clusters; NOTE: a value of -1000 indicates a missing entry
with open(os.path.join(ref_data_dir,'clusters_prevalences.json'), 'r') as clusters_prevs_f:
        clusters_2_prevs = json.load(clusters_prevs_f)
     
def cluster_2_prevs(cluster_id):    
    if cluster_id in clusters_2_prevs:
        pilot_prev = cluster_2_pilot_prev(cluster_id)
        if pilot_prev:
            clusters_2_prevs[cluster_id]['RDT']['0'] = pilot_prev
        else:
            clusters_2_prevs[cluster_id]['RDT']['0'] = -1000     
        return clusters_2_prevs[cluster_id]['RDT']
    else:
        return None
    
def cluster_2_pops(cluster_id):
    if cluster_id in clusters_2_prevs:
        return clusters_2_prevs[cluster_id]['pop']
    else:
        return None
    
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False   
    
def cluster_2_mean_pop(cluster_id):
    if cluster_id in clusters_2_prevs:
        pops =  clusters_2_prevs[cluster_id]['pop']
    
        mean_pop = 0.0
        count_pops = 0
        for pop in pops.values():
            if is_number(pop) and pop > 0:
                mean_pop = mean_pop + pop
                count_pops = count_pops + 1
        
        mean_pop = mean_pop / (count_pops + 0.0)
        
        return mean_pop
                
    else:
        return None


def hfca_2_pop(hfca_id):
    with open(os.path.join(ref_data_dir, 'hfcas_pops.json'), 'r') as hfcas_pops_f:
        hfcas_pops = json.load(hfcas_pops_f)
    if hfca_id in hfcas_pops:
        return hfcas_pops[hfca_id]
    else:
        return None


def cluster_2_pilot_prev(cluster_id):
    # get pilot round prevalences for all clusters
    with open(os.path.join(ref_data_dir,'clusters_pilot_prevalences.json'), 'r') as clusters_pilot_prevs_f:
        clusters_2_pilot_prevs = json.load(clusters_pilot_prevs_f)
    if cluster_id in clusters_2_pilot_prevs:
        return clusters_2_pilot_prevs[cluster_id]
    else:
        return None


def cluster_2_itn_traj(cluster_id):
    with open(os.path.join(ref_data_dir, 'itn_prior.json'), 'r') as itn_f:
        clusters_2_itn_trajs = json.load(itn_f)
    if cluster_id in clusters_2_itn_trajs:
        return clusters_2_itn_trajs[cluster_id]
    else:
        return None


def cluster_2_drug_cov(cluster_id):
    with open(os.path.join(ref_data_dir, 'drug_coverages.json'), 'r') as drug_cov_f:
        clusters_2_drug_covs = json.load(drug_cov_f)
    if cluster_id in clusters_2_drug_covs:
        return clusters_2_drug_covs[cluster_id]
    else:
        return None

def cluster_2_reinfection_rates(cluster_id):
    # get reinfection rate per pairs of rounds for each cluster
    with open(os.path.join(ref_data_dir,'reinfections.json'), 'r') as clusters_reinf_f:
        clusters_2_reinfections = json.load(clusters_reinf_f)
    if cluster_id in clusters_2_reinfections:
        return clusters_2_reinfections[cluster_id]
    else:
        return None    
    
    
def cluster_2_cc(cluster_id):
    # get reinfection rate per pairs of rounds for each cluster
    with open(os.path.join(ref_data_dir,'cc.json'), 'r') as clusters_cc_f:
        clusters_2_cc = json.load(clusters_cc_f)
    hfca_id = cluster_id.split('_')[0]
    if hfca_id in clusters_2_cc:
        return clusters_2_cc[hfca_id]
    else:
        return None    