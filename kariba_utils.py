import ast
import os
import json
import pandas as pd
import numpy as np
import random
import itertools as it
import math
from scipy.stats import scoreatpercentile as satp

from datetime import datetime,timedelta

from kariba_settings import hfca_2_pop, cc_correction_factor, get_fold_bins, cc_agg_fold, cc_agg_period, fold_start_date, fold_end_date, calib_node_pop, cluster_2_cc, cc_sim_start_date, cc_ref_start_date, cc_ref_end_date, tags_report_data_file, fit_terms_file, sim_data_dir, cc_penalty_model, fit_terms_types, hfca_id_2_cluster_ids, get_hfca_ids, get_cc_cluster_weight_factor, cluster_2_pops, sample_size_percentile, weighted_ccs_by_hfca_id_file, min_num_combos_per_hfca, load_prevalence_based_ccs, sample_size, cluster_2_mean_pop, get_cluster_category

from utils import warn_p, debug_p, feature_scale
from copy import deepcopy


def get_cc_penalty(fit_entry):
    
    if 'corr_folded' in cc_penalty_model:
        return fit_entry['fit_terms']['cc_penalty']['corr_folded']['penalty']
    
    if 'corr_not_folded' in cc_penalty_model:
        return fit_entry['fit_terms']['cc_penalty']['corr_not_folded']['penalty']
    
    if 'ls_folded_norm' in cc_penalty_model: 
        return fit_entry['fit_terms']['cc_penalty']['ls_norm']
    elif 'ls_norm_not_folded' in cc_penalty_model:
        return fit_entry['fit_terms']['cc_penalty']['ls_norm_not_folded']
                    
    if 'ls_no_norm' in cc_penalty_model: 
        return fit_entry['fit_terms']['cc_penalty']['ls_no_norm']

    debug_p('No clinical cases penalty found. This should not happen!')
    return None
        

def sim_key_2_group_key(sim_key):
    group_key_terms = sim_key.split('_')[2:]
    group_key = ''
    for i,gt in enumerate(group_key_terms):
        if i < len(group_key_terms)-1:
            group_key += gt + '_'
        else:
            group_key += gt
    return group_key


def get_sim_key(temp_h, const_h, itn_level, drug_level):
    return str(temp_h) + '_' + str(const_h) + '_' + get_sim_group_key(itn_level, drug_level)
    
    
def get_sim_group_key(itn_level, drug_level):
    return str(itn_level) + '_' + str(drug_level)


def sim_meta_2_itn_level(meta):
    itn_level_struct = ast.literal_eval(meta['add_ITN_mult'])
    itn_level = itn_level_struct[0][1][0][0][1]
    return itn_level
    
    
def sim_meta_2_drug_cov(meta):
    drug_coverage_level_struct = ast.literal_eval(meta['add_drug_multi_campaigns'])
    drug_coverage_level = drug_coverage_level_struct[0][1][0]['coverage']
    return drug_coverage_level
        
    
def sim_meta_2_temp_h(meta):
    x_temp_h = float(meta['x_Temporary_Larval_Habitat'])
    return x_temp_h
    
    
def sim_meta_2_const_h(meta):    
    const_h_struct = ast.literal_eval(meta['scale_larval_habitats_single'])
    const_h = float(const_h_struct[0][1][1])
    return const_h


def get_model_params(model):

    model_meta_data = model.get_meta()
    sim_meta_data = model_meta_data['sim_meta']
        
    # get best fit params
    temp_h = sim_meta_2_temp_h(sim_meta_data)
    const_h = sim_meta_2_const_h(sim_meta_data)
    itn_level = sim_meta_2_itn_level(sim_meta_data)
    drug_coverage_level = sim_meta_2_drug_cov(sim_meta_data)
    
    return temp_h, const_h, itn_level, drug_coverage_level


def combine_tags_reports(base_dirs, output_dir):

    tags_report_comb = {}
    for base_dir in base_dirs:
        print "Processing tags report in " + base_dir
                    
        # update tags information
        with open(os.path.join(base_dir, tags_report_data_file), 'r') as tags_report_f:
            tags_report = json.load(tags_report_f)
            
        for param, values in tags_report.iteritems():
            if param not in tags_report_comb:
                tags_report_comb[param] = []
            tags_report_comb[param] = set(tags_report_comb[param])
            tags_report_comb[param].update(values)
            tags_report_comb[param] = list(tags_report_comb[param])
            
        print "DONE"
        print 
        
    if not os.path.exists(output_dir):
         os.makedirs(output_dir)
        
    print "Writing tags to " + os.path.join(output_dir, tags_report_data_file)
    with open(os.path.join(output_dir, tags_report_data_file), 'w') as tags_report_f:
        json.dump(tags_report_comb, tags_report_f)            
    print "DONE" 


def cc_data_aggregate(model_clinical_cases, cluster_id, cc_cluster_weight_factor = None):
    
    # aggregate on a periodic basis to as many values as there are in the ref data
    ccs_ref_agg = cluster_2_cc(cluster_id)                                        
    ccs_model_agg = []
    ccs_ref_agg_cleaned = []
    
    hfca_id = cluster_id.split('_')[0]
    
    hfca_pop = hfca_2_pop(hfca_id)
    
    pop_norm_factor = 0.0
    
    if not cc_cluster_weight_factor:
        pop_norm_factor = cc_correction_factor*(hfca_pop + 0.0)/calib_node_pop
    else:
        pop_norm_factor = cc_cluster_weight_factor
        
    #pop_norm_factor = 1
    #debug_p('pop of hfca ' + hfca_id + ' is ' + str(hfca_pop))
    #debug_p('pop norm factor for cluster ' + cluster_id + ' is ' + str(pop_norm_factor))
    dates, cases = zip(*ccs_ref_agg)
    
    
    sim_start_date = cc_sim_start_date
    '''
    ref_start_date = cc_ref_start_date
    ref_end_date = cc_ref_end_date
    '''

    #ref_start_date = max(min(dates),cc_ref_start_date) 
    #ref_end_date = min(max(dates) ,cc_ref_end_date)
    ref_start_date = pd.to_datetime(min(dates)) 
    ref_end_date = min(pd.to_datetime(max(dates)), cc_ref_end_date)
    
    # note: assume the simulation has started more than 6 weeks before clinical cases collection;
    # this should always be the case for a well tempered simulation
    sim_ref_start_date = ref_start_date - timedelta(days = 6*7 - 1) 
         
    fold_bins = get_fold_bins()
          
    model_cc_start_idx = ((sim_ref_start_date - sim_start_date)/ np.timedelta64(1, 'D')).astype(int) - 1 # index in model clinical cases at sim_ref_start_date
    model_cc_end_idx = ((pd.to_datetime(ref_end_date) - pd.to_datetime(sim_start_date))/ np.timedelta64(1, 'D')).astype(int) # index in model clinical cases until ref_end_date
    
    cur_date = sim_ref_start_date
    
    ccs_model = []
    for i,idx in enumerate(range(model_cc_start_idx, model_cc_end_idx)):
        if i > 0:
            cur_date = cur_date+timedelta(days = 1) 

        ccs_model.append((cur_date, model_clinical_cases[idx]))   
 
    cases_period = 0.0
    periods = 0
    
    #debug_p('==============================================')
    #debug_p(cluster_id)
    #debug_p('before adding cases ' + str(fold_bins))
    
    '''
    with open('fold_bins_before_cases_'+ cluster_id+'.json','w') as d_ccs_f:
        json.dump(fold_bins, d_ccs_f, indent = 4)
    '''
    
    for i, (date, cases) in enumerate(ccs_model):
        cases_period = cases_period + cases
        if (i+1) % (6*7) == 0 or (i+1) == len(ccs_model):
            if periods < len(ccs_ref_agg):
                ccs_ref_agg[periods][0] = str(ccs_ref_agg[periods][0]) 
                if not ccs_ref_agg[periods][1] == 'nan': 
                    ccs_model_agg.append((str(date), cases_period*pop_norm_factor)) #timedelta operation ensures we record the end of the period as the clinical cases report date
                    fold_bins = update_fold_bins(date, cases_period*pop_norm_factor, ccs_ref_agg[periods][1], fold_bins)
                else:
                    ccs_model_agg.append((str(date), 'nan'))
                
                cases_period = 0.0
                periods = periods + 1
            else:
                break
    
    #debug_p('after adding cases ' + str(fold_bins))
    
    '''
    with open('fold_bins_after_cases_'+ cluster_id+'.json','w') as d_ccs_f:
        json.dump(fold_bins, d_ccs_f, indent = 4)
    
    with open('debug_ccs_model_'+ cluster_id+'.json','w') as d_ccs_f:
        json.dump(ccs_model_agg, d_ccs_f, indent = 3)
        
    with open('debug_ccs_ref_'+ cluster_id+'.json','w') as d_ccs_f:
        json.dump(ccs_ref_agg, d_ccs_f, indent = 3)
    '''
      
    return ccs_model_agg, ccs_ref_agg, fold_bins



def update_fold_bins(date, cases_model, cases_ref, fold_bins):
    
    jan_1 = pd.to_datetime( '1/1/' + str(date.year) ) # get jan 1 of cur_date's year 
    num_days_in_fold = ((pd.to_datetime(date) - jan_1)  / np.timedelta64(1, 'D')).astype(int) # get num days in fold since jan 1 
    
    if fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_model'] =='nan':
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_model'] = cases_model
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['num_updates'] = 1
    else:
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_model'] = fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_model'] + cases_model
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['num_updates'] = fold_bins[num_days_in_fold/(cc_agg_period * 7)]['num_updates'] + 1 
    
    if fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] =='nan':
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] = cases_ref
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['num_updates'] = 1
    else:
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] = fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] + cases_ref
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['num_updates'] = fold_bins[num_days_in_fold/(cc_agg_period * 7)]['num_updates'] + 1
    
    
    '''
    found = False
    for idx,bin in fold_bins.iteritems():
        
        if date in bin['dates']:
            
            if bin['cases_model'] == 'nan':
               bin['cases_model'] = cases_model
            else:
               bin['cases_model'] = bin['cases_model'] + cases_model
            
            if bin['cases_ref'] == 'nan':
               bin['cases_ref'] = cases_ref
            else:
               bin['cases_ref'] = bin['cases_ref'] + cases_ref
            
            #bin['cases_ref_entries'].append((date, cases_ref))
            
            found = True
            
            break
        
    if not found: 
        debug_p('Date ' + str(date) + ' was not found in fold_bins. This should not happen!')
    '''
        
    return fold_bins


def fold_nan_clean(fold_bins):
    
    ccs_folded_model = []
    ccs_folded_ref = []
    
    for fold_idx in sorted(fold_bins): # sort to preserve the bins order in time
        if not (fold_bins[fold_idx]['cases_model'] == 'nan' or fold_bins[fold_idx]['cases_ref'] == 'nan'):
           ccs_folded_model.append(fold_bins[fold_idx]['cases_model']/(fold_bins[fold_idx]['num_updates'] + 0.0)) # average by the number of updates / years data has been aggregated per bin
           ccs_folded_ref.append(fold_bins[fold_idx]['cases_ref']/(fold_bins[fold_idx]['num_updates'] + 0.0)) # average by the number of updates / years data has been aggregated per bin
           
    return ccs_folded_model, ccs_folded_ref
            

def cc_data_nan_clean(ccs_model_agg, ccs_ref_agg, cluster_id):
    
    ccs_model_agg_clean = [] 
    ccs_ref_agg_clean = []
    count_nans = 0
    
    for i,(date, cases_model) in enumerate(ccs_model_agg):
        cases_ref = ccs_ref_agg[i][1]
        if not cases_ref == 'nan':
            ccs_model_agg_clean.append(cases_model)
            ccs_ref_agg_clean.append(cases_ref)
        else:
            count_nans = count_nans + 1

    return ccs_model_agg_clean, ccs_ref_agg_clean 


def get_cc_model_ref_traces(model_clinical_cases, cluster_id, cc_cluster_weight_factor = None):
    
    ccs_model_agg, ccs_ref_agg, fold_bins = cc_data_aggregate(model_clinical_cases, cluster_id, cc_cluster_weight_factor)
    
    if cc_agg_fold:
        ccs_model_agg, ccs_ref_agg = fold_nan_clean(fold_bins)
    else:
        ccs_model_agg, ccs_ref_agg = cc_data_nan_clean(ccs_model_agg, ccs_ref_agg, cluster_id)
        
    return ccs_model_agg, ccs_ref_agg


def load_fit_terms(fit_terms_file_path):
    
    if os.path.exists(fit_terms_file_path):
        with open(fit_terms_file_path, 'r') as ft_f:
            fit_terms = json.load(ft_f)   
    else:
        fit_terms = None
        
    return fit_terms


# supply fit terms type that containing the fit term hierarchy path in fit_terms.json
# e.g. clinical cases penalty computed using rank correlation on folded clinical cases timeseries 
# would have type ['cc_penalty', 'corr_folded', 'penalty']

def normalize_fit_terms(fit_terms_file_path):
    
    if os.path.exists(fit_terms_file_path):
        with open(fit_terms_file_path, 'r') as ft_f:
            fit_terms = json.load(ft_f)
    else:
        error_loading_fit_terms()
    
    terms = {}

    for fit_terms_type, path in fit_terms_types.iteritems():
        terms[fit_terms_type] = []
    
    count = 0
    for cluster_id, fit_entries in fit_terms.iteritems():
        
        #if not 'min_terms' in fit_terms[cluster_id]: 
        fit_terms[cluster_id]['min_terms'] = {}
                                    
        #if not 'max_terms' in fit_terms[cluster_id]: 
        fit_terms[cluster_id]['max_terms'] = {}

        
        for sim_key, fit_entry in fit_entries.iteritems():
            
            if not (sim_key == 'min_terms'  or sim_key == 'max_terms'):

                for fit_terms_type, path in fit_terms_types.iteritems():
                    fit_entry_copy = deepcopy(fit_entry)
                    path_copy = deepcopy(path)
                        
                    terms[fit_terms_type].append(unroll_term(fit_entry_copy['fit_terms'], path_copy))
        
        for fit_terms_type, path in fit_terms_types.iteritems():
            min_term = np.min(terms[fit_terms_type])
            max_term = np.max(terms[fit_terms_type])
            
            path_copy = deepcopy(path)
            
            if not fit_terms_type in fit_terms[cluster_id]['min_terms']:
                fit_terms[cluster_id]['min_terms'][fit_terms_type] = {}
                
            if not fit_terms_type in fit_terms[cluster_id]['max_terms']:
                fit_terms[cluster_id]['max_terms'][fit_terms_type] = {}
                            
            fit_terms[cluster_id]['min_terms'][fit_terms_type] = roll_term(fit_terms[cluster_id]['min_terms'][fit_terms_type], path_copy, min_term)
            fit_terms[cluster_id]['max_terms'][fit_terms_type] = roll_term(fit_terms[cluster_id]['max_terms'][fit_terms_type], path_copy, max_term)
        count = count + 1
    
    '''
    debug_p( fit_terms.keys() )
    debug_p( fit_terms['81102_4']['min_terms'] )
    debug_p( fit_terms['81102_4']['max_terms'] )
    '''
            
    with open(fit_terms_file_path, 'w') as ft_f:
        json.dump(fit_terms, ft_f)
        
    return fit_terms
            

def unroll_term(fit_entry, path):
    if len(path) == 1: # if only one level in hierarchy left, return
        return fit_entry[path[0]]
    else:
        fit_entry = fit_entry[path[0]]
        path = path[1:]
        return unroll_term(fit_entry, path)
         

# can be implemented better w/ recursion
def roll_term(fit_term, path, extremal_term):
    path_copy = deepcopy(path)
    
    fit_term_nl = fit_term
    
    idx = 0
    while path[0] in fit_term_nl:
        fit_term_nl = fit_term_nl[path[0]]
        
        if len(path[1:]) == 1 and path[1] in fit_term_nl: # check if only one entry is left in path and if that entry is in dictionary update it with the corresponding value 
            fit_term_nl[path[1]] = extremal_term
            return fit_term
        
        path = path[1:]
        idx = idx + 1
        
                          
    while idx < len(path_copy):
         
        if idx + 1 == len(path_copy):
            fit_term_nl.update({path_copy[idx]:extremal_term})
        else: 
            fit_term_nl.update({path_copy[idx]:{}})
        
        fit_term_nl = fit_term_nl[path_copy[idx]]
        idx = idx + 1 

    return fit_term

def gazetteer_list_update(data, file_path):
    
    with open(file_path, 'r') as f_p:
        cur_data = json.load(f_p)
    
    found = False
    for entry in cur_data:
        if entry['model'] == data['model']:
            entry['select'].append(data['select'][1]) # append only the new model selection (skipping the initial empty selection)
            found = True 
            break
        
    if not found: # if model is not already in the list append new model
        cur_data.append(data)
    
    with open(file_path, 'w') as f_p:
        json.dump(cur_data, f_p, indent = 4)
          
def error_loading_fit_terms():
    
    debug_p('Could not load fit terms! Check whether fit terms file at ' + os.path.join(sim_data_dir, fit_terms_file) + ' is accessible.')
    
    raise ValueError('Could not load fit terms! Check whether fit terms file at ' + os.path.join(sim_data_dir, fit_terms_file) + ' is accessible.')

def get_prevalence_opt_region_sims(best_fits, all_fits, cluster_id):
    
    opt_region_sims = []
    
    cluster_record = best_fits[cluster_id]
    
    opt_group_key = cluster_record['group_key']

    error_points = {}
    for sim_key,fit_entry in all_fits[cluster_id].iteritems():

        if sim_key == 'min_terms' or sim_key == 'max_terms':
            continue

        mse = fit_entry['fit_terms']['mse']
        group_key = fit_entry['group_key']    
        
            
        if group_key not in error_points:
            error_points[group_key] = {
                                           'mse':[],
                                           'sim_key': []
                                        }
            
        error_points[group_key]['mse'].append(mse)
        error_points[group_key]['sim_key'].append(sim_key)
               
                
    for j,group_key in enumerate(error_points.keys()):

        if group_key == opt_group_key:     
                z = error_points[group_key]['mse']
                sim_keys = error_points[group_key]['sim_key']
                sample_space = zip(sim_keys, z) 
                
                opt_region_sims  = [(group_key, sim_key) for (sim_key, res) in sample_space if res <= satp(z, sample_size_percentile)]
                break
                
    return opt_region_sims

     
    
def get_prevalence_based_cc(best_fits, all_fits, calib_data):
    
    weighted_ccs_by_hfca_id_file_path = os.path.join(sim_data_dir, weighted_ccs_by_hfca_id_file)
    
    hfca_ids = get_hfca_ids()
    
    debug_p('Getting clinical cases samples based on prevalence optimal regions')
    
    weighted_ccs_model_agg_by_hfca = {} 
    for hfca_id in hfca_ids:
        hfca_id = str(hfca_id)
        weighted_ccs_model_agg_by_hfca[hfca_id] = {}
        
        cluster_ids = hfca_id_2_cluster_ids(hfca_id)
        
        for cluster_id in cluster_ids:
            
            cluster_cat = get_cluster_category(cluster_id)
            
            sims_opt_region = get_prevalence_opt_region_sims(best_fits, all_fits, cluster_id)
            
            # assume sample size is always less than the size of the population!
            sample_sims_opt_region = random.sample(sims_opt_region, sample_size) 
            
            for i, (sample_group_key, sample_sim_key) in enumerate(sample_sims_opt_region): 
                      
                sample_cc_trace = calib_data[cluster_cat][sample_group_key][sample_sim_key]
                cc_cluster_weight_factor = get_cc_cluster_weight_factor(cluster_id) # accounting for health seeking behavior data
                cc_cluster_weight_factor = (cluster_2_mean_pop(cluster_id)/(calib_node_pop + 0.0)) * cc_cluster_weight_factor # accounting for real cluster population (mean across all rounds)
                ccs_model_agg, ccs_ref_agg = get_cc_model_ref_traces(sample_cc_trace, cluster_id, cc_cluster_weight_factor)
                ccs_model_agg_unweighted, ccs_ref_agg_unweighted = get_cc_model_ref_traces(sample_cc_trace, cluster_id)
                
                if not cluster_id in weighted_ccs_model_agg_by_hfca[hfca_id]:  
                    weighted_ccs_model_agg_by_hfca[hfca_id][cluster_id] = {
                                                                           'weighted':[],
                                                                           'unweighted':[]
                                                                           }
                    
                weighted_ccs_model_agg_by_hfca[hfca_id][cluster_id]['weighted'].append(ccs_model_agg)
                weighted_ccs_model_agg_by_hfca[hfca_id][cluster_id]['unweighted'].append((sample_group_key,sample_sim_key,ccs_model_agg_unweighted))

    with open(weighted_ccs_by_hfca_id_file_path, 'w') as w_ccs_f:
        json.dump(weighted_ccs_model_agg_by_hfca, w_ccs_f, indent = 3)
        
    debug_p('DONE getting clinical cases samples based on prevalence optimal regions')
    
    debug_p('Saved clinical cases samples based on prevalence optimal regions to ' + weighted_ccs_by_hfca_id_file_path)
    
    return weighted_ccs_model_agg_by_hfca



def weighted_ccs_combos(hfca_id, weighted_ccs_model_agg_by_hfca):
    
    debug_p('Generating clinical cases combos')
    
    samples_by_cluster = {}
    
    num_clusters = len(weighted_ccs_model_agg_by_hfca) # number of clusters in the given hfca
    
    '''
    # determine number of sample clinical cases timeseries per cluster to satisfy the minimum number of combos (cartesian product cardinality) of clinical case time series per hfcas
    
    
    satisfied = False
    sample_size_ccs_per_cluster = 1
    
    while not satisfied:
        if math.pow(num_clusters, sample_size_ccs_per_cluster) >= min_num_combos_per_hfca:
            satisfied = True
            break
        sample_size_ccs_per_cluster = sample_size_ccs_per_cluster + 1
        debug_p('sample size per cluster is ' + str(sample_size_ccs_per_cluster))
    '''   
    
    # randomly sample clinical case time series per cluster
    for cluster_id, ccs in weighted_ccs_model_agg_by_hfca.iteritems():
        #sample_ccs = random.sample(ccs['weighted'], sample_size_ccs_per_cluster)        
        sample_ccs = random.sample(ccs['weighted'], sample_size)
        samples_by_cluster[cluster_id] = sample_ccs
        
    ccs_combos_hfca = []
    
    shuffled_ccs = samples_by_cluster.values() 
    np.random.shuffle(shuffled_ccs)
    
    # find a set of cartesian products on the randomly sampled clinical cases time series across clusters 
    for ccs_combo in it.product(*shuffled_ccs): # construct a cartesian product from the random samples of weighted clinical case traces within hfca_id
        ccs_combos_hfca.append(ccs_combo)
        if len(ccs_combos_hfca) >= min_num_combos_per_hfca:
            break
    
    debug_p('DONE generating clinical cases combos')
        
    return ccs_combos_hfca, samples_by_cluster