import ast
import os
import json
import pandas as pd
import numpy as np

from datetime import datetime,timedelta

from kariba_settings import hfca_2_pop, cc_correction_factor, get_fold_bins, cc_agg_fold, cc_agg_period, fold_start_date, fold_end_date, calib_node_pop, cluster_2_cc, cc_sim_start_date, cc_ref_start_date, cc_ref_end_date, tags_report_data_file, fit_terms_file, sim_data_dir, cc_penalty_model

from utils import warn_p, debug_p

def get_cc_penalty(fit_entry):
    
    if 'corr_folded' in cc_penalty_model:
        return fit_entry['fit_terms']['cc_penalty']['corr_folded']['penalty']
    
    if 'corr_not_folded' in cc_penalty_model:
        return fit_entry['fit_terms']['cc_penalty']['corr_not_folded']['penalty']
    
    else:
        debug_p('No clinical cases penalty found. This should not happen!')
        return None
        

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


def cc_data_aggregate(model_clinical_cases, cluster_id):
    
    # aggregate on a periodic basis to as many values as there are in the ref data
    ccs_ref_agg = cluster_2_cc(cluster_id)                                        
    ccs_model_agg = []
    ccs_ref_agg_cleaned = []
    
    hfca_id = cluster_id.split('_')[0]
    
    hfca_pop = hfca_2_pop(hfca_id)
    
    pop_norm_factor = cc_correction_factor*(hfca_pop + 0.0)/calib_node_pop
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
    else:
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_model'] = fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_model'] + cases_model
    
    if fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] =='nan':
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] = cases_ref
    else:
        fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] = fold_bins[num_days_in_fold/(cc_agg_period * 7)]['cases_ref'] + cases_ref
    
    
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
           ccs_folded_model.append(fold_bins[fold_idx]['cases_model'])
           ccs_folded_ref.append(fold_bins[fold_idx]['cases_ref'])
           
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


def get_cc_model_ref_traces(model_clinical_cases, cluster_id):
    
    ccs_model_agg, ccs_ref_agg, fold_bins = cc_data_aggregate(model_clinical_cases, cluster_id)
    
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

    
def error_loading_fit_terms():
    
    debug_p('Could not load fit terms! Check whether fit terms file at ' + os.path.join(sim_data_dir, fit_terms_file) + ' is accessible.')
    
    raise ValueError('Could not load fit terms! Check whether fit terms file at ' + os.path.join(sim_data_dir, fit_terms_file) + ' is accessible.')