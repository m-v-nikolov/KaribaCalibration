import ast
import os
import json
from datetime import datetime,timedelta

from kariba_settings import cluster_2_cc, cc_sim_start_date, cc_ref_start_date, cc_ref_end_date, tags_report_data_file

from utils import warn_p, debug_p

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
    #debug_p('model cc input cc_aggregate ' + str(len(model_clinical_cases)))
    # aggregate on a periodic basis to as many values as there are in the ref data
    ccs_ref_agg = cluster_2_cc(cluster_id)                               
    #debug_p('ref cc input cc_aggregate ' + str(len(ccs_ref_agg)))         
    ccs_model_agg = []
    ccs_ref_agg_cleaned = []
    
    sim_start_date = cc_sim_start_date
    ref_start_date = cc_ref_start_date
    ref_end_date = cc_ref_end_date
    
    ccs_model = []
    cur_date = sim_start_date         
           
    for i,value in enumerate(model_clinical_cases):
        if i > 0:
            cur_date = cur_date+timedelta(days = 1)
            
        if cur_date >= ref_start_date and cur_date <= ref_end_date:
            ccs_model.append(value)
            
            
    cases_period = 0.0
    periods = 0
    for i, value in enumerate(ccs_model):
        cases_period = cases_period + value
        if (i+1) % (6*7) == 0 or (i+1) == len(ccs_model):
            if periods < len(ccs_ref_agg): 
                if not ccs_ref_agg[periods] == 'nan': 
                    ccs_model_agg.append(cases_period)
                    ccs_ref_agg_cleaned.append(ccs_ref_agg[periods])
                    cases_period = 0.0
                    periods = periods + 1
            else:
                break
    
    
    #debug_p('model cc output cc_aggregate ' + str(len(ccs_model_agg)))
    #debug_p('ref cc output cc_aggregate ' + str(len(ccs_ref_agg_cleaned)))        
    
    return ccs_model_agg, ccs_ref_agg_cleaned 