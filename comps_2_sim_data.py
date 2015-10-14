import os
import json
import sys
import string

from utils import warn_p, debug_p

from dtk.utils.simulation.OutputParser import CompsDTKOutputParser #as parser
from dtk.utils.simulation.COMPSJavaInterop import Experiment, QueryCriteria, Client, Configuration, Priority
from dtk.utils.core.DTKSetupParser import DTKSetupParser
from dtk.utils.parsers.json2dict import json2dict

from kariba_settings import campaign_days, reports_channels, channels, objectives_channel_codes, calibration_data_file
from kariba_utils import get_sim_key, get_sim_group_key, sim_meta_2_itn_level, sim_meta_2_drug_cov, sim_meta_2_temp_h, sim_meta_2_const_h 


def comps_login():
    setup = DTKSetupParser()
    Client.Login(setup.get('HPC','server_endpoint'))
    
# merges a set of sweep calibration files in different base_dirs in a single base_dir
def combine_sweep_results(base_dirs, output_dir):
    
    calib_output = {}
    tags_report_comb = {}
    for base_dir in base_dirs:
        print "Processing calibration in " + base_dir
       
        # get this base_dir's calibration results  
        with open(os.path.join(base_dir, calibration_data_file), 'r') as calib_f:
            calib = json.load(calib_f)
        
        # update each sim group with the respective group sims in this base_dir's calibration results
        for sim_group, sim_keys in calib.iteritems():
            calib_output[sim_group].update(sim_keys)
            
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

    
    print "Writing calibration output to " + os.path.join(output_dir, calibration_data_file)    
    with open(os.path.join(output_dir, calibration_data_file), 'w') as comb_calib_f:
        json.dump(calib_output, comb_calib_f)
    print "DONE"
    
        
    print "Writing tags to " + os.path.join(output_dir, tags_report_data_file)
    with open(os.path.join(output_dir, tags_report_data_file), 'w') as tags_report_f:
        json.dump(tags_report_comb, tags_report_f)            
    print "DONE"    


# get the specified channels/ reports output from a specified experiment sweep 
def get_sweep_results(sim_meta_config_files, calib_file_path, tags_data_file_path):
    
    # Login to COMPs
    comps_login()
    
    
    # find total number of simulations across given experiment files
    num_sims = 0
    for sim_meta_config_file in sim_meta_config_files:
        with open(sim_meta_config_file) as metadata_file:
             metadata = json.loads(metadata_file.read())

        num_sims = num_sims +  len(metadata['sims'])
    
    

    # Download simulations locally
    # sample sim meta config file (like "C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalina_Sinazongwe_Calibration_e9979059-33f8-e411-93f9-f0921c16b9e7.json")
    #print 'Downloading simulations from experiment ' + str(sim_meta_config_files) + '...'

    
    # simulations tag data structure: accumulates sims meta information from sims tags
    tag_data = {
                        'ITN trajectory': [],\
                        'Drug coverage per round': [],\
                        'Temporary habitat scale': [],\
                        'Constant habitat scale': []\
                }
    
    
    # iterate through experiments
    calib_output = {}    
    
    # count processed sims to updated progress
    count = 0
    for sim_meta_config_file in sim_meta_config_files:
        
        
        # construct experiment directory structure
        with open(sim_meta_config_file) as metadata_file:
             metadata = json.loads(metadata_file.read())
    
        output_path = metadata['sim_root']
        exp_id = metadata['exp_id']
        exp_name = metadata['exp_name']
        
        sim_dir_map  = CompsDTKOutputParser.createSimDirectoryMap(exp_id)
        
        # get all successfully completed sims in experiment
        for sim_id, sim in metadata['sims'].items():
            
            # get path to the sim timeseries channels data
            timeseries_path = os.path.join(sim_dir_map[sim_id],'output', 'InsetChart.json')

            
            #get sim timeseries channels data; json2dict returns None if timeseries_path points to non-existing file, which is the case if the sim has not successfully finished 
            sim_output = json2dict(timeseries_path)

            # only download successfully completed  simulations
            if sim_output == None:
                continue

            
            # delete all but the specified channels
            for channel in sim_output['Channels'].keys():
                if not channel in channels:
                   del(sim_output['Channels'][channel])
                          
            # process specified reports
            report_channels_data = {}
            if not reports_channels == None:
               report_channels_data = process_reports(reports_channels, sim_dir_map, sim_id)
            
            # record sim meta information including sim tags
            tags_path = os.path.join(sim_dir_map[sim_id], 'tags.json')
            f = open(tags_path, 'r')
            tags = f.read()
            sim_meta =  ast.literal_eval(tags)
            append_tag_data(sim_meta, tag_data)
        
            
            # construct sim group key and sim key
            x_temp_h = sim_meta_2_temp_h(sim_meta)
            const_h = sim_meta_2_const_h(sim_meta)
            itn_level = sim_meta_2_itn_level(sim_meta)
            drug_coverage_level = sim_meta_2_drug_cov(sim_meta)
            
            sim_key = get_sim_key(x_temp_h, const_h, itn_level, drug_coverage_level)
            sim_group_key =  get_sim_group_key(itn_level, drug_coverage_level)
        

            # store sim channels data  
            if sim_group_key not in calib_output:
                calib_output[sim_group_key] = {}

            calib_output[sim_group_key][sim_key] = {
                                                    # can add/remove data entries depending on needs
                                                    'prevalence': sim_output['Channels']['New Diagnostic Prevalence']['Data'],
                                                    'reinfections': report_channels_data['reinfections'],
                                                    'meta':sim_meta,
                                                    'sim_id':sim_id
                                                    }
    '''    
    count = count + 1
    percent_complete = 100*count/(num_sims+0.0)
    sys.stdout.write('\r')
    sys.stdout.write('%2f %%' % percent_complete)
    #sys.stdout.write('%d' % count)
    sys.stdout.flush()
    '''    
        
    print ""
    print "Writing files..."
    
    with open(calib_file_path, 'w') as calib_f:
            json.dump(calib_output, calib_f)
            print str(len(calib_output)) + ' simulation results saved to ' + calib_file_path
            
    with open(tags_data_file_path, 'w') as tags_f:
            json.dump(tag_data, tags_f)
            print 'Meta data tags saved to ' + tags_data_file_path
    
    print ""
    
    return calib_f


# process specified reports 
def process_reports(reports, sim_dir_map, sim_id):
    
    reports_channels_data = {}    
    
    for report in reports:
        
        if report == 'reinfections':
            reports_channels_data['reinfections'] = process_reinfections_report(sim_dir_map, sim_id)
            
    return reports_channels_data
            
            
# extract re-infections data from patient drug survey reports             
def process_reinfections_report(sim_dir_map, sim_id):

    survey_report_output = {}
    reinfections = {}
    
    for i,day in enumerate(campaign_days):
        
        # get reports before and after campaign
        survey_day_prior = day - 5
        survey_day_after = day + 5
        survey_report_output[survey_day_after] = {} 
        survey_report_prior_path = os.path.join(sim_dir_map[sim_id],'output', 'MalariaSurveyJSONAnalyzer_Day_' + str(survey_day_prior) + '_0.json')
        survey_report_after_path = os.path.join(sim_dir_map[sim_id],'output', 'MalariaSurveyJSONAnalyzer_Day_' + str(survey_day_after) + '_0.json')
        
        survey_report_prior = json2dict(survey_report_prior_path)
        survey_report_after = json2dict(survey_report_after_path)
        
        # debug statements
        #print len(survey_report_prior['patient_array'])
        #print len(survey_report_after['patient_array'])
        
        
        # re-index reports data by patient id
        survey_report_prior_reindexed = {}
        # index patient by id which perhaps is a bit more sensible than the current reporter output
        for patient in survey_report_prior['patient_array']:
            survey_report_prior_reindexed[patient['id']] = patient
        
        survey_report_after_reindexed = {}
        for patient in survey_report_after['patient_array']:
            survey_report_after_reindexed[patient['id']] = patient
            
        survey_report_prior = survey_report_prior_reindexed
        survey_report_after = survey_report_after_reindexed
        

        # determine if each patient in report is reinfected for each pair of consecutive surveys 
        total_tested_patients = 0
        num_reinfected = 0
        count_not_found = 0
        
        for patient_id, patient in survey_report_prior.iteritems():
            treatment_prior = "" 
            for drugs in patient['treatment'][0]:
                if drugs != "": 
                    treatment_prior = drugs
                    
            if patient_id in survey_report_after:
                patient_after = survey_report_after[patient_id]
            else:
                continue
        
            treatment_after = ""
            for drugs in patient_after['treatment'][0]:
                if drugs != "": 
                    treatment_after = drugs
            
            
            if ' + '+treatment_prior in treatment_after and treatment_prior != '':
                treatment_after = string.replace(treatment_after, ' + ' + treatment_prior, '', maxreplace = 1)
            else:
                treatment_after = string.replace(treatment_after, treatment_prior, '', maxreplace = 1)
            #print "clean after treatment " + str(treatment_after)

            survey_report_output[survey_day_after][patient_id] = {} 
            survey_report_output[survey_day_after][patient_id]['initial_age'] = patient['initial_age']
            survey_report_output[survey_day_after][patient_id]['treatment'] = treatment_after
            
            if i > 0:
                survey_day_1 = campaign_days[i - 1] + 5
                survey_day_2 = survey_day_after
                if patient_id in survey_report_output[survey_day_1]:
                    treatment_1 = survey_report_output[survey_day_1][patient_id]['treatment']
                    treatment_2 = treatment_after
                    if treatment_1 != ''  and treatment_2 != '':
                        total_tested_patients = total_tested_patients + 1
                        
                        # DO NOT APPLY FOR MDA CAMAPIGNS WITHOUT MODIFICATION!!!!!!
                        if ('Artemether' in treatment_1) and ('Artemether' in treatment_2):
                                  num_reinfected = num_reinfected + 1
                                  
                        '''
                        # for MDA do something along the lines of 
                        if (('DHA' in patient['treatment'] and 'Vehicle' in patient['treatment']) or 'Artemether' in patient['treatment']) and\
                           (('DHA' in survey_report_output[second_survey_day][patient_id]['treatment'] and 'Vehicle' in survey_report_output[second_survey_day][patient_id]['treatment']) or 'Artemether' in survey_report_output[second_survey_day][patient_id]['treatment']):
                              num_reinfected = num_reinfected + 1
                        '''
                    else:
                        #print "patient " + str(patient_id) + " not found in second_survey day (" + str(second_survey_day) + ") report"
                        count_not_found = count_not_found + 1
        
        if i > 0:
            if total_tested_patients != 0:
                reinfections['round_' + str(i) + '_' + str(i+1)] = num_reinfected/(total_tested_patients + 0.0)
            else:
                reinfections['round_' + str(i) + '_' + str(i+1)] = 'nan'

    return reinfections


def append_tag_data(sim_meta, tag_data):
    
    itn_level_struct = ast.literal_eval(sim_meta['add_ITN_mult'])
    itn_level = itn_level_struct[0][1][0][0][1]
                    
    drug_coverage_level_struct = ast.literal_eval(sim_meta['add_drug_multi_campaigns'])
    drug_coverage_level = drug_coverage_level_struct[0][1][0]['coverage']
    
    x_temp_h = float(sim_meta['x_Temporary_Larval_Habitat'])
    
    const_h_struct = ast.literal_eval(sim_meta['scale_larval_habitats_single'])
    const_h = float(const_h_struct[0][1][1])
    
    if itn_level not in tag_data['ITN trajectory']: 
        tag_data['ITN trajectory'].append(itn_level) 
    if drug_coverage_level not in tag_data['Drug coverage per round']: 
        tag_data['Drug coverage per round'].append(drug_coverage_level)            
    if x_temp_h not in tag_data['Temporary habitat scale']: 
        tag_data['Temporary habitat scale'].append(x_temp_h)                
    if const_h not in tag_data['Constant habitat scale']: 
        tag_data['Constant habitat scale'].append(const_h)