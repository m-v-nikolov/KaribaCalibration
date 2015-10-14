import pandas as pd
from pandas import *
import numpy as np
import random as rand
import operator

import json
import os
import math
import copy
import sys 


def parse_linked_individuals(all_records_file_path = 'individuals_ids_rounds_rdts_gw_snz.csv', rounds_files_folder = '.'):

    df_surv = pd.io.parsers.read_csv(all_records_file_path)
    df_surv_dict = df_surv.to_dict()
    
    df_surv_person_dict = df_surv['PersonID'].to_dict()
    df_surv_person_dict_rev = dict(zip(df_surv_person_dict.values(), df_surv_person_dict.keys()))
    
    df_surv_rdtpos_dict = df_surv['rdtpos'].to_dict()
    df_surv_round_dict = df_surv['round'].to_dict()
    df_surv_cluster_dict = df_surv['cluster.id'].to_dict()
    
    processed_ids = []
    
    df_rnds_links = {}
    
    for i in range (1,6):
    
        rnd_file = os.path.join(rounds_file_folder, 'round_'+str(i)+'_'+str(i+1)+'.csv')
        
        df_rnds_links[i] = {}
        df_rnds_links[i]['file'] = pd.io.parsers.read_csv(rnd_file).to_dict()
        
        df_rnds_links[i]['fwd'] = {}
        df_rnds_links[i]['fwd']['wave1']  = dict(zip(df_rnds_links[i]['file']['wave1_id'].values(), df_rnds_links[i]['file']['wave1_id'].keys()))
        df_rnds_links[i]['fwd']['wave2']  = dict(zip(df_rnds_links[i]['file']['wave2_id'].values(), df_rnds_links[i]['file']['wave2_id'].keys()))
        
        df_rnds_links[i]['bwd'] = {}
        df_rnds_links[i]['bwd']['wave2'] = dict(zip(df_rnds_links[i]['file']['wave2_id'].values(), df_rnds_links[i]['file']['wave2_id'].keys()))    
        df_rnds_links[i]['bwd']['wave1'] = dict(zip(df_rnds_links[i]['file']['wave1_id'].values(), df_rnds_links[i]['file']['wave1_id'].keys()))
    
    num_entries = len(df_surv_person_dict)
    num_processed = 0
    for idx, person_id_orig in df_surv_person_dict.iteritems():
        
        if person_id_orig in processed_ids:
            num_processed  = num_processed  + 1
            percent_complete = 100*num_processed /(num_entries+0.0)
            sys.stdout.write('\r')
            sys.stdout.write('%2f %%' % percent_complete)
            sys.stdout.flush()
            continue
        
        rnd_orig = df_surv_round_dict[idx]
        rnd = rnd_orig
        person_id = person_id_orig
        
        while rnd < 6:
            
            df_rnd_link = df_rnds_links[rnd]['fwd']['wave1']
    
            if person_id in df_rnd_link:
                link_idx = df_rnd_link[person_id]
                person_id_next_rnd = df_rnds_links[rnd]['file']['wave2_id'][link_idx]
                
                # there are ids in the linked files that are not in the file w/ all longitudinal data?
                if  person_id_next_rnd in df_surv_person_dict_rev:
                    df_surv_dict['PersonID'][df_surv_person_dict_rev[person_id_next_rnd]] = person_id_orig
                person_id = person_id_next_rnd
                processed_ids.append(person_id_next_rnd)
                
                if person_id_orig not in processed_ids:
                    processed_ids.append(person_id_orig)
    
                rnd = rnd + 1
    
            else:
                break
        
        person_id = person_id_orig
        rnd = rnd_orig
        
        while rnd > 1:
            
            rnd = rnd - 1
            df_rnd_link = df_rnds_links[rnd]['bwd']['wave2']
       
            if person_id in df_rnd_link:
                link_idx = df_rnd_link[person_id]
                person_id_prev_rnd = df_rnds_links[rnd]['file']['wave1_id'][link_idx]
                
                # Q: there are ids in the linked files that are not in the file w/ all longitudinal data?
                # A: Yes
                if  person_id_prev_rnd in df_surv_person_dict_rev:
                    df_surv_dict['PersonID'][df_surv_person_dict_rev[person_id_prev_rnd]] = person_id_orig
                person_id = person_id_prev_rnd
                processed_ids.append(person_id_prev_rnd)
                
                if person_id_orig not in processed_ids:
                    processed_ids.append(person_id_orig)
    
            else:
                break
             
        num_processed  = num_processed  + 1
        percent_complete = 100*num_processed /(num_entries+0.0)
        sys.stdout.write('\r')
        sys.stdout.write('%2f %%' % percent_complete)
        sys.stdout.flush()
    
    
    print "Saving linked individuals to ./linked.csv"
    df_surv = pd.DataFrame.from_dict(df_surv_dict)
    df_surv.to_csv("linked.csv")


def analyze_reinfections(linked_individuals_file = 'linked.csv'):

    df_linked = DataFrame.from_csv(linked_individuals_file)
    
    df_linked_by_cluster = df_linked.groupby(['cluster.id'])
    
    clusters_to_proc = len(df_linked_by_cluster)
    
    reinfections = {} 
    for cluster_id, records in df_linked_by_cluster:
        gr_records_by_person = records.groupby(['PersonID'])
    
        reinfections[cluster_id] = {'total':0}

        sys.stdout.write('Processing cluster '  +  cluster_id + '\n')
        sys.stdout.write('\r')    
        num_people = len(gr_records_by_person)
        count = 0
        for person_id, df_linked_records in gr_records_by_person: 
        
            size = len(df_linked_records)
            if size > 1:
                dict_linked_records = df_linked_records.to_dict()
                rounds = sorted(dict_linked_records['round'].items(), key = operator.itemgetter(1))
    
                idxs_sorted, rounds_sorted = zip(*rounds)                 
                for i in range(0,len(idxs_sorted)-1):
                    
                    if not 'reinf_' + str(i+1) + '_' + str(i+2) in reinfections[cluster_id]:
                        reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)] = {'total': 0.0, 'reinf':0.0}
                    
                    if dict_linked_records['rdtpos'][idxs_sorted[i]] == 1 and dict_linked_records['rdtpos'][idxs_sorted[i+1]] == 1:
                        reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['reinf'] = reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['reinf'] + 1
                    reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['total'] = reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['total'] + 1
            
            count = count + 1
            percent_complete = 100*count/(num_people+0.0)
            sys.stdout.write('%2f %%' % percent_complete)
            sys.stdout.write('\r')
            sys.stdout.flush()
        
        sys.stdout.flush()
        print 'Processed cluster ' + cluster_id
        print 'Remaining clusters ' + str(clusters_to_proc) 
            
        clusters_to_proc = clusters_to_proc - 1
    
    print 'Writing reinfections data to file...'
    
    with open('reinfections.json','w') as r_f:
        json.dump(reinfections, r_f)
     
    print 'DONE'
    
def in_bin(bin, x):
    if x >= bin[0] and x < bin[1]:
        return True
    else:
        return False    
    
def analyze_reinfections_hfca_age(linked_individuals_file = 'linked.csv'):

    df_linked = DataFrame.from_csv(linked_individuals_file)
    
    df_linked_by_hfca = df_linked.groupby(['HealthFacCatchm'])
    
    hfcas_to_proc = len(df_linked_by_hfca)
    
    hist_bin_bnds = [0,5,10,15,20,25,30,40,50,60,70]
    bins = []
    bins_dict = {}
    reinfections_by_age_hfca_csv = "hfca_id"
    for i in range(len(hist_bin_bnds) - 1):
        bins.append((hist_bin_bnds[i], hist_bin_bnds[i+1]))
        bins_dict[str(hist_bin_bnds[i])+'_'+ str(hist_bin_bnds[i+1])] = 0
        reinfections_by_age_hfca_csv = reinfections_by_age_hfca_csv + "," + str(hist_bin_bnds[i])+'_'+ str(hist_bin_bnds[i+1])  
    
    reinfections_by_age_hfca_csv = reinfections_by_age_hfca_csv + "total_tested," + "\n"
    
    reinfections = {}
    
    for hfca_id, records in df_linked_by_hfca:
        gr_records_by_person = records.groupby(['PersonID'])
    
        reinfections[hfca_id] = {'total_tested':0}
        reinfections[hfca_id] = {'age_hist':bins_dict}

        sys.stdout.write('Processing cluster '  +  hfca_id + '\n')
        sys.stdout.write('\r')    
        num_people = len(gr_records_by_person)
        count = 0
        for person_id, df_linked_records in gr_records_by_person: 
        
            size = len(df_linked_records)
            if size > 1:
                dict_linked_records = df_linked_records.to_dict()
                rounds = sorted(dict_linked_records['round'].items(), key = operator.itemgetter(1))
    
                idxs_sorted, rounds_sorted = zip(*rounds)                 
                for i in range(0,len(idxs_sorted)-1):
                    
                    if not 'reinf_' + str(i+1) + '_' + str(i+2) in reinfections[hfca_id]:
                        reinfections[hfca_id]['reinf_' + str(i+1) + '_' + str(i+2)] = {'total_tested': 0.0, 'reinf':0.0, 'age_hist':bins_dict}
                    
                    if dict_linked_records['rdtpos'][idxs_sorted[i]] == 1 and dict_linked_records['rdtpos'][idxs_sorted[i+1]] == 1 and (i+1 != 3 and i+2 != 4):
                        reinfections[hfca_id]['reinf_' + str(i+1) + '_' + str(i+2)]['reinf'] = reinfections[hfca_id]['reinf_' + str(i+1) + '_' + str(i+2)]['reinf'] + 1
                        for bin in bins:
                            if in_bin(bin, dict_linked_records['age'][idxs_sorted[i]]):
                                reinfections[hfca_id]['age_hist'][str(bin[0])+'_'+ str(bin[1])] = reinfections[hfca_id]['age_hist'][str(bin[0])+'_'+ str(bin[1])] + 1
                                reinfections[hfca_id]['reinf_' + str(i+1) + '_' + str(i+2)]['age_hist'][str(bin[0])+'_'+ str(bin[1])] = reinfections[hfca_id]['reinf_' + str(i+1) + '_' + str(i+2)]['age_hist'][str(bin[0])+'_'+ str(bin[1])] + 1 
                                break
                    if (i+1 != 3 and i+2 != 4):            
                        reinfections[hfca_id]['reinf_' + str(i+1) + '_' + str(i+2)]['total_tested'] = reinfections[hfca_id]['reinf_' + str(i+1) + '_' + str(i+2)]['total_tested'] + 1
                        reinfections[hfca_id]['total_tested'] = reinfections[hfca_id]['total_tested'] + 1
                        
    
            count = count + 1
            percent_complete = 100*count/(num_people+0.0)
            sys.stdout.write('%2f %%' % percent_complete)
            sys.stdout.write('\r')
            sys.stdout.flush()
        
        
        reinfections_by_age_hfca_csv = reinfections_by_age_hfca_csv + str(hfca_id)
         
        for bin, count in reinfections[hfca_id]['age_hist'].iteritems():
            reinfections_by_age_hfca_csv = reinfections_by_age_hfca_csv + "," + str(count)
        reinfections_by_age_hfca_csv = reinfections_by_age_hfca_csv + "," + reinfections[hfca_id]['total_tested'] + "\n"
        
         
        
        sys.stdout.flush()
        print 'Processed HFCA ' + hfca_id
        print 'Remaining HFCAs ' + str(hfcas_to_proc) 
            
        hfcas_to_proc = hfcas_to_proc - 1
    
    print 'Writing reinfections data to file...'
    
    with open('reinfections.json','w') as r_f:
        json.dump(reinfections, r_f)
    print 'DONE'
    
    print 'Writing reinfections histogram (by age) to file...'
    with open('reinfections_by_age_hfca.csv','w') as rcsv_f:
        rcsv_f.write(reinfections_by_age_hfca_csv)
    print 'DONE'


def get_reinfection_rate_per_cluster(reinfections_file = 'reinfections.json'):

    with open(reinfections_file,'r') as r_f:
        reinfections = json.load(r_f)
    
    reinfection_rates = {}
    for cluster_id, reinfection in reinfections.iteritems():
        #print cluster_id
        rounds_recorded = 0
        reinf_frac = 0.0
        reinfection_rates[cluster_id] = {}
        for i in range(0,len(reinfection.keys())-1):  
            if ('reinf_' + str(i+1) + '_' + str(i+2) in reinfections[cluster_id]) and (i+1 != 3 and i+2 != 4):
                rounds_recorded = rounds_recorded + 1 
                reinf_frac = reinf_frac + reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['reinf']/(reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['total'] + 0.0)
                reinfection_rates[cluster_id][str(i+1) + '_' + str(i+2)] = reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['reinf']/(reinfections[cluster_id]['reinf_' + str(i+1) + '_' + str(i+2)]['total'] + 0.0)
                
        if rounds_recorded != 0:
            #print reinf_frac/(rounds_recorded + 0.0)
            #print '=========='
            reinfection_rates[cluster_id]['avg_rate'] = reinf_frac/(rounds_recorded + 0.0)
        else:
            #print 'No reinfection data for cluster ' + cluster_id
            reinfection_rates[cluster_id]['avg_rate'] = None
    #print reinfection_rates
    return reinfection_rates