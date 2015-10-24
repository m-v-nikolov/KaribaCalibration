from scipy.stats import spearmanr
import math
import json
from datetime import datetime,timedelta, date

from utils import is_number, val_scale, feature_scale, warn_p, debug_p

from sim_data_2_models import sim_channels_2_model, sim_report_channels_model_format

from kariba_settings import load_cc_penalty, load_prevalence_mse, load_reinf_penalty, cc_correction_factor, cc_agg_fold, cc_weight, reports_channels, calib_node_pop, cluster_2_pops, hfca_2_pop, cluster_2_reinfection_rates, cluster_2_cc, cc_penalty_model, cc_sim_start_date, cc_ref_start_date, cc_ref_end_date,\
    reinf_weight
from kariba_utils import get_cc_model_ref_traces, sim_meta_2_itn_level, sim_meta_2_drug_cov, sim_meta_2_temp_h, sim_meta_2_const_h, error_loading_fit_terms


class KaribaModel:
    
    def __init__(self, model, sim_data, cluster_id, reinfection_penalty = 0.0, reinfection_penalty_weight = 0.0,  clinical_cases_penalty = 0.0, clinical_cases_penalty_weight = 0.0, all_fits = None):

        self.cluster_id = cluster_id
        self.reinfection_penalty = reinfection_penalty
        self.reinfection_penalty_weight = reinfection_penalty_weight
        self.ref_reinfection_num_points = 0
        
        self.rho = None
        self.p_val = None
        
        self.clinical_cases_penalty = clinical_cases_penalty
        self.clinical_cases_penalty_weight = clinical_cases_penalty_weight
        self.ref_clinical_cases_num_points = 0
        self.sim_data = sim_data
        
        self.all_fits = all_fits
        
        self.ref_avg_reinfection_rate = 0.0
        self.sim_avg_reinfection_rate = 0.0 
        
        #debug_p('model id during kariba conversion prior model assignment ' + str(model.get_model_id()))
        self.model = model
        
        model_meta = self.model.get_meta()
        sim_key = model_meta['sim_key']
        
        #debug_p('model id during kariba conversion after model assignment ' + str(self.model.get_model_id()))
        
        # get reinfection rates from sim data, compute reinfection penalty and model penalties
        
        if not reinf_weight == 0:
            model_report_channels = sim_report_channels_model_format(reports_channels, self.sim_data)
            if not load_reinf_penalty:
                self.set_reinfection_penalty(model_report_channels['reinfections'], self.cluster_id)
            else:
                if self.all_fits:
                    self.reinfection_penalty = self.all_fits[self.cluster_id][sim_key]['reinf_penalty']
                    self.reinfection_penalty_weight = reinf_weight
                else:
                    error_loading_fit_terms('reinfection penalty')
                    
        
        if not cc_weight == 0:
            
            if not load_cc_penalty:
                if 'ls_norm' in cc_penalty_model: 
                    self.set_clinical_cases_penalty_by_ls(self.sim_data['cc'], self.cluster_id)
                if 'ls_no_norm' in cc_penalty_model: 
                    self.set_clinical_cases_penalty_by_ls_no_norm(self.sim_data['cc'], self.cluster_id)
                if 'corr' in cc_penalty_model:
                    self.set_clinical_cases_penalty_by_corr(self.sim_data['cc'], self.cluster_id)
                    
                    
            else:
                
                if self.all_fits:
                
                    if 'ls_norm' in cc_penalty_model: 
                        self.clinical_cases_penalty = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['ls_norm']
                    elif 'ls_norm_not_folded' in cc_penalty_model:
                        self.clinical_cases_penalty = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['ls_norm_not_folded']
                    
                    if 'ls_no_norm' in cc_penalty_model: 
                        self.clinical_cases_penalty = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['ls_no_norm']
                    if 'corr_folded' in cc_penalty_model:
                        self.clinical_cases_penalty = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['corr_folded']['penalty']
                        self.rho = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['corr_folded']['rho']
                        self.p_val = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['corr_folded']['p_val']
                    if 'corr_not_folded' in cc_penalty_model:
                        self.clinical_cases_penalty = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['corr_not_folded']['penalty']
                        self.rho = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['corr_not_folded']['rho']
                        self.p_val = self.all_fits[self.cluster_id][sim_key]['fit_terms']['cc_penalty']['corr_not_folded']['p_val']
                        
                    self.clinical_cases_penalty_weight = cc_weight
                else:
                    error_loading_fit_terms('clinical cases penalty')
                
        else:
            self.clinical_cases_penalty = 0.0
            self.clinical_cases_penalty_weight = 0.0
            
                    
        self.set_model_penalties()

    
    # only non None if rank correlation method cc_penalty is used
    def get_rho(self):
        return self.rho
    
    # only non None if rank correlation method cc_penalty is used
    def get_p_val(self):
        return self.p_val
    
    def set_reinfection_penalty(self, model_reinfection_rates, cluster_id):
        
        ref_reinfection_rates = cluster_2_reinfection_rates(cluster_id)
        
        if ref_reinfection_rates:
            cluster_pops = cluster_2_pops(cluster_id)
            reinfection_feature = []
            pop_feature = []
            total_pop = 0.0
            # find max and min values of reinfection rates feature
            
            count_reinf = 0
            for i in range(0,5):    
                if ('reinf_' + str(i+1) + '_' + str(i+2) in model_reinfection_rates) and (i+1 != 3 and i+2 != 4):
        
                    cluster_pop = get_cluster_pop_per_rnd_pair(i+1, i+2)
                    total_pop = total_pop + ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['total']
                    
                    if cluster_pop:
                        pop_feature = pop_feature.append(ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['total']/cluster_pop)

                    ref_reinfection_rate = ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['reinf']/(ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['total'] + 0.0)
                    model_reinfection_rate = model_reinfection_rates['round_' + str(i+1) + '_' + str(i+2)]
                    if(is_number(ref_reinfection_rate) and is_number(model_reinfection_rate)):
                        reinfection_feature.append(ref_reinfection_rate)  
                        reinfection_feature.append(model_reinfection_rate)
                        
                        self.sim_avg_reinfection_rate = self.sim_avg_reinfection_rate + model_reinfection_rate 
                        self.ref_avg_reinfection_rate = self.ref_avg_reinfection_rate + ref_reinfection_rate
                        
                        count_reinf = count_reinf + 1
                        
            if count_reinf != 0:
                self.sim_avg_reinfection_rate = self.sim_avg_reinfection_rate / (count_reinf + 0.0)
                self.ref_avg_reinfection_rate = self.ref_avg_reinfection_rate / (count_reinf + 0.0)
                        
            max_reinf_val = None
            min_reinf_val = None
            if reinfection_feature:        
                max_reinf_val = max(reinfection_feature)
                min_reinf_val = min(reinfection_feature)
            else: # no data observed; penalty is set to 0.0
                self.reinfection_penalty = 0.0
                return
                
            max_pop_val = None
            min_pop_val = None
            if pop_feature:
                max_pop_val = max(pop_feature)
                min_pop_val = min(pop_feature)
            else: # no data observed; penalty is set to 0.0
                self.reinfection_penalty = 0.0
                return
            
                    
            # compute square error between reference and model scaled reinfection features to use as a penalty if there are more than threshold number of people linked
            num_linked_threshold = 40
            
            se_reinfection_rates = []
            self.reinfection_penalty = 0.0
            self.ref_reinfection_num_points = 0.0
            
            for i in range(0,5):        
                # do feature scaling
                if ('reinf_' + str(i+1) + '_' + str(i+2) in model_reinfection_rates) and (i+1 != 3 and i+2 != 4) and ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['total'] > num_linked_threshold:
                    cluster_pop = get_cluster_pop_per_rnd_pair(i+1, i+2)
                    if cluster_pop:
                        ref_reinfection_rate = ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['reinf']/(ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['total'] + 0.0)
                        model_reinfection_rate = model_reinfection_rates['round_' + str(i+1) + '_' + str(i+2)]
                        if(is_number(ref_reinfection_rate) and is_number(model_reinfection_rate)):
                            self.ref_reinfection_num_points = self.ref_reinfection_num_points + 1
            
                            # weight square error se for this round pair proportional to the number of linked people for this round pair over the total number of linked people at this cluster for all rounds
                            # also multiple by a weight in [0,1] depending on how close the number of linked people for this round pair is to the known population of the cluster at these rounds;
                            # the closer the number of linked people the closer the weight to 1; the round pair with closest number of linked people is weighted the most 
                            rnd_pair_weight = (val_scale(ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['total']/(cluster_pop + 0.0), max_pop_val, min_pop_val))*ref_reinfection_rates['reinf_' + str(i+1) + '_' + str(i+2)]['total']/total_pop
                            
                            se = pow(val_scale(ref_reinfection_rate, max_reinf_val, min_reinf_val) - val_scale(model_reinfection_rate, max_reinf_val, min_reinf_val),2)
                            self.reinfection_penalty = self.reinfection_penalty + rnd_pair_weight*se
                        
           
            # weight the reinfection penalty at this cluster based on how much data is available; number of potentially available reinfection measurements is 
            # max_ref_reinfection_points in kariba_settings.py
            
            #debug_p('reinfection penalty ' + str(self.reinfection_penalty))
            
            self.reinfection_penalty_weight = self.ref_reinfection_num_points/(max_ref_reinfection_points + 0.0)
            #debug_p('reinfection penalty weight ' + str(self.reinfection_penalty_weight))
            
            #debug_p('weighted reinfection penalty ' + str(self.reinfection_penalty*self.reinfection_penalty_weight))
            
            return 
        
        else: # no reinfection data found so set penalty to 0
            self.reinfection_penalty = 0.0
            return    



    def get_sim_avg_reinfection_rate(self):
        self.sim_avg_reinfection_rate
    
    def get_ref_avg_reinfection_rate(self):
        self.ref_avg_reinfection_rate

    
    def set_clinical_cases_penalty_by_corr(self, model_clinical_cases, cluster_id):

        ccs_model_agg, ccs_ref_agg = get_cc_model_ref_traces(model_clinical_cases, cluster_id)
        

        '''
        cc_debug_agg = {}
        cc_debug_agg['model'] = ccs_model_agg
        cc_debug_agg['ref'] = ccs_ref_agg
         
        with open('cc_debug_agg_'+cluster_id+'.json' ,'w') as ccd_f:
            json.dump(cc_debug_agg, ccd_f, indent=3)
        '''
        
        '''
        cc_debug_agg_clean = {}
        cc_debug_agg_clean['model_clean'] = ccs_model_agg
        cc_debug_agg_clean['ref_clean'] = ccs_ref_agg
         
        
        with open('cc_debug_agg_clean'+cluster_id+'.json' ,'w') as ccd_f:
            json.dump(cc_debug_agg_clean, ccd_f, indent=3)
        '''
    
        rho, p = spearmanr(ccs_ref_agg, ccs_model_agg)
        
        self.clinical_cases_penalty = 1 - rho
        #debug_p('clinical cases penalty ' + str(self.clinical_cases_penalty))
        
        self.clinical_cases_penalty_weight = cc_weight
        #debug_p('weighted clinical cases penalty ' + str(self.clinical_cases_penalty_weight*self.clinical_cases_penalty)) 
        
        self.rho = rho
        self.p_val = p
        
        if rho > 0.75:
            debug_p('clinical cases rho ' + str(rho))
            debug_p('clinical cases p-value ' + str(p))
            debug_p('clinical cases penalty ' + str(self.clinical_cases_penalty))
            debug_p('weighted clinical cases penalty ' + str(self.clinical_cases_penalty_weight*self.clinical_cases_penalty))
        


        
    def set_clinical_cases_penalty_by_ls(self, model_clinical_cases, cluster_id):

        ccs_model_agg, ccs_ref_agg = get_cc_model_ref_traces(model_clinical_cases, cluster_id)                            
    
        max_ccs_sim = max(ccs_model_agg)
        min_ccs_sim = min(ccs_model_agg)        
        ccs_model_agg_sc = feature_scale(ccs_model_agg, max_ccs_sim, min_ccs_sim)
        
        
        max_ccs_ref = max(ccs_ref_agg)
        min_ccs_ref = min(ccs_ref_agg)
        ccs_ref_agg_sc = feature_scale(ccs_ref_agg, max_ccs_sim, min_ccs_sim)
        
        sse = 0.0
        for i, value in enumerate(ccs_ref_agg_sc):
            se = math.pow(value - ccs_model_agg_sc[i], 2)
            sse = sse + se
            
        rmse = math.sqrt(sse/(len(ccs_ref_agg_sc)+0.0))
    
        #debug_p('clinical cases sum of square errors ' + str(rmse))
        
        self.clinical_cases_penalty = rmse
        self.clinical_cases_penalty_weight = cc_weight
        #debug_p('clinical cases penalty ' + str(self.clinical_cases_penalty))
        
        #self.clinical_cases_penalty_weight = 100
        #debug_p('weighted clinical cases penalty ' + str(self.clinical_cases_penalty_weight*self.clinical_cases_penalty))
        
        
    def set_clinical_cases_penalty_by_ls_no_norm(self, model_clinical_cases, cluster_id):
                                       
        ccs_model_agg, ccs_ref_agg = get_cc_model_ref_traces(model_clinical_cases, cluster_id)
        
        sse = 0.0
        for i, value in enumerate(ccs_ref_agg):
            se = math.pow(value - ccs_model_agg[i], 2)
            sse = sse + se
            
        rmse = math.sqrt(sse/(len(ccs_ref_agg)+0.0))
    
        #debug_p('clinical cases sum of square errors ' + str(rmse))
        
        self.clinical_cases_penalty = rmse
        #debug_p('clinical cases penalty ' + str(self.clinical_cases_penalty))
        
        self.clinical_cases_penalty_weight = 100
        #debug_p('weighted clinical cases penalty ' + str(self.clinical_cases_penalty_weight*self.clinical_cases_penalty))
        
        
    def set_model_penalties(self):
        for obj in self.model.get_objectives():
            
            obj_name = obj.get_name()
            
            if (obj_name == 'prevalence'):
                
                #debug_p('clinical cases to reinfection penalty ratios' + str(self.clinical_cases_penalty*self.clinical_cases_penalty_weight/(self.reinfection_penalty*self.reinfection_penalty_weight)))
                
                #prevalence_model_penalty = self.reinfection_penalty*self.reinfection_penalty_weight + self.clinical_cases_penalty*self.clinical_cases_penalty_weight
                
                prevalence_model_penalty = self.clinical_cases_penalty*self.clinical_cases_penalty_weight
                
                #debug_p('prevalence model penalty ' + str(prevalence_model_penalty))
                
                obj.set_model_penalty(prevalence_model_penalty) 
            else:
                obj.set_model_penalty(0.0)
 
    
    def get_cluster_pop_per_rnd_pair(rnd_1, rnd_2):
        if not cluster_pops[rnd_1] == -1000 and not cluster_pops[rnd_2] == -1000: 
              cluster_pop = (cluster_pops[rnd_1] + cluster_pops[rnd_2])/2.0 
        elif not cluster_pops[rnd_1] == -1000:
            cluster_pop = cluster_pops[rnd_1]
        elif not cluster_pops[rnd_2] == -1000:
            cluster_pop = cluster_pops[rnd_2]
        else:
            cluster_pop = None 
        return cluster_pop
    
    
    def fit_entry(self):
        
        model_meta = self.model.get_meta()
        sim_key = model_meta['sim_key']
        
        temp_h = sim_meta_2_temp_h(model_meta['sim_meta'])
        const_h = sim_meta_2_const_h(model_meta['sim_meta'])
        itn_level = sim_meta_2_itn_level(model_meta['sim_meta'])
        drug_cov = sim_meta_2_drug_cov(model_meta['sim_meta'])
        
        if self.all_fits:
            fit_terms = self.all_fits[self.cluster_id][sim_key]['fit_terms']
        else:
            fit_terms = {}
                        
        if not 'cc_penalty' in fit_terms:
            fit_terms['cc_penalty'] = {}

        if 'ls_norm' in cc_penalty_model: 
            fit_terms['cc_penalty']['ls_norm'] = self.clinical_cases_penalty
        elif 'ls_norm_not_folded' in cc_penalty_model:
            fit_terms['cc_penalty']['ls_norm_not_folded'] = self.clinical_cases_penalty
            
        if 'ls_no_norm' in cc_penalty_model: 
            fit_terms['cc_penalty']['ls_no_norm'] = self.clinical_cases_penalty
             
        if 'corr_folded' in cc_penalty_model:
            if not 'corr_folded' in fit_terms['cc_penalty']:
                fit_terms['cc_penalty']['corr_folded'] = {}
            fit_terms['cc_penalty']['corr_folded']['penalty'] = self.clinical_cases_penalty
            fit_terms['cc_penalty']['corr_folded']['rho'] = self.get_rho()
            fit_terms['cc_penalty']['corr_folded']['p_val'] = self.get_p_val()
            
        if 'corr_not_folded' in cc_penalty_model:
            if not 'corr_not_folded' in fit_terms['cc_penalty']:
                fit_terms['cc_penalty']['corr_not_folded'] = {}
            fit_terms['cc_penalty']['corr_not_folded']['penalty'] = self.clinical_cases_penalty
            fit_terms['cc_penalty']['corr_not_folded']['rho'] = self.get_rho()
            fit_terms['cc_penalty']['corr_not_folded']['p_val'] = self.get_p_val()
            
        fit_terms['reinf_penalty'] = self.reinfection_penalty
        
        fit_terms['mse'] = self.get_mse()
            
        fit_entry = {}
        fit_entry[model_meta['sim_key']] = {
                                             'group_key': model_meta['group_key'],
                                             'sim_id':model_meta['sim_id'],
                                             'fit_val': self.get_fit_val(),
                                             'rho_val' : self.get_rho(), # from most recent run
                                             'p_val' : self.get_p_val(), # from most recent run
                                             'x_temp_h': temp_h,
                                             'const_h': const_h,
                                             'fit_terms':fit_terms,
                                             'itn_level': itn_level,
                                             'drug_cov': drug_cov             
                                             }
        
        return fit_entry
    

    # note: we use composition here instead of inheriting from Model; hence all methods of models that would normally be inherited are made available in KaribaModel
    
    def get_objectives(self):
        return self.model.get_objectives()
    
    def set_objectives(self, objectives):
        self.model.set_objectives(objectives)
        
    def get_meta(self):
        return self.model.get_meta()
    
    def set_meta(self, meta):
        self.model.set_meta(meta)
     
    def set_model_id(self, model_id):
        self.model.set_model_id(model_id)
        
    def get_model_id(self):
        return self.model.get_model_id()
        
    def set_fit_val(self, fit_val):
        self.model.set_fit_val(fit_val)
        
    def get_fit_val(self):
        return self.model.get_fit_val()
    
    def set_mse(self, mse):
        self.model.set_mse(mse)
    
    def get_mse(self):
        return self.model.get_mse()
    
    def get_cached_mse(self):
        
        if self.all_fits:
            return self.all_fits[self.cluster_id][sim_key]['fit_terms']['mse']
        else:
            error_loading_fit_terms('mse')
        
         
        
    def add_objective(self, name, m_points, weight = 0.0, m_points_weights = [], fit_penalty = 0.0):
        self.model.add_objective(name, m_points, weight, m_points_weights, fit_penalty)
        
    def to_dict(self):
        return self.model.to_dict()