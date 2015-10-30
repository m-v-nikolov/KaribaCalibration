import json
import os
import shutil as sh

from utils import warn_p, debug_p, verbose_p, json_list_update

from kariba_settings import ref_data_dir, cc_weight, reinf_weight, sim_data_dir, kariba_viz_dir, d3js_src_dir, d3js_src_files, gazetteer_params_file_name, gazetteer_params_header_file_name, gazetteer_base_map_file_name, gazetteer_sim_mn_base_map_file_name, tags_report_data_file,\
    cc_penalty_model, load_cc_penalty, load_prevalence_mse, load_reinf_penalty


class VizConfig():

    def __init__(self, best_fits, sweep_name):
        self.best_fits = best_fits
        self.sweep_name = sweep_name
        self.base_map = None
        self.sweep_map = None


    def update_d3js(self):
        
        for src_file in d3js_src_files:
            src_file_path = os.path.join(d3js_src_dir, src_file)
            dst_file_path = os.path.join(kariba_viz_dir, src_file)
            sh.copyfile(src_file_path, dst_file_path)
        

    def generate_gazetteer(self):
        
        tags_data_file_path = os.path.join(sim_data_dir, self.sweep_name, tags_report_data_file)
    
        gazetteer_file_path = os.path.join(kariba_viz_dir, gazetteer_params_file_name)
        
        with open(tags_data_file_path, 'r') as t_f:
            tags_data = json.load(t_f)

        # extract penalty model w/o weights
        x = cc_penalty_model.index('_cc_w')
        penalty_model = cc_penalty_model[0:x]
        
        
        # if penalties and mse have been pre-cached and are only loaded then assume model re-weighting 
        # find the corresponding existing model and adjust entry, adding reweighting
        if load_cc_penalty and load_prevalence_mse:
            with open(gazetteer_file_path, 'r') as g_f:
                gazetteer_entries = json.load(g_f)
                
                # if preloading penalty model was successful, than assume the corresponding penalty model is already in gazetteer
                # (add exception throw in the future if penalty model is not in gazetteer
                for entry in gazetteer_entries:
                    debug_p(entry['model'])
                    debug_p(penalty_model)
                    if entry['model'] == penalty_model:
                        model_weight = cc_penalty_model[x+1:]
                        entry['select'].append({'name':model_weight, 'value':self.sweep_name})
                        break
                
                # save updated gazetteer entries
                with open(gazetteer_file_path, 'w') as g_f:
                    json.dump(gazetteer_entries, g_f, indent = 4)
                    
        # otherwise, create a new model entry 
        else:
            gazetteer_nav_entry = {}
            gazetteer_nav_entry['model'] = penalty_model 
            gazetteer_nav_entry['params'] = self.get_model_params_gazetteer(tags_data, penalty_model)
            gazetteer_nav_entry['select'] = [{'name':'', 'value':'unselect'}]
            gazetteer_nav_entry['select'].append({'name':model_weight, 'value':self.sweep_name}) 
            
            if os.path.exists(gazetteer_file_path):
                json_list_update(gazetteer_nav_entry, gazetteer_file_path)
            else:
                with open(gazetteer_file_path, 'w') as g_f:
                    json.dump([gazetteer_nav_entry], g_f, indent = 4)
        


    def get_model_params_gazetteer(self, tags_data, penalty_model):
        
        entry_str = ""
            
        for param, values in tags_data.iteritems():
            
            gazetteer_nav_entry['sweep_name'] = self.sweep_name
        
            if len(values) < 6:
                entry_str = entry_str + "{"
                count = 0
                for v in values:
                    entry_str = entry_str + str(v)
                    if count < len(values)-1:
                       entry_str = entry_str  + ','
                    else:
                       entry_str = entry_str  + '}'
                    count = count + 1
                    
            else:
                # assuming numerical values
                max_val = max(values) 
                min_val = min(values)
                entry_str = entry_str + '[' + str(min_val) + ', ' + str(max_val) + ']'
        
            entry_str = entry_str + "<br />"
            
        entry_str = entry_str + penalty_model + "<br />"
        
        return entry_str
    


    def generate_gazetteer_header(self):
        
        tags_data_file_path = os.path.join(sim_data_dir, self.sweep_name, tags_report_data_file)

        with open(tags_data_file_path, 'r') as t_f:
            tags_data = json.load(t_f)
        
        gazeteer_header_file_path = os.path.join(kariba_viz_dir, gazetteer_params_header_file_name)
        
        gazeteer_header = ""
        
        for param, values in tags_data.iteritems():
            
            gazeteer_header = gazeteer_header + param + ":<br />"
            
        gazeteer_header = gazeteer_header + "Penalty model:<br />"
        gazeteer_header = gazeteer_header + "Penalty weight:<br />"
            
        with open(gazeteer_header_file_path, 'w') as gh_f:
            json.dump(gazeteer_header, gh_f)              
            


    def generate_gazetteer_map(self):
        
        gazetteer_base_map_file_path = os.path.join(ref_data_dir, gazetteer_base_map_file_name)
        
        with open(gazetteer_base_map_file_path ,'r') as map_f:
            self.base_map = json.load(map_f)
            
        self.sweep_map = []
        for cluster_id, cluster_record in self.best_fits.iteritems():
            cluster_map_record = self.get_cluster_map_record(cluster_id)
            
            if cluster_map_record:
                #cluster_map_record['sim_avg_reinfection_rate'] = cluster_record['sim_avg_reinfection_rate']
                #cluster_map_record['ref_avg_reinfection_rate'] = cluster_record['ref_avg_reinfection_rate']
                cluster_map_record['temp_h'] = cluster_record['habs']['temp_h']
                cluster_map_record['const_h'] = cluster_record['habs']['const_h']
                
                mn_sim_map_file_path = os.path.join(kariba_viz_dir,  gazetteer_sim_mn_base_map_file_name)
                if os.path.exists(mn_sim_map_file_path):
                    with open(mn_sim_map_file_path ,'r') as mn_map_f:
                        mn_map = json.load(mn_map_f)
                        cluster_map_record['RDT_mn_sim'] = mn_map[cluster_id]['RDT_sim']
                
                cluster_map_record['itn_level'] = cluster_record['ITN_cov']
                cluster_map_record['drug_coverage'] = cluster_record['MSAT_cov']
                cluster_map_record['fit_value'] = cluster_record['fit_value']
                cluster_map_record['RDT_sn_sim'] = cluster_record['prevalence']
                
                self.sweep_map.append(cluster_map_record)
                 
            else:
                warn_p('No cluster map record found in base map. Skipping generating sim entries in cluster map record.')
                
            
        gazetteer_sweep_map_file_path = os.path.join(kariba_viz_dir, self.sweep_name + '_' + gazetteer_base_map_file_name)
        
        with open(gazetteer_sweep_map_file_path ,'w') as map_f:
            json.dump(self.sweep_map, map_f, indent = 4)
             
            
    
    def get_cluster_map_record(self, cluster_id):
        for cluster_map_record in self.base_map:
            if cluster_map_record["FacilityName"] == cluster_id:
                return cluster_map_record
        return None

    
    
    
    