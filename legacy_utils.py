import sys
import os
import json
import gc
import ast

from kariba_settings import debug_flag, warnings_flag, verbosity_flag
from kariba_utils import get_sim_key, get_sim_group_key

if __name__ == '__main__':
    
    sweep =  {              
          
            "Gwembe_pilot":{
                          'Gwembe_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Gwembe_1_node_MSAT_0.7_w_pilot_calib_a84a10ba-ef48-e511-93f8-f0921c16849c.json'],
                          'Gwembe_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Gwembe_1_node_MSAT_0.55_w_pilot_calib_bc4d78c1-ee48-e511-93f8-f0921c16849c.json'],
                          'Gwembe_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Gwembe_1_node_MSAT_0.35_w_pilot_calib_330dee5a-ee48-e511-93f8-f0921c16849c.json']
                          },

            "Munumbwe_no_pilot":{
                               'Munumbwe_1_node_MSAT_0.7_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.7_w_NO_pilot_calib_5706b2fd-ed48-e511-93f8-f0921c16849c.json'],
                               'Munumbwe_1_node_MSAT_0.55_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.55_w_NO_pilot_calib_d16f43a0-ed48-e511-93f8-f0921c16849c.json'],
                               'Munumbwe_1_node_MSAT_0.35_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.35_w_NO_pilot_calib_55816a27-ed48-e511-93f8-f0921c16849c.json']
                               },
            
            "Munumbwe_pilot":{
                            'Munumbwe_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.7_w_pilot_calib_0e55abc7-ec48-e511-93f8-f0921c16849c.json'],
                            'Munumbwe_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.55_w_pilot_calib_b6195a52-ec48-e511-93f8-f0921c16849c.json'],
                            'Munumbwe_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Munumbwe_1_node_MSAT_0.35_w_pilot_calib_af4a9ce1-eb48-e511-93f8-f0921c16849c.json']
                            },
           
            "Lukonde_pilot":{
                           'Lukonde_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Lukonde_1_node_MSAT_0.7_w_pilot_calib_b324f46a-eb48-e511-93f8-f0921c16849c.json'],
                           'Lukonde_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Lukonde_1_node_MSAT_0.55_w_pilot_calib_44b497de-ea48-e511-93f8-f0921c16849c.json'],
                           'Lukonde_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Lukonde_1_node_MSAT_0.35_w_pilot_calib_7072677b-ea48-e511-93f8-f0921c16849c.json']                           
                           },
          
           "Sinamalima_no_pilot":{
                                 'Sinamalima_1_node_MSAT_0.7_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.7_w_NO_pilot_calib_d1297ec8-e848-e511-93f8-f0921c16849c.json'],
                                 'Sinamalima_1_node_MSAT_0.55_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.55_w_NO_pilot_calib_71ff6f4c-e848-e511-93f8-f0921c16849c.json'],
                                 'Sinamalima_1_node_MSAT_0.35_w_NO_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.35_w_NO_pilot_calib_98ee99d3-e748-e511-93f8-f0921c16849c.json']
                                 },
      
           'Sinamalima_pilot':{
                              'Sinamalima_1_node_MSAT_0.55_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.55_w_pilot_calib_90937da8-e648-e511-93f8-f0921c16849c.json'],
                              'Sinamalima_1_node_MSAT_0.7_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_node_MSAT_0.7_w_pilot_calib_8beb772d-e748-e511-93f8-f0921c16849c.json'],
                              'Sinamalima_1_node_MSAT_0.35_w_pilot':['C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node\\simulations\\Sinamalima_1_nodeMSAT_0.35_w_pilot_calib_a433fc45-e648-e511-93f8-f0921c16849c.json']
                              }
          }
    
    
    legacy_base_dir = "C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\1node"
    updated_base_dir = "C:\\Users\\Mnikolov\\Zambia-raw\\dtk-scripts\\calibration\\sim_data"
    
    for category in sweep:
        
        legacy_calib_file_path = os.path.join(legacy_base_dir, category, 'calibration.json')
        updated_calib_file_path = os.path.join(updated_base_dir, category, 'calibration.json')
        
        with open(legacy_calib_file_path, 'r') as calib_f:
            calib_data = json.load(calib_f)
            
        cc = {}
        for sim_id, sim in calib_data.iteritems():
        
            sim_output = sim['output']
            
            x_temp_h = float(sim['meta']['x_Temporary_Larval_Habitat'])
            
            const_h_struct = ast.literal_eval(sim['meta']['scale_larval_habitats_single'])
            const_h = const_h_struct[0][1][1]
            
            itn_level_struct = ast.literal_eval(sim['meta']['add_ITN_mult'])
            itn_level = itn_level_struct[0][1][0][0][1]
                   
            drug_coverage_level_struct = ast.literal_eval(sim['meta']['add_drug_multi_campaigns'])
            drug_coverage_level = drug_coverage_level_struct[0][1][0]['coverage']
            
            sim_key = get_sim_key(x_temp_h, const_h, itn_level, drug_coverage_level)
            sim_group_key =  get_sim_group_key(itn_level, drug_coverage_level)
         
        
            sim_data = sim_output['Channels']['New Clinical Cases']['Data']
            
            
            if not sim_group_key in cc:
                cc[sim_group_key] = {}

            cc[sim_group_key][sim_key] = sim_data
            
        
        del calib_data
        gc.collect()


        with open(updated_calib_file_path, 'r') as calib_f:
            calib_data = json.load(calib_f)
            
        count_sims = 0
        for group_key, sims in cc.iteritems():
            for sim_key, cc_data in sims.iteritems():
                calib_data[group_key][sim_key]['cc'] = cc_data
                count_sims = count_sims + 1
                
        

        print category
        print count_sims
        
        print 'saving to updated calibration_cc.json'
        
        augmented_calibration = os.path.join(updated_base_dir, category, 'calibration_cc.json')
        
        with open(augmented_calibration, 'w') as calib_f:
            json.dump(calib_data, calib_f)
            
        print 'done'
            
        del calib_data
        gc.collect()            