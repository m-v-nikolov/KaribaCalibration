import copy

from utils import is_number, val_scale, warn_p, debug_p
from sim_data_2_models import calib_data_2_models_list
from surv_data_2_ref import surv_data_2_ref as d2f


from kariba_settings import load_prevalence_mse, category_2_clusters as c2c, cluster_2_prevs as c2p, objectives_channel_codes, cluster_2_itn_traj, cluster_2_drug_cov, rdt_max
from kariba_model import KaribaModel
from kariba_utils import get_sim_group_key, get_model_params


from fitting_set import FittingSet
from fit import Fit


class KaribaFit:
    
    def __init__(self, category, calib_data, fit_terms = None):
        self.category = category
        self.calib_data = calib_data
        self.fit_terms = fit_terms


    def fit(self):
        
        models_list_prime = calib_data_2_models_list(self.calib_data)
                
        best_fits = {}
        all_fits = {}
        #all_fits = {'fit':{'min_residual':float('inf')}, }
        
        all_fits['min_residual'] = float('inf')
        all_fits['max_residual'] = 0.0
        
        
        
        all_fits['models'] = {}
        
        debug_p('category ' + self.category)
        
        for idx,cluster_id in enumerate(c2c(self.category)):
        
            models_list = copy.deepcopy(models_list_prime)
            
            print "Processing cluster " + cluster_id + "."
            debug_p('Processing cluster ' + cluster_id + " in " + self.category + ".")
            
            itn_traj = cluster_2_itn_traj(cluster_id)
            drug_cov = cluster_2_drug_cov(cluster_id)
            
            # prune models to the ones matching prior data
            cluster_models = []
            for model in models_list:
                model_meta = model.get_meta()
                if model_meta['group_key'] == get_sim_group_key(itn_traj, drug_cov):
                    #debug_p('model id before kariba conversion ' + str(model.get_model_id()))
                    group_key = model_meta['group_key']
                    sim_key = model_meta['sim_key']

                    model = KaribaModel(model, self.calib_data[group_key][sim_key], cluster_id, all_fits = self.fit_terms)
                    
                    #model = kariba_model
                    #debug_p('model id after kariba conversion ' + str(model.get_model_id()))
                    cluster_models.append(model)
                
            surv_data = {}
            all_ref_objs_found = True
            for channel_code in objectives_channel_codes:
                if channel_code == 'prevalence':
                    prev_data = c2p(cluster_id)
                    if prev_data:
                        surv_data[channel_code] = prev_data
                    else:
                        msg = 'Prevalence objective reference data was not found!\n Skipping cluster ' + cluster_id + ' fit!'
                        print msg
                        all_ref_objs_found = False
                else:
                    msg = "Channel objective" + channel_code + " not implemented yet!\nSetting objective reference data to None."
                    warn_p(msg)
                    surv_data[channel_code] = None
            
            # one of the reference objective channels was not found; skipping cluster fit!
            if not all_ref_objs_found:
                continue
                        
            ref = d2f(surv_data)
            
            # adjust highest possible fit to account for RDT+ model in dtk not reflecting reality at the upper end
            obj_prev = ref.get_obj_by_name('prevalence')
            d_points = obj_prev.get_points()
            obj_prev.set_points([min(point, rdt_max) for point in d_points])
            
            
            fitting_set = FittingSet(cluster_id, cluster_models, ref)
            
            if load_prevalence_mse:
                fit = Fit(fitting_set, type = 'mmse_distance_cached')
            else:
                fit = Fit(fitting_set)
            
            best_fit_model = fit.best_fit_mmse_distance()
            
            min_residual = fit.get_min_residual()
            max_residual = fit.get_max_residual()
            
            if min_residual  < all_fits['min_residual']:
                all_fits['min_residual'] = min_residual 
                
            if max_residual  > all_fits['max_residual']:
                all_fits['max_residual'] = max_residual
            
            if best_fit_model: 
            
                temp_h, const_h, itn_level, drug_coverage_level = get_model_params(best_fit_model)
                best_fit_meta = best_fit_model.get_meta()
                best_fits[cluster_id] = {}
                best_fits[cluster_id]['habs'] = {}
                best_fits[cluster_id]['habs']['const_h'] = const_h 
                best_fits[cluster_id]['habs']['temp_h'] = temp_h
                best_fits[cluster_id]['ITN_cov'] = itn_level
                best_fits[cluster_id]['category'] = self.category
                best_fits[cluster_id]['MSAT_cov'] = drug_coverage_level
                best_fits[cluster_id]['sim_id'] = best_fit_meta['sim_id']
                best_fits[cluster_id]['sim_key'] = best_fit_meta['sim_key'] 
                best_fits[cluster_id]['group_key'] = best_fit_meta['group_key']
                best_fits[cluster_id]['fit_value'] = best_fit_model.get_fit_val()
                best_fits[cluster_id]['sim_avg_reinfection_rate'] = best_fit_model.get_sim_avg_reinfection_rate()
                best_fits[cluster_id]['ref_avg_reinfection_rate'] = best_fit_model.get_ref_avg_reinfection_rate()
                best_fits[cluster_id]['prevalence'] = best_fit_model.get_objective_by_name('prevalence').get_points()
            
                # redundancy; to be refactored via FitEntry class                
                best_fits[cluster_id]['fit'] = {}
                best_fits[cluster_id]['fit']['value'] = best_fit_model.get_fit_val()
                best_fits[cluster_id]['fit']['temp_h'] = temp_h
                best_fits[cluster_id]['fit']['const_h'] = const_h
                best_fits[cluster_id]['fit']['ITN_cov'] = itn_level
                best_fits[cluster_id]['fit']['MSAT_cov'] = drug_coverage_level
                best_fits[cluster_id]['fit']['sim_id'] = best_fit_meta['sim_id']
                best_fits[cluster_id]['fit']['sim_key'] = best_fit_meta['sim_key']
                
                
                best_fits[cluster_id]['mse'] = {}
                best_fits[cluster_id]['mse']['value'] = fit.get_min_mses()['prevalence']['value'] # get mmse for objective prevalence
                best_fit_mse_model = fit.get_min_mses()['prevalence']['model']
                temp_h, const_h, itn_level, drug_coverage_level = get_model_params(best_fit_mse_model)
                model_meta_data = best_fit_mse_model.get_meta()
                best_fits[cluster_id]['mse']['temp_h'] = temp_h
                best_fits[cluster_id]['mse']['const_h'] = const_h
                best_fits[cluster_id]['mse']['ITN_cov'] = itn_level
                best_fits[cluster_id]['mse']['MSAT_cov'] = drug_coverage_level
                best_fits[cluster_id]['mse']['sim_id'] = model_meta_data['sim_id']
                best_fits[cluster_id]['mse']['sim_key'] = model_meta_data['sim_key']
                
                best_fits[cluster_id]['cc_penalty'] = {}
                best_fits[cluster_id]['cc_penalty']['value'] = fit.get_min_penalties()['prevalence']['value'] # get clinical penalty for objective prevalence; at present this is just the clinical cases penalty; if reinfection is considered the code needs to be adjusted
                best_fit_cc_penalty_model = fit.get_min_penalties()['prevalence']['model']
                temp_h, const_h, itn_level, drug_coverage_level = get_model_params(best_fit_cc_penalty_model)
                model_meta_data = best_fit_cc_penalty_model.get_meta()
                best_fits[cluster_id]['cc_penalty']['temp_h'] = temp_h
                best_fits[cluster_id]['cc_penalty']['const_h'] = const_h
                best_fits[cluster_id]['cc_penalty']['ITN_cov'] = itn_level
                best_fits[cluster_id]['cc_penalty']['MSAT_cov'] = drug_coverage_level
                best_fits[cluster_id]['cc_penalty']['sim_id'] = model_meta_data['sim_id']
                best_fits[cluster_id]['cc_penalty']['sim_key'] = model_meta_data['sim_key']
                  
    
                rho = best_fit_model.get_rho()
                p_val = best_fit_model.get_p_val()
                
                if rho and p_val :
                    best_fits[cluster_id]['rho'] = rho
                    best_fits[cluster_id]['p_val'] = p_val
                    
                    debug_p('rho' + str(rho))
                    debug_p('p_val' + str(p_val)) 
                
                
            else:
                msg = "something went wrong and the best fit for " + cluster_id + " could not be found."
                warn_p(msg)
                
            
            all_fits['models'][cluster_id] = cluster_models
            #all_fits['models'][cluster_id] = fit.get_fitting_set_models()
            
            print str(idx+1) + " clusters have been processed."
            debug_p( str(idx+1) + " clusters have been processed in category " + self.category)
            
            '''
            if idx > 0:
                break 
            '''      
        return best_fits, all_fits