import math
from utils import debug_p, warn_p

class Fit:
    
    def __init__(self, fitting_set, type = 'mmse_distance'):
        self.fitting_set = fitting_set
        self.type = type
        self.min_residual = float('inf')
        self.max_residual = 0.0
        self.min_penalties = None
        self.min_mses = None
        
    
    def get_type(self):
        return self.type
    
    def set_type(self, type):
        self.type = type
    
    def get_min_residual(self):
        return self.min_residual
    
    def get_max_residual(self):
        return self.max_residual
    
    def get_min_penalties(self):
        return self.min_penalties
    
    def get_min_mses(self):
        return self.min_mses
    
    def fit(self):
        
        best_fit = {}
        
        if self.type == 'mmse_distance':
            best_fit = self.best_fit_mmse_distance()
        elif self.type == 'mmse_distance_cached':
            best_fit = self.best_fit_mmse_distance(cached = True)
        else:
            msg = "unrecognized fit function type " + self.type + "!\nReturning None."
            warn_p(msg)
            return None
        
        return best_fit
            
    
    
    def best_fit_mmse_distance(self, cached = False):
        
        models_list = self.fitting_set.get_models_list()
        
        ref = self.fitting_set.get_ref()
        
        
        min_distance = float('inf')
        min_mses = {}
        min_penalties = {}
        
        
        # assume all models have the same objectives and initialize min obj penalty to inf
        # initialize model minimizer to first model in list
        for obj in models_list[0].get_objectives():
            min_penalties[obj.get_name()] = {} 
            min_penalties[obj.get_name()]['value'] = float('inf')
            min_penalties[obj.get_name()]['model'] = models_list[0]
             
            min_mses[obj.get_name()] = {}
            min_mses[obj.get_name()]['value'] = float('inf')
            min_mses[obj.get_name()]['model'] = models_list[0]
            

        best_fit = None
        
        debug_p(len(models_list))
        debug_p('======================================================================================')
        for model in models_list:
            
            distance = None
            mse = 0.0
            for obj in model.get_objectives():
                m_points = obj.get_points()     # model values for this obj
                d_points = ref.get_obj_by_name(obj.get_name()).get_points()  # reference data values for this obj
                
                #debug_p('m_points len: ' + str(len(m_points)))
                #debug_p('d_points len: ' + str(len(d_points)))
                
                points_weights = obj.get_points_weights()
        
                if not cached:
                    mse = self.mse(m_points, d_points, points_weights)
                else:
                    mse = model.get_cached_mse()
                
                if not mse == None:
                    
                    if mse <= min_mses[obj.get_name()]['value']:
                       min_mses[obj.get_name()]['value'] = mse
                       min_mses[obj.get_name()]['model'] = model
                       
                    if obj.get_model_penalty() <= min_penalties[obj.get_name()]:
                       min_penalties[obj.get_name()]['value'] = obj.get_model_penalty()
                       min_penalties[obj.get_name()]['model'] = model  
                    
                    if distance == None:
                        #debug_p('obj weight: ' + str(obj.get_weight()))
                        #debug_p('obj model penalty: ' + str(obj.get_model_penalty()))
                        distance  = obj.get_weight()*(mse + obj.get_model_penalty())
                    else:
                        #debug_p('obj weight: ' + str(obj.get_weight()))
                        #debug_p('distance ' + str(distance))
                        distance  = distance + obj.get_weight()*(mse + obj.get_model_penalty())
            
            model.set_mse(mse)
            
            if distance:
                
                # a bit redundant since we also find min_distance; will need to adjust
                if distance < self.min_residual:
                    self.min_residual = distance
                    
                if distance > self.max_residual:
                    self.max_residual = distance
                    
                if distance <= min_distance:
                   
                   debug_p('current best distance ' + str(min_distance))
                   if best_fit:
                       #debug_p('current best fit model ' + str(best_fit.get_model_id()))
                       debug_p('current best fit model mmse ' + str(best_fit.get_fit_val()))
                   #debug_p('model ' + str(model.get_model_id()))
                   debug_p('improving fit ' + str(distance))
                   debug_p('fit difference ' + str(distance - min_distance))
                   
                   min_distance = distance
                   model.set_fit_val(min_distance)
                   best_fit = model
                else:
                   model.set_fit_val(distance)
        
        
        #debug_p('best_fit ' + str(best_fit))
        
        self.min_mses = min_mses
        self.min_penalties = min_penalties
        
        return best_fit
                    
        
    def mse(self, m_points, d_points, points_weights):
        
        if d_points:
            
            num_obs = len(d_points)
            
            if not len(m_points) == num_obs:
                msg = "number of points in model does not match num of points in reference data!\nReturning None."
                warn_p(msg)
                return None  
            
            if num_obs == 0:
                msg = "no observations provided!\nReturning None."
                warn_p(msg)
                return None
            
            mse = 0.0
            non_nan_obs = 0
            for idx, m_p in enumerate(m_points):
                d_p = d_points[idx]
                if not (d_p == 'nan' or m_p == 'nan'): 
                    mse = mse + math.pow(m_p - d_p, 2)*points_weights[idx]
                    non_nan_obs = non_nan_obs + 1
            if non_nan_obs > 0:
                mse = mse / non_nan_obs
            else:
                msg = "only nan observations provided!\nReturning None"
                warn_p(msg)
                return None
            
            return mse
        
        
        msg = "no reference data provided for fit!\nReturning None"
        warn_p(msg)
        
        return None
        
    def get_fitting_set(self):
        return self.fitting_set
    
    def get_fitting_set_models(self):
        return self.fitting_set.get_models_list()    