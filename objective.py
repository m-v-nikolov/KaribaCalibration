from utils import warn_p, debug_p, verbose_p

'''
class Objective: specifies a single model objective properties: 
- [] m_points: array of model's objective values ordered w.r.t. a reference data set values for the same objective 
- float weight: relative objective weight; reasonably, the set of all model objectives' weights should add up to one
- [] m_points_weights: array of weights corresponding to each of the m _points  
'''

class Objective:

    def __init__(self, name, points = [], points_weights = [], weight = 0.0, fit_penalty = 0.0):
        
        self.name = name
        
        if not points:
            msg = "No points provided. Cannot initialize empty objective."
            verbose_p(msg)
            return
        
        self.points = points
            
        if not points_weights:
            msg = "No points weights provided. Defaulting to equal weights."
            verbose_p(msg)
            self.points_weights = len(self.points)*[(1.0/(len(self.points) + 0.0))]
        else:
            self.points_weights = points_weights
            
        self.weight = weight
        
        self.fit_penalty = fit_penalty
            
    
    def get_name(self):
        return self.name
    
    
    def set_name(self, name):
        self.name = name
        
        
    def get_points(self):
        return self.points
    
    
    def set_points(self, points = []):
        if not points:
            msg = "No points provided. Cannot initialize empty objective."
            verbose_p(msg) 
            return
        self.points = points
        
    
    def get_weight(self):
        return self.weight
    
    
    def set_weight(self, weight):
        self.weight = weight
        
        
    def get_model_penalty(self):
        return self.fit_penalty
    
    
    def set_model_penalty(self, fit_penalty):
        self.fit_penalty = fit_penalty
         
    
    def get_points_weights(self):
        return self.points_weights
    
    
    def set_points_weights(self, points_weights = []):
        if not points_weights:
            msg = "No points weights provided. Defaulting to equal weights."
            verbose_p(msg)
            self.set_equal_points_weights()
        else:
            self.points_weights = points_weights
            
    def set_equal_points_weights(self):
        self.points_weights = len(self.points)*[(1.0/(len(self.points) + 0.0))]