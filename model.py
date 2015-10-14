from utils import warn_p, debug_p

from objective import Objective

'''
class Model: specifies model's properties: 
- objectives
- model meta info: (e.g. simulation ids, tags, etc.) 
'''

class Model:

    def __init__(self, model_id, objectives = None, meta = None):
        
        self.objectives = objectives if objectives is not None else []
        self.meta = meta if meta is not None else {}
        self.model_id = model_id
        #debug_p('creating model with id ' + str(self.model_id))
        #debug_p('num objectives on instantiation ' + str(len(self.objectives)))
        #debug_p('num objectives passed on instantiation ' + str(len(objectives)))
        
        
        self.fit_val = float('inf')
    
    
    def get_objectives(self):
        return self.objectives
    
    def set_objectives(self, objectives):
        self.objectives = objectives
        
    
    def get_meta(self):
        return self.meta
    
    def set_meta(self, meta):
        self.meta = meta
        
    
    def add_objective(self, name, m_points, weight = None, m_points_weights = [], fit_penalty = 0.0):
        obj = Objective(name, points = m_points, weight = weight, points_weights = m_points_weights, fit_penalty = fit_penalty)
        
        #debug_p('adding objective to model ' + str(self.get_model_id()))
        
        self.objectives.append(obj)
        #debug_p('model id (in model) ' + str(self.get_id()))
        #debug_p('num objectives ' + str(len(self.objectives)))
        if weight == None:
            # setting all objectives to equal weight
            eq_weight = 1/(len(self.objectives) + 0.0)
            for obj in self.objectives:
                obj.set_weight(eq_weight)
        
    
    def get_model_id(self):
        return self.model_id
    
    def set_model_id(self, model_id):
        self.model_id = model_id
    
        
    def set_fit_val(self, fit_val):
        self.fit_val = fit_val
    
    def get_fit_val(self):
        return self.fit_val    
        
        
    def to_dict(self):
        model_dict = {'meta': model.get_meta()}
           
        for obj in self.objectives:
            model_dict[obj.get_name()] = {
                                             'm_points':obj.get_m_points(),
                                             'weight':obj.get_weight(),
                                             'm_points_weights':obj.get_m_points_weights(),
                                             'fit_penalty':obj.get_fit_penalty(),
                                             }
        return model_dict 