from utils import warn_p, debug_p

from objective import Objective

'''
class Reference: specifies observed/reference data properties: 
- objectives []: data columns of interest
- meta {}: reference meta data (e.g. name of data source, collection/compilation date, etc..)
'''
class Reference:
    
    def __init__(self, objectives = None, meta = None):
        self.meta = meta if meta is not None else {} 
        self.objectives = objectives if objectives is not None else [] 
    
    
    def get_objectives(self):
        return self.objectives
    
    
    def set_objectives(self, objectives):
        self.objectives = objectives
        
    
    def get_meta(self):
        return meta
    

    def set_meta(self, meta = {}):
        self.meta = meta
        
        
    def add_objective(self, name, d_points):
        obj = Objective(name, points = d_points)
        self.objectives.append(obj)
        
    def get_obj_by_name(self, name):
        
        #debug_p('Obj query name: ' + name)
        
        for obj in self.objectives:
            #debug_p('Obj result name: ' + obj.get_name())
            if name == obj.get_name():
                return obj
            
        # objective has not been found; return None
        #debug_p('Obj has not been found!')
        return None 
    
    def to_dict(self):
        ref_dict = {'meta':self.meta}
        
        for obj in self.objectives:
            ref_dict[obj.get_name()] = {'d_points':obj.get_points()}
        
        return ref_dict