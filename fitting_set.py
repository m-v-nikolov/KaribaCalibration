from utils import warn_p, debug_p

'''
class FittingSet provides access to models and reference data for fitting

'''

class FittingSet:

    def __init__(self, models_list_name, models_list = None, ref = None):
        self.models_list_name = models_list_name
        self.models_list = models_list if models_list is not None else []
        self.ref = ref
        

    def get_models_list_name(self):
        return self.models_list_name
    
    
    def set_models_list_name(self, models_list_name):
        self.models_list_name = models_list_name
        
    

    '''
    # return a list of models in the form:
    
    [
        # model 1
        {
          'meta': {} # dictionary of model's meta info (sim_id, tags, etc.),
        
          'obj_1': {
                      'm_points':[], # array of model values order-matching reference values 
                      'obj1_weight':float, # a 0-to-1 weight of this objective; reasonably accumulated weights across all objectives should add up to 1 but that is not enforced
                      'm_points_weights':[], # array of model values weights order-matching m_points array 
                      'obj1_fit_penalty',:float # penalty of misfit to obj_1
                  },
          
          
          'obj_2': {
                      'm_points':[], 
                      'obj2_weight':float,
                      'm_points_weights':[],
                      'obj2_fit_penalty':float,
                  },
                  
            ...
            
          'obj_n': {
                      'm_points':[], 
                      'obj2_weight':float,
                      'm_points_weights':[],
                      'obj2_fit_penalty':float,
                  }
        },
        
        # model 2
        {
          'obj_1': {
                      'm_points':[],  
                      'obj1_weight':float, 
                      'm_points_weights':[],  
                      'obj1_fit_penalty':float,
                  },
          
          
          'obj_2': {
                      'm_points':[], 
                      'obj2_weight':float,
                      'm_points_weights':[],
                      'obj2_fit_penalty':float,
                  },
                  
            ...
            
          'obj_n': {
                      'm_points':[], 
                      'obj2_weight':float,
                      'm_points_weights':[],
                      'obj2_fit_penalty':float,
                  }
        },
        
        ...
    ]
    '''
    
    def get_models_list_dict(self):
       
       models_set = []
       for model in self.models:
           models_set.append(model.to_dict())
           
       return models_set

    
    def get_models_list(self):
        return self.models_list
    
    
    def get_ref(self):
        return self.ref
    
    
    
    '''
    # return reference data points for a set of objectives in the form:
    
    {
          'meta': {} # dictionary of ref's meta info (source, collection time stamps etc..),
        
          'obj_1': {
                      'd_points':[] # array of data values 
                   },
          
          
          'obj_2': {
                      'd_points':[] 
                   },
                  
            ...
            
          'obj_n': {
                      'd_points':[] 
                   }
    }
    '''
    
    def get_ref_to_dict(self):
        return self.ref.to_dict()