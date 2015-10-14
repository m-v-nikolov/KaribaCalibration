import sys
import glob
import shutil
import json

from kariba_settings import debug_flag, warnings_flag, verbosity_flag

def is_number(s):
    
    try:
        float(s)
        return True
    except ValueError:
        return False   
        
    
def feature_scale(feature_vals, max_val = None, min_val = None):
    
    num_vals = len(feature_vals)

    if not max_val:
        max_val = max(feature_vals)
    if not min_val:
        min_val = min(feature_vals)

    if num_vals > 1 and (not max_val - min_val == 0):            
        feature_vals_sc = [(val - min_val + 0.0)/(max_val - min_val) for val in feature_vals]
    
    elif num_vals == 1 or max_val - min_val == 0:
        feature_vals_sc = feature_vals
        
    else:
        feature_vals_sc = None
    
    return feature_vals_sc



def json_list_update(data, file_path):
    
    with open(file_path, 'r') as f_p:
        cur_data = json.load(f_p)
        
    cur_data.append(data)
    
    with open(file_path, 'w') as f_p:
        json.dump(cur_data, f_p, indent = 4)
        


def json_dict_update(data, file_path):
    
    with open(file_path, 'r') as f_p:
        cur_data = json.load(f_p)
        
    cur_data.update(data)
    
    with open(file_path, 'w') as f_p:
        json.dump(cur_data, f_p, indent = 4)
        


def val_scale(val, max_val, min_val):
    
    if not max_val - min_val == 0:
        return (val - min_val + 0.0)/(max_val - min_val)
    else:
        return 0
    
def copy_all(src_dir, dst_dir):
    for f in glob.glob(os.path.join(src_dir, '*.*')):
        shutil.copy(f, dst_dir)
    

def debug_p(msg):
    if debug_flag:
        sys.stderr.write('DEBUG: ' + str(msg) + '\n')
        
def verbose_p(msg):
    if verbosity_flag:
       print 'VERBOSE: ' + str(msg)
        
def warn_p(msg):
    if warnings_flag:
       print 'WARNING: ' + str(msg)
    