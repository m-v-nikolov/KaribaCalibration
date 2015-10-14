import pandas as pd
import json
import numpy as np
from nodenames import *
    

calib_coverages = [0.35, 0.55, 0.7]

msat_cov_file = 'linked_and_covered_proportions_byHCCA.csv'

df = pd.io.parsers.read_csv(msat_cov_file)

msat_coverages = {}

for node_id, cluster_id in node_names.iteritems():
    
    hfca = cluster_id.split('_')[0]
    
    hfca_df = df[df['hcca'] == int(hfca)]
    
    obs_avg_cov = hfca_df['avg_coverage'].to_dict().values()[0]
    
    cidx = np.argmin(np.absolute(np.subtract(calib_coverages,obs_avg_cov)))
    
    msat_coverages[cluster_id] = calib_coverages[cidx]
    

with open('drug_coverages.json','w') as cluster_covs_f:
    json.dump(msat_coverages, cluster_covs_f)