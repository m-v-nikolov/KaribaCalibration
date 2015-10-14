import pandas as pd
import json
from nodenames import *

pilot_prev_file = 'reference_data_scripts/pilot_clusters.csv'
df = pd.io.parsers.read_csv(pilot_prev_file)

pilot_prevs = {}

for cluster_id in node_names.values():

    cluster_df = df[df['cluster_id'] == cluster_id]
    
    #print cluster_df['rdtpos'].to_dict()
    
    if cluster_df['rdtpos'].to_dict():
        obs_rdtpos = cluster_df['rdtpos'].to_dict().values()[0]
        pilot_prevs[cluster_id] = obs_rdtpos
     
with open('clusters_pilot_prevalences.json','w') as clusters_pilot_prevs_f:
    json.dump(pilot_prevs, clusters_pilot_prevs_f)