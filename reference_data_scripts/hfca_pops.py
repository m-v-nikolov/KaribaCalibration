import pandas as pd
import json
from nodenames import *
import numpy as np

df = pd.read_csv('prevalence_net-use_service-use_by_new_cluster_waves1-6.csv')

by_cluster = df.groupby('cluster.id')
cluster_ids = node_names.values() # cluster_ids

#population = by_cluster['pop'].max().to_dict()

hfcas_pops = {}
found_in_data = 0
for cluster_id, cluster_record in by_cluster:
    #print str(name)
    if str(cluster_id) in cluster_ids:
        pops = cluster_record.to_dict()['pop']
        hfca_id = cluster_id.split('_')[0]
        if len(pops.values()) > 0:
            print cluster_id
            mean_pop = np.mean([pop if pop != -1000 else 0 for pop in pops.values()])
            print 'mean pop ' + str(mean_pop)
            
            if not hfca_id in hfcas_pops:
                hfcas_pops[hfca_id] = mean_pop
            else:
                hfcas_pops[hfca_id] = hfcas_pops[hfca_id] + mean_pop
        else:
            print cluster_ids + " not in prevalence data"
            continue
    else:
        continue

print 'HFCAs w/ records in data: '
print len(hfcas_pops)
print

with open('hfcas_pops.json','w') as hfcas_pops_f:
    json.dump(hfcas_pops, hfcas_pops_f, indent = 3)               
    