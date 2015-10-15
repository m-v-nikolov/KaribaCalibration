import pandas as pd
import json
from nodenames import *
    
prevalence_info_file = 'prevalence_net-use_service-use_by_new_cluster_waves1-6.csv'
df=pd.read_csv(prevalence_info_file)
by_cluster = df.groupby('cluster.id')
clusters = node_names.values() # cluster_ids

#population = by_cluster['pop'].max().to_dict()

cluster_prevalences = {}
found_in_data = 0
for name, cluster_record in by_cluster:
    #print str(name)
    if str(name) in clusters:
        found_in_data = found_in_data + 1 
        cluster = {}
        record = cluster_record[['cluster.id','longitude','latitude', 'rdtpos', 'pop', 'round']]
        pops = record.to_dict()['pop']
        rdts = record.to_dict()['rdtpos']
        rounds = record.to_dict()['round']
        cluster['RDT'] = {}
        cluster['pop'] = {}
        for idx in sorted(rdts):
            cluster['RDT'][str(rounds[idx])] = rdts[idx]
            cluster['pop'][str(rounds[idx])] = int(pops[idx])
    
        if not cluster['RDT'] == {}: 
            cluster_prevalences[name] = cluster
        else:
            print cluster + " not in prevalence data"
            continue
    else:
        continue

print 'Clusters w/ records in data: '
print len(cluster_prevalences)
print

with open('clusters_prevalences.json','w') as cluster_prevs_f:
    json.dump(cluster_prevalences, cluster_prevs_f)