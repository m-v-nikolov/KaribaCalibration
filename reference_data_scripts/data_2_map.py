import pandas as pd
import json
from nodenames import *
    
prevalence_info_file = 'prevalence_net-use_service-use_by_new_cluster_waves1-6.csv'
df=pd.read_csv(prevalence_info_file)
by_cluster = df.groupby('cluster.id')
clusters = node_names.values() # cluster_ids

#population = by_cluster['pop'].max().to_dict()

clusters_map = []
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
        
        cluster['RDT_obs'] = len(rdts)*[0]
        cluster['Population'] = len(pops)*[0]
        cluster['Longitude'] = record.to_dict()['longitude'].values()[0]
        cluster['Latitude'] = record.to_dict()['latitude'].values()[0]
        cluster['FacilityName'] = name
        
        for idx in sorted(rdts):
            cluster['RDT_obs'][int(rounds[idx])-1] = rdts[idx]
            cluster['Population'][int(rounds[idx])-1] = int(pops[idx])
    
        if not cluster['RDT_obs'] == {}: 
            clusters_map.append(cluster)
        else:
            print cluster + " not in prevalence data"
            continue
    else:
        continue

print 'Clusters w/ records in data: '
print len(clusters_map)
print

with open('map.json','w') as clusters_map_f:
    json.dump(clusters_map, clusters_map_f, indent = 4)