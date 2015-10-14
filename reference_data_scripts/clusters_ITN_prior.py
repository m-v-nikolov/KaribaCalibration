import pandas as pd
import json
from nodenames import *
from dtk.utils.parsers.json2dict import json2dict
import operator

    
net_info_file = 'prevalence_net-use_service-use_by_new_cluster_waves1-6.csv'
nets_calib = json2dict('cluster_tags_net_usage_single_node_updated.json')

df = pd.io.parsers.read_csv(net_info_file)
grouped_rounds = df.groupby('cluster.id')

net_usage = {}

for cluster_id, group in grouped_rounds:
     pop_cluster = group['pop'].to_dict()
     net_usage_cluster = group.UsedLastNight.to_dict()
     for k,v in net_usage_cluster.iteritems():
         if str(v) == 'nan':
             vals = np.array(net_usage_cluster.values())
             net_usage_cluster[k] = vals[~np.isnan(vals)].mean()
         
     rounds_idx = group.round.to_dict()
     
     
     
     for k,v in rounds_idx.iteritems():
         net_usage_cluster[str(v)] = net_usage_cluster.pop(k)
         pop_cluster[str(v)] = pop_cluster.pop(k)
    
     rnd_max_pop = max(pop_cluster.iteritems(), key=operator.itemgetter(1))[0]
     
     net_usage_obs = net_usage_cluster[rnd_max_pop]
     
     print cluster_id
     print net_usage_obs
     
     min_diff = 100
     opt_coverage = 'lowest'
     for coverage, campaign in nets_calib.iteritems():
         print coverage
         print campaign['best_match']['eff_dists'][-1][1]
         print campaign['best_match']['eff_dists'][-1][1] - net_usage_obs
         if (abs(campaign['best_match']['eff_dists'][-1][1] - net_usage_obs) <= min_diff):
             min_diff = abs(campaign['best_match']['eff_dists'][-1][1] - net_usage_obs)
             opt_coverage = coverage
     print '============================='
     net_usage[cluster_id] = opt_coverage
     

with open('itn_prior.json', 'w') as itn_f:
    json.dump(net_usage, itn_f, indent = 3)