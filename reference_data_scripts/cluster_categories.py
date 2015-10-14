import pandas as pd
from pandas import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import json
import os
import math
import operator

from nodenames import *


# group households per cluster
df = pd.io.parsers.read_csv('hh_altitudes.csv')
grouped_hh = df.groupby('cluster.id')

# generate altitude bins of width 100 and find the number of households per cluster 
# falling in each bins

ref_alts = [600, 800, 1000, 1300]
ref_nodes = ['Sinamalima', 'Munumbwe','Lukonde', 'Gwembe']

no_pilot_hfcas = [
                 '80201',
                 '80205',
                 '80208',
                 '80210',
                 '80212',
                 '81103',
                 '81107',
                 '81110',
                 '81112',
                 '81113',
                 '81116'
                 ]

cluster_categories = {
                      'Sinamalima_pilot':[],
                      'Sinamalima_no_pilot':[],
                      'Lukonde_pilot':[],
                      'Lukonde_no_pilot':[],
                      'Munumbwe_pilot':[],
                      'Munumbwe_no_pilot':[],
                      'Gwembe_pilot':[]
                      }

count_clusters = 0
for name, group in grouped_hh:
    #if not name in node_names.values():
    #    continue
    count_clusters = count_clusters + 1
    # get the altitudes of households in the cluster
    group_alt = DataFrame({'altitudes':group.alt.to_dict().values()})

    # label the 12 bins of altitudes from 1 to 12
    labels = range(1,11)
    #bin the households w.r.t. altitude
    group_alt['alt_freq'] = pd.cut(group_alt.altitudes, range(300, 1400, 100), right=False, labels = labels)

    # find the bin with largest number of households
    group_dict = group_alt.groupby('alt_freq').size().to_dict()
    
    alts = range(300, 1400, 100)

    num_hhs = 0
    for alt_label,hhs in group_dict.iteritems():
        if not str(hhs) == 'nan':
            num_hhs = num_hhs + hhs 
    
    avg_alt = 0.0    
    for alt_label,hhs in group_dict.iteritems():
        if not str(hhs) == 'nan':
            weight = hhs/(num_hhs + 0.0)
            alt = alts[alt_label]
            avg_alt = avg_alt + weight * alt 
    
    # if a bin has no households pandas places a value of nan; this is not well handled by python's max function
    # hence, we clean the data a bit to replace nan's
    '''
    inverse = [(float("-inf"), key) if str(value) == 'nan' else (value, key) for key, value in group_dict.items()]
    avg_alt = alts[max(inverse)[1]]
    ''' 
    idx = np.argmin(np.absolute(np.subtract(ref_alts,avg_alt)))
    node = ref_nodes[idx]
    hfca_id = name.split('_')[0]
    pilot = 'no_pilot' if hfca_id in no_pilot_hfcas else 'pilot'
    
    print np.absolute(np.subtract(ref_alts,avg_alt))
    print name
    print avg_alt
    print node
    print idx
    print hfca_id
    print pilot
    print '================================='
    
    cluster_category = node + '_' + pilot
    cluster_categories[cluster_category].append(name)
        
    
    # if a bin has no households pandas places a value of nan; this is not well handled by python's max function
    # hence, we clean the data a bit to replace nan's
    #inverse = [(float("-inf"), key) if str(value) == 'nan' else (value, key) for key, value in group_dict.items()]
    #print inverse
    #cluster_tags[name] = str(max(inverse)[1])

print count_clusters
print cluster_categories

'''
with open( os.path.join( "cluster_categories.json"), "w" ) as cc_f:
        json.dump(cluster_categories,cc_f, indent = 2)
'''