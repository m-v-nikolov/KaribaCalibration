import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime,timedelta
import seaborn as sns
import json

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False  


#from statsmodels.sandbox.tsa.movstat import va

#df = pd.read_csv('../simulations data/cc_hfca/cc_by_clusters_and_rounds_high_migration_ITN_data_MDA_high_5rnd_clusters.csv')

df = pd.read_csv('cc_data.csv', parse_dates=[3])
#print(df.describe())

clinic_names = sorted(pd.unique(df.shortname.ravel()))
#print(clinic_names)

# Gwembe + Sinazongwe (ordered roughly by descending altitude)
clinics=[ #'Gwembe', 'Munyumbwe', 'Chabbobboma'
        'Gwembe HAHC',
        'Lukonde RHC',
        'Bbondo RHC',
        'Nyanga_Chaamwe RHC',
        'Munyumbwe RHC',
        'Chasanga HP',
        'Sinafala RHC',
        'Chabbobboma RHC',
        'Luumbo RHC',
        'Sinamalima RHC',
        'Chiyabi RHC',
        'Chipepo RHC',
         ]

clinic_2_hfca_id = { 
                      'Gwembe HAHC':80206,
                      'Lukonde RHC':80207,
                      'Bbondo RHC':80201,
                      'Nyanga_Chaamwe RHC':80210,
                      'Chasanga HP':80204,
                      'Munyumbwe RHC':80209,
                      'Sinafala RHC':80211,
                      'Sinamalima RHC':81111,
                      'Chiyabi RHC':81102,
                      'Chipepo RHC':80203,
                      'Chabbobboma RHC':80202,
                      'Luumbo RHC':80208
                     }

'''
sns.set_style('white')
plt.figure('RDT positive cases',figsize=(15,11),facecolor='w')
'''

ccs_agg = {}
ccs_orig = {}
for i,clinic in enumerate(clinics):
    
    #ax=plt.subplot(len(clinics),1,i+1)
    #ax.set_frame_on(False)
    if clinic not in clinic_names:
        print('No clinic called %s' % clinic)
        continue
    # for rapid reporting data
    by_clinic=df[(df.shortname == clinic) & (df.dataelement == 'RDT positive cases')].sort('startdate')
    #print(by_clinic.head())
    tt=[pd.to_datetime(date) for date in by_clinic.startdate]
    
    yy=by_clinic.value.tolist()
    ccs_agg[str(clinic_2_hfca_id[clinic])] = (len(yy)/6 + 1)*[0]
    ccs_orig[str(clinic_2_hfca_id[clinic])] = yy
    
    cases_per_period = 0.0
    periods = 0
    nans_present = 0
    
    if clinic_2_hfca_id[clinic] == 81111:
        yy[0] = 0 # removing a NaN value
        print tt[0]
        print tt[-1]
    
    for i,value in enumerate(yy):
        
        # WE NEED A MORE THROUGH WAY TO DETECT NaNs
        
        if not is_number(value):
            value = 'nan'
            nans_present = nans_present + 1
        else:
            cases_per_period = cases_per_period + value
            
        if (i+1) % 6 == 0 or (i+1) == len(yy):
            if cases_per_period == 0.0 and nans_present > 0:
               cases_per_period = 'nan'
            
            ccs_agg[str(clinic_2_hfca_id[clinic])][periods] = cases_per_period
            cases_per_period = 0.0
            periods = periods + 1
            nans_present = 0
    
    #plt.bar(tt,yy,width=6,color='gray',edgecolor='gray', linewidth=0)
    
    #plt.ylim(0,350)
    
    #ax.get_yaxis().set_visible(False)
    
    '''
    ax.text(0.5,0.6,clinic.split()[0],
        horizontalalignment='center',
        fontweight='bold',
        color='#444444',
        transform=ax.transAxes)

    ax.get_xaxis().tick_bottom()
    if i != len(clinics)-1:
        ax.get_xaxis().set_ticklabels([])
    else:
        pass

    plt.xlim(datetime(2011,5,1),datetime(2015,8,25))
    '''
    
with open('cc.json','w') as cc_f:
    json.dump(ccs_agg, cc_f, indent = 2)

with open('cc_orig.json','w') as cc_f:
    json.dump(ccs_orig, cc_f, indent = 2)

'''
plt.tight_layout()
plt.subplots_adjust(hspace=0.15)
plt.savefig('RDT pos clinical cases data 150825.pdf', format='PDF')
plt.show()
'''