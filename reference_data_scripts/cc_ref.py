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

    # the following lines ensure that all missing dates are filled in w/ nans
    # and all dates are sorted chronologically    
    tt_yy = zip(tt, yy)
    tt_yy.sort(key=lambda(x): x[0])
    
    tt.sort()
    
    start_date = min(tt)
    end_date = max(tt)
    
    num_days = ((end_date - start_date)/np.timedelta64(1,'D')).astype(int)
    tt_full = [start_date + timedelta(days = d) for d in range(0, num_days+1)]

    tt_yy_full = len(tt_full)*[0]
    
    for i,t in enumerate(tt_full):
        if t in tt:
            if clinic_2_hfca_id[clinic] == 81111 and i == 0: # manually removing a strange NA value at the first date of 81111 reports
                tt_yy_full[i] = (str(t),'nan')
            else:
                tt_yy_full[i] = (str(t), tt_yy[tt.index(t)][1])
        else: 
            tt_yy_full[i] = (str(t),'nan')
        
    
    ccs_agg[str(clinic_2_hfca_id[clinic])] = (len(tt_yy_full)/(6*7) + 1)*[0]
    ccs_orig[str(clinic_2_hfca_id[clinic])] = tt_yy_full
    
    cases_per_period = 0.0
    periods = 0    
    all_nans = True
    num_weekly_reports = 0 # somewhat redundant with all_nans, but will live with it for now
    for i, (date, cases) in enumerate(tt_yy_full):
        
        if not cases == 'nan':
            cases_per_period = cases_per_period + cases
            num_weekly_reports = num_weekly_reports + 1
            all_nans = False
        
        if (i+1) % (6*7) == 0 or (i+1) == len(tt_yy_full):
            if all_nans:
                ccs_agg[str(clinic_2_hfca_id[clinic])][periods] = (str(date),'nan')
            else:
                ccs_agg[str(clinic_2_hfca_id[clinic])][periods] = (str(date),cases_per_period/(num_weekly_reports + 0.0)) # average by the number of present weeks per bin
            
            cases_per_period = 0.0
            num_weekly_reports = 0
            periods = periods + 1
            all_nans = True
    
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