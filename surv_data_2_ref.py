from utils import warn_p, debug_p

from reference import Reference
from kariba_settings import objectives_channel_codes

       
def surv_data_2_ref(surv_data):
        
    # instantiate an empty reference
    ref = Reference()
    
    # add meta data to reference
    # at present just a timestamp
    from time import gmtime, strftime
    ref.set_meta({'timestamp':strftime("%Y-%m-%d %H:%M:%S", gmtime())})
    
    for channel_code in objectives_channel_codes:
        
        # add reference objective corresponding to each model objective
        if channel_code == 'prevalence':
            #debug_p('surv_data ' + str(surv_data['prevalence']))
            d_points = prevalence_surv_data_2_d_points(surv_data['prevalence'])
        else:
            msg = "Channel " + channel_code + " not implemeneted yet!\nSetting reference data to None."
            warn_p(msg)
            d_points = None
        
        #debug_p('adding objective ' + channel_code)
        #debug_p('num d_points ' + str(len(d_points)))
        
        ref.add_objective(channel_code, d_points)
        
    return ref
        
        
def prevalence_surv_data_2_d_points(prev_surv_data):
    
    if prev_surv_data:
        prev_ts = len(prev_surv_data)*[0]
        
        for rnd,prev in prev_surv_data.iteritems():
            if not prev == -1000:
                prev_ts[int(rnd)] = prev;
            else:
                # if no data in reference impute to 'nan'
                prev_ts[int(rnd) ] = 'nan'    
        
        return prev_ts
    else:
        msg = 'no reference data for prevalence available!\n Returning none.'
        warn_p(msg)
        return None