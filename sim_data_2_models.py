import json

from utils import warn_p, debug_p

from model import Model

from kariba_settings import channels_sample_points, reports_channels, objectives_channel_codes 

def sim_channels_2_model(sim_data, sim_key, group_key,):
    
    # instantiate an empty model
    model = Model(sim_key, objectives = [], meta = {})
    #debug_p('model id ' + str(model.get_model_id()))
    #debug_p('model num objectives before adding obj' + str(len(model.get_objectives())))
    
    # add meta data to model
    model.set_meta({'sim_id': sim_data['sim_id'], 'sim_meta':sim_data['meta'], 'sim_key':sim_key, 'group_key':group_key})

    
    for channel_code in objectives_channel_codes:
        
        # add model objective
        m_points = []
        for sample_point in channels_sample_points[channel_code]:
            m_points.append(sim_data[channel_code][sample_point])
        
        # assuming equal weights of objectives
        #debug_p('adding obj to simulation ' + str(sim_key))
        model.add_objective(channel_code, m_points)
        
    #debug_p('model num objectives after adding obj' + str(len(model.get_objectives())))
    
    return model
    

def sim_report_channels_model_format(reports_channels, sim_data):
    
    report_channels_data = {}
    for report_channel in reports_channels:
        if report_channel == 'reinfections':
           report_channels_data['reinfections'] = get_sim_report_reinfections(sim_data)
        else:
            msg = "Channel " + report_channel + " not implemeneted yet!\nSetting report channels data to None."
            warn_p(msg)
            report_channels_data[report_channel] = None
    
    return report_channels_data   


# returns a dictionary of round pairs reinfection rates
def get_sim_report_reinfections(sim_data):
    return sim_data['reinfections']


# parse sim calibration data into a model list
def calib_data_2_models_list(calib_data):
        
    models_list = []
    for group_key, sims in calib_data.iteritems():
        for sim_key, sim_data in sims.iteritems():
            #debug_p('creating model for simulation ' + str(sim_key))
            model = sim_channels_2_model(sim_data, sim_key, group_key)
            
            models_list.append(model)
            
    return models_list     