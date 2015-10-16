import os
import sys

import numpy as np
import pandas as pd

from scipy import interpolate
from scipy.interpolate import Rbf
from scipy.interpolate import interp2d

from datetime import datetime,timedelta

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.mlab import griddata
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MultipleLocator
import matplotlib.gridspec as gridspec
import brewer2mpl as b2mpl

from utils import warn_p, debug_p

from surv_data_2_ref import surv_data_2_ref as d2f

from kariba_settings import cc_subopt_traces_plots_dir, opt_marker, opt_marker_size, markers, subopt_plots_threshold, cc_penalty_model, hfca_id_2_facility, cluster_2_prevs as c2p, traces_plots_dir, traces_base_file_name, cc_traces_plots_dir, cc_traces_base_file_name, err_surfaces_plots_dir, err_surfaces_base_file_name, sim_data_dir, calibration_data_file, tags_data_file, channels_sample_points, objectives_channel_codes, reports_channels, channels, cc_sim_start_date, cc_ref_start_date, cc_ref_end_date
from kariba_utils import cc_data_aggregate, cc_data_nan_clean

class PlotUtils():
    
    def __init__(self, best_fits, all_fits, calib_data, root_sweep_dir, category):
        
        self.best_fits = best_fits
        self.all_fits = all_fits
        self.calib_data = calib_data
        self.root_sweep_dir = root_sweep_dir
        self.category = category
        
        self.fit_entries_2_markers = {}
     
    def get_ref(self, cluster_id):

        surv_data = {}
        all_ref_objs_found = True
        for channel_code in objectives_channel_codes:
            if channel_code == 'prevalence':
                prev_data = c2p(cluster_id)
                if prev_data:
                    surv_data[channel_code] = prev_data
                else:
                    msg = 'Prevalence objective reference data was not found!\n Skipping plotting cluster ' + cluster_id + ' fit!'
                    print msg
                    all_ref_objs_found = False
            else:
                msg = "Channel objective" + channel_code + " not implemented yet!\nSetting objective reference data to None for plotting."
                warn_p(msg)
                surv_data[channel_code] = None
        
        # one of the reference objective channels was not found; skipping cluster fit!
        if not all_ref_objs_found:
            ref = None      
        else:  
            ref = d2f(surv_data)
            
        return ref       

    
    def plot_calib_prev_traces(self):
        
        count = 0
        for cluster_id, cluster_record in self.best_fits.iteritems():
            
            debug_p('Plotting clinical cases trace for cluster ' + cluster_id + ' in category ' + self.category)
        
            #fig = plt.figure(cluster, figsize=(11, 4), dpi=100, facecolor='white')
            fig = plt.figure(cluster_id, figsize=(9.2, 4), dpi=100, facecolor='white')
            
            gs = gridspec.GridSpec(1, 2)
            ymax = 16
        
            scale_int = np.array(range(0,ymax+1))
            pal = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=ymax+1)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pal)
                         
            ax = plt.subplot(gs[0:2])
            
    
            ax.set_ylim(0,1)
            ax.set_xlim(2180,3010)    
            
            opt_sim_key = cluster_record['sim_key']
            opt_group_key = cluster_record['group_key']
            
            opt_const_h = cluster_record['habs']['const_h']
            opt_x_temp_h = cluster_record['habs']['temp_h']
            opt_itn = cluster_record['ITN_cov']
            opt_drug = cluster_record['MSAT_cov']
            
            opt_prev_trace = self.calib_data[opt_group_key][opt_sim_key]['prevalence']
            
            ax.plot(range(2180, 3000), opt_prev_trace[2179:2999], alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn))
            
            ref = self.get_ref(cluster_id)
            
            if not ref:
                # no reference object could be constructed; no data found; return w/o plotting
                return
            
            obs_prevs = ref.to_dict()['prevalence']['d_points']
            
            label_obs_prev_shown = False
            label_sim_prev_shown = False
            for i,prev in enumerate(obs_prevs):
                
                if prev != 'nan':
                    if not label_obs_prev_shown:
                        label_obs_prev = 'Observed prevalence'
                        label_obs_prev_shown = True
                    else: 
                        label_obs_prev = None
                         
                    ax.scatter(channels_sample_points['prevalence'][i], prev, c = 'red', facecolor = 'red', marker='o', s = 40, label = label_obs_prev)
                    
                if not label_sim_prev_shown:
                    label_sim_prev = 'Simulated prevalence wrt rnds.'
                    label_sim_prev_shown = True
                else: 
                    label_sim_prev = None
                
                ax.scatter(channels_sample_points['prevalence'][i], opt_prev_trace[channels_sample_points['prevalence'][i]], c = 'black', facecolor = 'none', marker='o', s = 60, label = label_sim_prev)
                                        
            plt.xlabel('Time (days)', fontsize=8)
            plt.ylabel('Prevalence (population fraction)', fontsize=8)
            plt.legend(loc=1, fontsize=8)
            plt.title('Prevalence timeseries', fontsize = 8, fontweight = 'bold', color = 'black')
            plt.gca().tick_params(axis='x', labelsize=8)
            plt.gca().tick_params(axis='y', labelsize=8)
        
            plt.tight_layout()
            output_plot_file_path = os.path.join(self.root_sweep_dir, traces_plots_dir, traces_base_file_name + cluster_id + '.png')
            plt.savefig(output_plot_file_path, format='png')
            plt.close()
    
            count = count + 1
    
    
    def plot_calib_err_surfaces(self):
        
        count = 0
        for cluster_id, cluster_records in self.all_fits.iteritems():
            
            debug_p('Plotting error surface for cluster ' + cluster_id + ' in category ' + self.category)
        
            fig = plt.figure(cluster_id, figsize=(4.35, 4), dpi=100, facecolor='white')
            
            gs = gridspec.GridSpec(1, 1)
         
            error_points = {}
            
            title_ITN = 'ITN distribution: '
            title_drug_coverage = 'drug coverage: '
            
            
            for fit_entry in cluster_records:
        
                x_temp_h = fit_entry['x_temp_h']
                const_h = fit_entry['const_h']
                fit_val = fit_entry['fit_val']

                itn_level = fit_entry['itn_level']
                drug_coverage_level = fit_entry['drug_cov']
                group_key = fit_entry['group_key']    
                
                    
                if group_key not in error_points:
                    error_points[group_key] = {
                                                   'x_temp_h':[],
                                                   'const_h':[],
                                                   'fit_val':[],
                                                   'title': title_ITN + itn_level + "; " + title_drug_coverage + str(drug_coverage_level),
                                                   'itn_level':itn_level,
                                                   'drug_coverage':drug_coverage_level
                                                   }
                error_points[group_key]['x_temp_h'].append(x_temp_h)
                error_points[group_key]['const_h'].append(const_h)
                error_points[group_key]['fit_val'].append(fit_val)
                    
                          
            ymax = 10
            
            scale_int = np.array(range(1,10))
            pal = cm = plt.get_cmap('nipy_spectral') 
            cNorm  = colors.Normalize(vmin=0, vmax=ymax)
            #scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pal)
            scalarMap = b2mpl.get_map('Spectral', 'Diverging', 5).mpl_colors
            
            opt_group_key = self.best_fits[cluster_id]['group_key']
            
            for i,group_key in enumerate(error_points.keys()):
                
                itn_level = error_points[group_key]['itn_level']
                drug_coverage = error_points[group_key]['drug_coverage']
                
                if group_key == opt_group_key: 
        
                    plt.subplot(gs[0,0])
                    x = error_points[group_key]['x_temp_h']
                    y = error_points[group_key]['const_h']
                    z = error_points[group_key]['fit_val']
                    #print len(z)
                    res = 125
                    ttl = error_points[group_key]['title']
                    
                    min_x = np.min(x)
                    min_y = np.min(y)
                    
                    max_x = np.max(x)
                    max_y = np.max(y)
                
                    #f = interpolate.interp2d(x, y, z)
                
                    xi = np.linspace(min_x, max_x , res)
                    yi = np.linspace(min_y, max_y , res)
                    
                    zi = griddata(x,y,z,xi,yi)
                    
                    #xig, yig = np.meshgrid(xi, yi)
                    #zig = f(xi,yi)

                    #rbf = Rbf(x, y, z, epsilon=2)
                    #zig = rbf(xig, yig)
                
                    #rmse_pl = plt.pcolor(xi, yi, zi, cmap=plt.get_cmap('RdYlGn_r'), vmin=0.475, vmax=0.8)
                    rmse_pl = plt.pcolor(xi, yi, zi, cmap=plt.get_cmap('RdYlGn_r'))
                    #rmse_pl = plt.contourf(xi,yi,zi,15,cmap=plt.cm.hot)
                    #rmse_pl.cmap.set_over('black')
                    #rmse_pl.cmap.set_under('grey')
                    cb = plt.colorbar(rmse_pl)

                    cb.set_label('Calibration-simulation distance', fontsize=8)
                    cb.ax.tick_params(labelsize=8)    
                    #plt.scatter(x, y, 10, z, cmap=plt.get_cmap('hot'))
                    #plt.title(ttl, fontsize = 8, fontweight = 'bold', color = 'white', backgroundcolor = scalarMap.to_rgba(scale_int[itn_levels_2_sbplts[itn_level]]))
                    plt.title(ttl, fontsize = 8, fontweight = 'bold', color = 'black')
                    plt.xlabel('All habitats scale', fontsize=8)
                    plt.ylabel('Constant habitat scale', fontsize=8)
                    plt.xlim(min_x+0.1, max_x+0.1)
                    plt.ylim(min_y+0.1, max_y+0.1)
                    #plt.ylim(0.01, 14)
                    plt.gca().tick_params(axis='x', labelsize=8)
                    plt.gca().tick_params(axis='y', labelsize=8)
                
                #plt.subplot(gs[itn_levels_2_sbplts[best_fit_itn_level], 0])
                plt.subplot(gs[0,0])
                
                cluster_record = self.best_fits[cluster_id]
                opt_const_h = cluster_record['habs']['const_h']
                opt_x_temp_h = cluster_record['habs']['temp_h']
                opt_itn = cluster_record['ITN_cov']
                opt_drug = cluster_record['MSAT_cov']
                plt.scatter(opt_x_temp_h, opt_const_h, c = 'black', marker = 'd', s = 40, facecolor='none', zorder=100, label='Best fit')
                plt.legend(bbox_to_anchor=(0., 1, 1., .1), loc=3, ncol=2, mode="expand", borderaxespad=0., fontsize=8)
                    
                plt.tight_layout()
                output_plot_file_path = os.path.join(self.root_sweep_dir, err_surfaces_plots_dir, err_surfaces_base_file_name + cluster_id +'.png')
                plt.savefig(output_plot_file_path, format='png')
                plt.close()
        
            count = count + 1
                
                
    def plot_calib_cc_traces_clusters(self):
        
        for cluster_id, cluster_record in self.best_fits.iteritems():
            
            debug_p('Plotting clinical cases trace for cluster ' + cluster_id + ' in category ' + self.category)
            
            fig = plt.figure(cluster_id, figsize=(9.2, 4), dpi=100, facecolor='white')
            
            opt_sim_key = cluster_record['sim_key']
            opt_group_key = cluster_record['group_key']
            
            opt_cc_trace = self.calib_data[opt_group_key][opt_sim_key]['cc']
            
            ccs_model_agg, ccs_ref_agg = cc_data_aggregate(opt_cc_trace, cluster_id)
            ccs_model_agg, ccs_ref_agg  = cc_data_nan_clean(ccs_model_agg, ccs_ref_agg, cluster_id)
            #debug_p('model length ' + str(len(ccs_model_agg)))
            #debug_p('ref length ' + str(len(ccs_ref_agg)))
            
            hfca_id = cluster_id.split('_')[0]
            
            facility = hfca_id_2_facility(hfca_id)
            
            ref_start_date = cc_ref_start_date
            ref_end_date = cc_ref_end_date
            
            '''
            cur_date = ref_start_date         
            dates = [cur_date]       
            
            for i,value in enumerate(ccs_ref_agg):
                if i > 0:
                    cur_date = cur_date+timedelta(days = 6*7)
                    dates.append(cur_date)
            '''
            gs = gridspec.GridSpec(1, 4)
            ymax = 16
        
            scale_int = np.array(range(0,ymax+1))
            pal = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=ymax+1)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pal)
                         
            ax = plt.subplot(gs[0:4])
            
            #ax.set_ylim(1000)
            
            opt_sim_key = cluster_record['sim_key']
            opt_group_key = cluster_record['group_key']
            
            opt_const_h = cluster_record['habs']['const_h']
            opt_x_temp_h = cluster_record['habs']['temp_h']
            opt_itn = cluster_record['ITN_cov']
            opt_drug = cluster_record['MSAT_cov']
            
            # the following code only relevant for rank correlation cc penalty fit
            opt_rho = None
            opt_p_val = None
            if 'rho' in cluster_record:
                opt_rho = cluster_record['rho']
            if 'p_val' in cluster_record:
                opt_p_val = cluster_record['p_val']
            
            
            '''
            mod_dates, mod_cases = zip(*ccs_model_agg)
            ref_dates, ref_cases = zip(*ccs_ref_agg)
            
            debug_p('model dates to print ' + str(mod_dates))
            debug_p('model cases to print ' + str(mod_cases))

            debug_p('ref dates to print ' + str(ref_dates))
            debug_p('ref cases to print ' + str(ref_cases))
            '''
            '''
            if opt_rho and opt_p_val:
                ax.plot(pd.to_datetime(mod_dates), mod_cases, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn) + ', rho=' + str(opt_rho) + ', p-val=' + str(opt_p_val), marker = opt_marker, markersize = opt_marker_size)
            else:
                ax.plot(pd.to_datetime(mod_dates), mod_cases, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn), marker = opt_marker, markersize = opt_marker_size)
            ax.bar(pd.to_datetime(ref_dates), ref_cases, width=12,color='red',edgecolor='red', linewidth=0, label = 'Observed in ' + facility)
            '''
            
            if opt_rho and opt_p_val:
                ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn) + ', rho=' + str(opt_rho) + ', p-val=' + str(opt_p_val), marker = opt_marker, markersize = opt_marker_size)
            else:
                ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn), marker = opt_marker, markersize = opt_marker_size)
            ax.plot(range(0, len(ccs_ref_agg)), ccs_ref_agg, alpha=1, linewidth=2.0, c = 'red', label = 'Observed in ' + facility)
            
            #ax.bar(range(0, len(ccs_ref_agg)), ccs_ref_agg, width=12,color='red',edgecolor='red', linewidth=0, label = 'Observed in ' + facility)
            
            
            '''
            opt_prev_trace = self.calib_data[opt_group_key][opt_sim_key]['prevalence']
            if opt_rho and opt_p_val:
                ax.plot(dates, ccs_model_agg, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn) + ', rho=' + str(opt_rho) + ', p-val=' + str(opt_p_val))
            else:
                ax.plot(dates, ccs_model_agg, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn))
            ax.bar(dates, ccs_ref_agg, width=12,color='red',edgecolor='red', linewidth=0, label = 'Observed in ' + facility)
            #ax.plot(dates, ccs_ref_agg, alpha=1, linewidth=2.0, c = 'red', label = 'Observed in ' + facility)
            '''
            
            plt.xlabel('Time (6-week bins)', fontsize=8)
            plt.ylabel('Clinical cases', fontsize=8)
            plt.legend(loc=1, fontsize=8)
            plt.title('Clinical cases timeseries', fontsize = 8, fontweight = 'bold', color = 'black')
            plt.gca().tick_params(axis='x', labelsize=8)
            plt.gca().tick_params(axis='y', labelsize=8)
            plt.tight_layout()
            output_plot_file_path = os.path.join(self.root_sweep_dir, cc_traces_plots_dir, cc_traces_base_file_name + cluster_id + '.png')
            plt.savefig(output_plot_file_path, format='png')
            plt.close()
            
            
    
    def plot_calib_cc_traces_clusters_opt_neigh(self):
        
        for cluster_id, cluster_record in self.best_fits.iteritems():
            
            debug_p('Plotting clinical cases trace for cluster ' + cluster_id + ' in category ' + self.category)
            
            fig = plt.figure(cluster_id, figsize=(9.2, 4), dpi=100, facecolor='white')
            
            opt_sim_key = cluster_record['sim_key']
            opt_group_key = cluster_record['group_key']
            
            opt_cc_trace = self.calib_data[opt_group_key][opt_sim_key]['cc']
            
            ccs_model_agg, ccs_ref_agg = cc_data_aggregate(opt_cc_trace, cluster_id)
            
            #debug_p('model length ' + str(len(ccs_model_agg)))
            #debug_p('ref length ' + str(len(ccs_ref_agg)))
            
            hfca_id = cluster_id.split('_')[0]
            
            facility = hfca_id_2_facility(hfca_id)
            
            gs = gridspec.GridSpec(1, 4)
            ymax = 16
        
            scale_int = np.array(range(0,ymax+1))
            pal = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=ymax+1)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pal)
                         
            ax = plt.subplot(gs[0:4])
            
            #ax.set_ylim(1000)
            
            opt_const_h = cluster_record['habs']['const_h']
            opt_x_temp_h = cluster_record['habs']['temp_h']
            opt_itn = cluster_record['ITN_cov']
            opt_drug = cluster_record['MSAT_cov']
            opt_fit_value = cluster_record['fit_value']
            
            # the following code only relevant for rank correlation cc penalty fit
            opt_rho = None
            opt_p_val = None
            if 'rho' in cluster_record:
                opt_rho = cluster_record['rho']
            if 'p_val' in cluster_record:
                opt_p_val = cluster_record['p_val']
            
            
            mod_dates, mod_cases = zip(*ccs_model_agg)
            ref_dates, ref_cases = zip(*ccs_ref_agg)
            
            if opt_rho and opt_p_val:
                ax.plot(mod_dates, mod_cases, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn) + ', rho=' + str(opt_rho) + ', p-val=' + str(opt_p_val), marker = opt_marker, markersize = opt_marker_size)
            else:
                ax.plot(mod_dates, mod_cases, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn), marker = opt_marker, markersize = opt_marker_size)
            ax.bar(ref_dates, ref_cases, width=12,color='red',edgecolor='red', linewidth=0, label = 'Observed in ' + facility)
            #ax.plot(dates, ccs_ref_agg, alpha=1, linewidth=2.0, c = 'red', label = 'Observed in ' + facility)
            
            count_traces = 0 
            for fit_entry in self.all_fits[cluster_id]:
                
                sim_key = cluster_record['sim_key']
                group_key = cluster_record['group_key']
                fit_val = fit_entry['fit_value']
            
                if sim_key == opt_sim_key and fit_val > opt_fit_value + opt_fit_value*subopt_plots_threshold: 
                # do not plot optimal traces since we've already plotted it ;also do not plot too suboptimal traces
                    continue
            
                cc_trace = self.calib_data[group_key][sim_key]['cc']
            
                ccs_model_agg, ccs_ref_agg = cc_data_aggregate(cc_trace, cluster_id)
                
                # the following code only relevant for rank correlation cc penalty fit
                rho = None
                p_val = None
                if 'rho' in cluster_record:
                    rho = cluster_record['rho']
                if 'p_val' in cluster_record:
                    p_val = cluster_record['p_val']
                    
                
                const_h = fit_entry['const_h']
                x_temp_h = fit_entry['temp_h']
                itn = fit_entry['ITN_cov']
                drug = fit_entry['MSAT_cov']
                
                
                mod_dates, mod_cases = zip(*ccs_model_agg)
                ref_dates, ref_cases = zip(*ccs_ref_agg)
                
                if not sim_key in self.fit_entries_2_markers:
                    marker = markers[count_traces % len(markers)]
                    self.fit_entries_2_markers[sim_key] = marker
                else:
                    marker = self.fit_entries_2_markers[sim_key]
                    
                if rho and p_val:
                    ax.plot(mod_dates, mod_cases, alpha=0.75, linewidth=2.0, marker = marker, label = 'eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h) + 'rho=' + str(rho) + ', p-val=' + str(p_val))
                else:
                    ax.plot(mod_dates, mod_cases, alpha=0.75, linewidth=2.0, marker = marker, label = 'eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h)) 
                ax.bar(ref_dates, ref_cases, width=12,color='red',edgecolor='red', linewidth=0, label = 'Observed in ' + facility)
                
            
            plt.xlabel('Time (6-week bins)', fontsize=8)
            plt.ylabel('Clinical cases', fontsize=8)
            plt.legend(loc=1, fontsize=8)
            plt.title('Clinical cases timeseries', fontsize = 8, fontweight = 'bold', color = 'black')
            plt.gca().tick_params(axis='x', labelsize=8)
            plt.gca().tick_params(axis='y', labelsize=8)
            plt.tight_layout()
            output_plot_file_path = os.path.join(self.root_sweep_dir, cc_subopt_traces_plots_dir, cc_traces_base_file_name + cluster_id + '.png')
            plt.savefig(output_plot_file_path, format='png')
            plt.close()