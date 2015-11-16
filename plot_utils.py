import os
import sys

import numpy as np
import pandas as pd
import math
import json

from scipy import interpolate
from scipy.interpolate import Rbf
from scipy.interpolate import interp2d
from scipy.interpolate import spline
from scipy.stats import scoreatpercentile as satp

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

from kariba_settings import opt_marker, opt_marker_size, markers, subopt_plots_threshold, cc_penalty_model, cc_agg_fold, hfca_id_2_facility, cluster_2_prevs as c2p, traces_plots_dir, traces_base_file_name, cc_traces_plots_dir, cc_traces_base_file_name, err_surfaces_plots_dir, err_surfaces_base_file_name, sim_data_dir, calibration_data_file, tags_data_file, channels_sample_points, objectives_channel_codes, reports_channels, channels, cc_sim_start_date, cc_ref_start_date, cc_ref_end_date, cc_penalty_model,\
    cc_weight, hfca_id_2_cluster_ids, sample_size, cc_num_fold_bins, weighted_cc_traces_base_file_name, weighted_cc_traces_plots_dir, sample_size_percentile, weighted_ccs_by_hfca_id_file, get_cluster_category
from kariba_utils import get_cc_model_ref_traces, get_cc_penalty, sim_key_2_group_key

class PlotUtils():
    
    def __init__(self, best_fits, all_fits, calib_data, residuals, root_sweep_dir, category):
        
        self.best_fits = best_fits
        self.all_fits = all_fits
        self.residuals = residuals
        #debug_p(self.residuals)
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
            
            debug_p('Plotting prevalence trace for cluster ' + cluster_id + ' in category ' + self.category)
        
            #fig = plt.figure(cluster, figsize=(11, 4), dpi=100, facecolor='white')
            fig = plt.figure(cluster_id, figsize=(9.2, 4), dpi=100, facecolor='white')
            
            gs = gridspec.GridSpec(1, 2)
            ymax = 16
        
            scale_int = np.array(range(0,ymax+1))
            pal = cm = plt.get_cmap('jet') 
            cNorm  = colors.Normalize(vmin=0, vmax=ymax+1)
            scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pal)
                         
            ax = plt.subplot(gs[0:2])

            ax.set_xlim(2150,3010)    
            
            ref = self.get_ref(cluster_id)
            
            if not ref:
                # no reference object could be constructed; no data found; return w/o plotting
                return
            
            
            with open(os.path.join(sim_data_dir,weighted_ccs_by_hfca_id_file), 'r') as wccs_f:
                cluster_ccs_samples = json.load(wccs_f)
            
            hfca_id = cluster_id.split('_')[0]
            
            count_labels = 0
            for (group_key, sim_key, ccs) in cluster_ccs_samples[hfca_id][cluster_id]['unweighted']:
                 prev_trace = self.calib_data[group_key][sim_key]['prevalence']
                 if count_labels == 0:
                     ax.plot(range(2150, 3000), prev_trace[2149:2999], alpha=0.35, linewidth=0.5, color = 'magenta', label = 'Opt 5-percentile samples for cluster ' + cluster_id)
                     count_labels += 1
                 else:
                    ax.plot(range(2150, 3000), prev_trace[2149:2999], alpha=0.35, linewidth=0.5, color = 'magenta', marker = None)



            opt_sim_key = cluster_record['sim_key']
            opt_group_key = cluster_record['group_key']
            
            opt_const_h = cluster_record['habs']['const_h']
            opt_x_temp_h = cluster_record['habs']['temp_h']
            opt_itn = cluster_record['ITN_cov']
            opt_drug = cluster_record['MSAT_cov']
            opt_fit_value = cluster_record['fit_value']
            
            opt_prev_trace = self.calib_data[opt_group_key][opt_sim_key]['prevalence']
            
            ax.plot(range(2150, 3000), opt_prev_trace[2149:2999], alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: prevalence + clinical cases: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn))
            
            
            sim_key = cluster_record['mse']['sim_key']
            # avoid doing that and add group_key to corresponding terms in best_fits
            group_key = sim_key_2_group_key(sim_key)

            const_h = cluster_record['mse']['const_h']
            x_temp_h = cluster_record['mse']['temp_h']
            itn = cluster_record['mse']['ITN_cov']
            drug = cluster_record['mse']['MSAT_cov']
            
            prev_trace_by_prev = self.calib_data[group_key][sim_key]['prevalence']
            
            ax.plot(range(2150, 3000), prev_trace_by_prev[2149:2999], alpha=1, linewidth=2.0, c = 'magenta', label = 'Best fit: prevalence: eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h) + ', drug cov.' + str(drug) + ', ITN dist. = '+str(itn))
            
            
            sim_key = cluster_record['cc_penalty']['sim_key']
            group_key = sim_key_2_group_key(sim_key)
            
            const_h = cluster_record['cc_penalty']['const_h']
            x_temp_h = cluster_record['cc_penalty']['temp_h']
            itn = cluster_record['cc_penalty']['ITN_cov']
            drug = cluster_record['cc_penalty']['MSAT_cov']
            
            prev_trace_by_cc = self.calib_data[group_key][sim_key]['prevalence']
            
            ax.plot(range(2150, 3000), prev_trace_by_cc[2149:2999], alpha=1, linewidth=2.0, c = 'blue', label = 'Best fit: clinical cases: eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h) + ', drug cov.' + str(drug) + ', ITN dist. = '+str(itn))

            
            obs_prevs = ref.to_dict()['prevalence']['d_points']
            
            label_obs_prev_shown = False
            label_sim_prev_shown = False
            max_obs_prev = 0.0


            for i,prev in enumerate(obs_prevs):
                
                if prev != 'nan':
                    if max_obs_prev < prev:
                        max_obs_prev = prev
                        
                    if not label_obs_prev_shown:
                        label_obs_prev = 'Observed prevalence'
                        label_obs_prev_shown = True
                    else: 
                        label_obs_prev = None
                         
                    ax.scatter(channels_sample_points['prevalence'][i], prev, c = 'red', facecolor = 'red', marker='o', s = 40, label = label_obs_prev, zorder=200)
                    
                if not label_sim_prev_shown:
                    label_sim_prev = 'Best fit simulated prevalence at rnds.'
                    label_sim_prev_shown = True
                else: 
                    label_sim_prev = None
                
                ax.scatter(channels_sample_points['prevalence'][i], opt_prev_trace[channels_sample_points['prevalence'][i]], c = 'black', facecolor = 'none', marker='o', s = 60, label = label_sim_prev, zorder=150)
                
                
            '''
            count_traces = 0 
            for sim_key,fit_entry in self.all_fits[cluster_id].iteritems():

                if sim_key == 'min_terms' or sim_key == 'max_terms':
                    continue
                
                group_key = fit_entry['group_key']
                fit_val = fit_entry['fit_val']
            
                const_h = fit_entry['const_h']  
                x_temp_h = fit_entry['x_temp_h']
            
                #if sim_key == opt_sim_key or fit_val > opt_fit_value + opt_fit_value*subopt_plots_threshold or count_traces > 10:
                #if fit_entry['fit_terms']['mse'] <= satp(z, sample_size_percentile):
                #do not plot optimal traces since we've already plotted it ;also do not plot too many suboptimal traces
                #    continue
                
                prev_trace = self.calib_data[group_key][sim_key]['prevalence']
                marker = self.get_marker(sim_key, count_traces)
                #ax.plot(range(2180, 3000), prev_trace[2179:2999], alpha=0.75, linewidth=0.5,  marker = marker, markersize = 0.5*opt_marker_size, label = 'eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h))
                #ax.plot(range(2180, 3000), prev_trace[2179:2999], alpha=0.75, linewidth=0.5,  marker = marker, markersize = 0.5*opt_marker_size)
                ax.plot(range(2150, 3000), prev_trace[2149:2999], alpha=0.75, linewidth=0.5)
                
                
                for i,prev in enumerate(obs_prevs):                    
                    ax.scatter(channels_sample_points['prevalence'][i], prev_trace[channels_sample_points['prevalence'][i]], marker = marker, c = 'black', facecolor = 'none', s = 30)
                
                count_traces = count_traces + 1 
            '''
            
            
            ax.set_ylim(0,min(max(max_obs_prev, max(opt_prev_trace))+0.1,1))            
            plt.xlabel('Time (days)', fontsize=8)
            plt.ylabel('Prevalence (population fraction)', fontsize=8)
            plt.legend(loc=1, fontsize=8)
            plt.title('Prevalence timeseries', fontsize = 8, fontweight = 'bold', color = 'black')
            plt.gca().tick_params(axis='x', labelsize=8)
            plt.gca().tick_params(axis='y', labelsize=8)
        
            plt.tight_layout()
            output_plot_file_path = os.path.join(self.root_sweep_dir, traces_plots_dir, traces_base_file_name + cluster_id + '.png')
            plt.savefig(output_plot_file_path, dpi = 300, format='png')
            plt.close()
    
            count = count + 1
    
    
    # err_surface_types specifies which residual surfaces to plot (e.g. clinical cases penalty, prevalence mse, or combined residuals)
    # err_surface_type is a dictionary: {err_surface_type:'plot_title'}
    def plot_calib_err_surfaces(self, err_surface_types):
        
        count = 0
        
        min_residual = self.residuals['min']
        max_residual = self.residuals['max']
        
        for cluster_id, cluster_record in self.best_fits.iteritems():
            
            debug_p('Plotting error surface for cluster ' + cluster_id + ' in category ' + self.category)
        
            fig_width = len(err_surface_types)*4.35
            fig = plt.figure(cluster_id, figsize=(fig_width, 4), dpi=300, facecolor='white')
            #debug_p('error surface types length' + str(len(err_surface_types)))
            #debug_p('fig width' + str(fig_width))
            
            gs = gridspec.GridSpec(1, 3)
         
            error_points = {}
            
            title_ITN = 'ITN distribution: '
            title_drug_coverage = 'drug coverage: '
            
            opt_fit = {}
            
            opt_sim_key = cluster_record['sim_key']
            opt_group_key = cluster_record['group_key']
            opt_const_h = cluster_record['habs']['const_h']
            opt_x_temp_h = cluster_record['habs']['temp_h']
            
            for err_surface_type in err_surface_types:    
                opt_fit[err_surface_type] = {} 
                opt_fit[err_surface_type]['const_h'] = cluster_record[err_surface_type]['const_h']
                opt_fit[err_surface_type]['temp_h'] = cluster_record[err_surface_type]['temp_h']
                opt_fit[err_surface_type]['value'] = cluster_record[err_surface_type]['value']
                
                if err_surface_type == 'cc_penalty':
                    opt_fit[err_surface_type]['value'] = opt_fit[err_surface_type]['value']*(math.pow(cc_weight, -1)) # opt_fit of penalties contains the weighted value; hence we reverse the weighting
                    
                 
            opt_neigh_fits = []
            
            for sim_key,fit_entry in self.all_fits[cluster_id].iteritems():
        
                if sim_key == 'min_terms' or sim_key == 'max_terms':
                    continue
        
                x_temp_h = fit_entry['x_temp_h']
                const_h = fit_entry['const_h']
                fit_val = fit_entry['fit_val']
                mse = fit_entry['fit_terms']['mse']
                cc_penalty = get_cc_penalty(fit_entry)

                itn_level = fit_entry['itn_level']
                drug_coverage_level = fit_entry['drug_cov']
                group_key = fit_entry['group_key']    
                
                    
                if group_key not in error_points:
                    error_points[group_key] = {
                                                   'x_temp_h':[],
                                                   'const_h':[],
                                                   'fit':[],
                                                   'cc_penalty':[],
                                                   'mse':[],
                                                   'title': title_ITN + itn_level + "; " + title_drug_coverage + str(drug_coverage_level),
                                                   'itn_level':itn_level,
                                                   'drug_coverage':drug_coverage_level
                                                }
                    
                error_points[group_key]['x_temp_h'].append(x_temp_h)
                error_points[group_key]['const_h'].append(const_h)
                error_points[group_key]['fit'].append(fit_val)
                error_points[group_key]['mse'].append(mse)
                error_points[group_key]['cc_penalty'].append(cc_penalty)
                       
                          
            ymax = 10
            
            scale_int = np.array(range(1,10))
            pal = cm = plt.get_cmap('nipy_spectral') 
            cNorm  = colors.Normalize(vmin=0, vmax=ymax)
            #scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=pal)
            scalarMap = b2mpl.get_map('Spectral', 'Diverging', 5).mpl_colors
            lgds = []
            for i,(err_surface_type, err_surface_style) in enumerate(err_surface_types.iteritems()):
            
                for j,group_key in enumerate(error_points.keys()):
                    
                    itn_level = error_points[group_key]['itn_level']
                    drug_coverage = error_points[group_key]['drug_coverage']
                    
                    # currently assume opt_group_key is the same for all err_surface_types
                    if group_key == opt_group_key: 
            
                        #debug_p('plot at position (0, ' + str(i) + ') in grid')
                        plt.subplot(gs[0,i])
                        x = error_points[group_key]['x_temp_h']
                        y = error_points[group_key]['const_h']
                        z = error_points[group_key][err_surface_type]
                        #print len(z)
                        res = 125
                        ttl = err_surface_style['title']
                        
                        min_x = np.min(x)
                        min_y = np.min(y)
                        min_z = np.min(z)
                        
                        max_x = np.max(x)
                        max_y = np.max(y)
                        max_z = np.max(z)
                        
                        
                    
                        #f = interpolate.interp2d(x, y, z)
                    
                        xi = np.linspace(min_x, max_x , res)
                        yi = np.linspace(min_y, max_y , res)
                        
                        zi = griddata(x,y,z,xi,yi)
                        
                        #xig, yig = np.meshgrid(xi, yi)
                        #zig = f(xi,yi)
    
                        #rbf = Rbf(x, y, z, epsilon=2)
                        #zig = rbf(xig, yig)
                    
                        blevels = self.get_colorbar_ticks(min_z, max_z, z)
                        num_colors = len(blevels)-1
                        from matplotlib.colors import BoundaryNorm
                        
                        cmap2 = self.custom_cmap(num_colors, mincol='DarkBlue', midcol='CornflowerBlue', maxcol='w')
                        cmap2.set_over('0.7') # light gray
                        
                        bnorm = BoundaryNorm(blevels, ncolors = num_colors, clip = False)
                    
                        #rmse_pl = plt.contourf(xi,yi,zi,15,cmap=plt.cm.hot)
                        #rmse_pl = plt.pcolor(xi, yi, zi, cmap=plt.get_cmap('RdYlGn_r'), vmin=0.475, vmax=0.8)
                        #rmse_pl = plt.pcolor(xi, yi, zi, cmap=plt.get_cmap('RdYlGn_r'))
                        #rmse_pl = plt.pcolor(xi, yi, zi, cmap=plt.get_cmap('cool'), vmin = min_residual, vmax = max_residual)
                        #rmse_pl = plt.contourf(xi, yi, zi, cmap=plt.get_cmap('cool'), vmin = min_z, vmax = max_z, levels = blevels, norm = bnorm, extend = 'both')
                        rmse_pl = plt.contourf(xi, yi, zi, cmap=cmap2, vmin = min_z, vmax = max_z, levels = blevels, norm = bnorm, extend = 'both')
                        #rmse_pl = plt.contourf(xi,yi,zi,15,cmap=plt.get_cmap('paired'))
                        #rmse_pl.cmap.set_over('black')
                        #rmse_pl.cmap.set_under('grey')
                        
                        max_blevel_in_sample = 0
                        for blevel in blevels:
                            if blevel <= satp(z, sample_size_percentile) and blevel > max_blevel_in_sample:
                               max_blevel_in_sample = blevel 
                                
                        pc = plt.contour(xi,yi,zi, levels=[max_blevel_in_sample], colors='r', linewidth=0, alpha=0.5)
                        
                        b_per_s = pc.collections[0].get_paths()
                        count_labels = 0
                        for per in range(len(b_per_s)):
                            b_per_s_x = b_per_s[per].vertices[:,0]
                            b_per_s_y = b_per_s[per].vertices[:,1]
                            if count_labels == 0:
                                plt.fill(b_per_s_x,b_per_s_y, 'magenta', linestyle='solid', alpha=0.3, label = 'Opt 5-percentile: ' + err_surface_style['title'])
                                count_labels = count_labels + 1
                            else:
                                plt.fill(b_per_s_x,b_per_s_y, 'magenta', linestyle='solid', alpha=0.3)
                            
                        
                        cb = plt.colorbar(rmse_pl, ticks = blevels, spacing='uniform')
    
                        cb.set_label(ttl + ' residual', fontsize=8)
                        cb.ax.tick_params(labelsize=8)    
                        #plt.scatter(x, y, 10, z, cmap=cmap2,  vmin = min_z, vmax = max_z, norm = bnorm)
                        
                        level1_opt_neighs_label = False
                        level2_opt_neighs_label = False
                        level3_opt_neighs_label = False
                        # plot all optimal markers on each surface
                
                        for (opt_err_surface_type, opt_err_surface_style) in err_surface_types.iteritems():
                            plt.scatter(opt_fit[opt_err_surface_type]['temp_h'], opt_fit[opt_err_surface_type]['const_h'], c = 'red', marker = opt_err_surface_style['marker'], s = 60, facecolor='none', edgecolor='black', zorder=100, label= opt_err_surface_style['title'] + ' best fit')
                        
                        '''
                        for k,fit_val in enumerate(z):

                                if fit_val < opt_fit[err_surface_type]['value'] + opt_fit[err_surface_type]['value']*subopt_plots_threshold:
                                    
                                    if not level1_opt_neighs_label:
                                        label = '< opt + 0.1opt'
                                        level1_opt_neighs_label = True
                                    else:
                                        label = None
                                    
                                    #plt.scatter(x[k], y[k], 10, fit_val, marker = 'd',  linewidth = 0.75, color = 'green', label = label)
                                    
                                    
                                elif fit_val < opt_fit[err_surface_type]['value'] + 2*opt_fit[err_surface_type]['value']*subopt_plots_threshold:
                                    
                                    if not level2_opt_neighs_label:
                                        label = '< opt + 0.2opt'
                                        level2_opt_neighs_label = True
                                    else:
                                        label = None
        
                                    #plt.scatter(x[k], y[k], 10, fit_val, marker = 'o', linewidth = 0.75, color = 'blue', label = label)
                                    
                                    
                                elif fit_val < opt_fit[err_surface_type]['value'] + 3*opt_fit[err_surface_type]['value']*subopt_plots_threshold:
                                    
                                    if not level3_opt_neighs_label:
                                        label = '< opt + 0.3opt'
                                        level3_opt_neighs_label = True
                                    else:
                                        label = None                            
                                    
                                    #plt.scatter(x[k], y[k], 10, fit_val, marker = 'x',  linewidth = 0.75, color = 'red',  label = label)
                        '''
                            
                        #plt.title(ttl, fontsize = 8, fontweight = 'bold', color = 'white', backgroundcolor = scalarMap.to_rgba(scale_int[itn_levels_2_sbplts[itn_level]]))
                        plt.title(ttl, fontsize = 8, fontweight = 'bold', color = 'black')
                        plt.xlabel('All habitats scale', fontsize=8)
                        plt.ylabel('Constant habitat scale', fontsize=8)
                        plt.xlim(min_x+0.1, max_x+0.1)
                        plt.ylim(min_y+0.1, max_y+0.1)
                        #plt.ylim(0.01, 14)
                        plt.gca().tick_params(axis='x', labelsize=8)
                        plt.gca().tick_params(axis='y', labelsize=8)
                        
                        '''
                        count_traces = 0
                        
                        # NEED TO update to new FIT_ENTRY DATA STRUCT IF REUSED
                        for fit_entry in opt_neigh_fits:
                            
                            x_temp_h = fit_entry['x_temp_h']
                            const_h = fit_entry['const_h'] 
                            sim_key = fit_entry['sim_key']
                            
                            marker = self.get_marker(sim_key, count_traces)
        
                            plt.scatter(x_temp_h, const_h, c = 'black', marker = marker, s = 20, facecolor='none', zorder=100)
                            
                            count_traces = count_traces + 1
                        '''
                    
                #plt.subplot(gs[itn_levels_2_sbplts[best_fit_itn_level], 0])
                #debug_p('plot optimal at position (0, ' + str(i) + ') in grid')
                plt.subplot(gs[0,i])
                
                cluster_record = self.best_fits[cluster_id]
                opt_itn = cluster_record['ITN_cov']
                opt_drug = cluster_record['MSAT_cov']
                
                #plt.annotate(opt_fit_value, opt_x_temp_h, opt_const_h)
    
                #lgds.append(plt.legend(bbox_to_anchor=(0., 1, 1., .1), loc=2, ncol=1, borderaxespad=0., fontsize=8))
                lgds.append(plt.legend(ncol=1,loc='upper center', bbox_to_anchor=(0.,-0.15), borderaxespad=0., fontsize=8, mode='expand'))
                
                    
            plt.tight_layout()
            output_plot_file_path = os.path.join(self.root_sweep_dir, err_surfaces_plots_dir, err_surfaces_base_file_name + cluster_id +'.png')
            plt.savefig(output_plot_file_path, dpi = 300, format='png', bbox_extra_artists=lgds, bbox_inches='tight')
            plt.close()
            
            count = count + 1

                
    def plot_calib_cc_traces_clusters(self):
        
        for cluster_id, cluster_record in self.best_fits.iteritems():
            
            debug_p('Plotting clinical cases trace for cluster ' + cluster_id + ' in category ' + self.category)
            
            fig = plt.figure(cluster_id, figsize=(9.2, 4), dpi=100, facecolor='white')
            
            opt_sim_key = cluster_record['sim_key']
            opt_group_key = cluster_record['group_key']
            
            opt_cc_trace = self.calib_data[opt_group_key][opt_sim_key]['cc']
            
            ccs_model_agg, ccs_ref_agg = get_cc_model_ref_traces(opt_cc_trace, cluster_id)
            
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
            
            '''
            mod_dates, mod_cases = zip(*ccs_model_agg)
            ref_dates, ref_cases = zip(*ccs_ref_agg)
            '''
                
            if opt_rho and opt_p_val:
                ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn) + ', rho=' + str(opt_rho) + ', p-val=' + str(opt_p_val), marker = opt_marker, markersize = opt_marker_size)
            else:
                ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn), marker = opt_marker, markersize = opt_marker_size)
            ax.plot(range(0, len(ccs_ref_agg)), ccs_ref_agg, alpha=1, linewidth=2.0, c = 'red', label = 'Observed in ' + facility)    
                
            '''
            if opt_rho and opt_p_val:
                ax.plot(mod_dates, mod_cases, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn) + ', rho=' + str(opt_rho) + ', p-val=' + str(opt_p_val), marker = opt_marker, markersize = opt_marker_size)
            else:
                ax.plot(mod_dates, mod_cases, alpha=1, linewidth=2.0, c = 'black', label = 'Best fit: eff. constant=' + str(opt_const_h*opt_x_temp_h) + ', all='+str(opt_x_temp_h) + ', drug cov.' + str(opt_drug) + ', ITN dist. = '+str(opt_itn), marker = opt_marker, markersize = opt_marker_size)
            ax.bar(ref_dates, ref_cases, width=12,color='red',edgecolor='red', linewidth=0, label = 'Observed in ' + facility)
            #ax.plot(dates, ccs_ref_agg, alpha=1, linewidth=2.0, c = 'red', label = 'Observed in ' + facility)
            '''
            count_traces = 0 
            for sim_key, fit_entry in self.all_fits[cluster_id].iteritems():
                
                if sim_key == 'min_terms' or sim_key == 'max_terms':
                    continue
                
                group_key = fit_entry['group_key']
                fit_val = fit_entry['fit_val']
            
                if sim_key == opt_sim_key or fit_val > opt_fit_value + opt_fit_value*subopt_plots_threshold or count_traces > 10: 
                # do not plot optimal traces since we've already plotted it ;also do not plot too suboptimal traces
                    continue
                
                cc_trace = self.calib_data[group_key][sim_key]['cc']
            
                ccs_model_agg, ccs_ref_agg = get_cc_model_ref_traces(cc_trace, cluster_id)
                
                # the following code only relevant for rank correlation cc penalty fit
                rho = None
                p_val = None
                if 'rho' in fit_entry:
                    rho = fit_entry['rho']
                if 'p_val' in fit_entry:
                    p_val = fit_entry['p_val']
                    
                
                const_h = fit_entry['const_h']
                x_temp_h = fit_entry['x_temp_h']
                itn = fit_entry['itn_level']
                drug = fit_entry['drug_cov']
                
                '''
                mod_dates, mod_cases = zip(*ccs_model_agg)
                ref_dates, ref_cases = zip(*ccs_ref_agg)
                '''
                
                marker = self.get_marker(sim_key, count_traces)
                
                if rho and p_val:
                    #ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=0.75, linewidth=0.5,  marker = marker, markersize = 0.5*opt_marker_size, label = 'eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h) + 'rho=' + str(rho) + ', p-val=' + str(p_val))
                    #ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=0.75, linewidth=0.5,  marker = marker, markersize = 0.5*opt_marker_size)
                    ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=0.75, linewidth=0.5)
                else:
                    #ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=0.75, linewidth=0.5, marker = marker, markersize = 0.5*opt_marker_size, label = 'eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h))
                    ax.plot(range(0, len(ccs_model_agg)), ccs_model_agg, alpha=0.75, linewidth=0.5)
                #ax.plot(range(0, len(ccs_ref_agg)), ccs_ref_agg, alpha=1, linewidth=1.0, c = 'red', label = 'Observed in ' + facility)
                
                count_traces = count_traces + 1    
                
                '''    
                if rho and p_val:
                    ax.plot(mod_dates, mod_cases, alpha=0.75, linewidth=2.0, marker = marker, label = 'eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h) + 'rho=' + str(rho) + ', p-val=' + str(p_val))
                else:
                    ax.plot(mod_dates, mod_cases, alpha=0.75, linewidth=2.0, marker = marker, label = 'eff. constant=' + str(const_h*x_temp_h) + ', all='+str(x_temp_h)) 
                ax.bar(ref_dates, ref_cases, width=12,color='red',edgecolor='red', linewidth=0, label = 'Observed in ' + facility)
                '''
            
            plt.xlabel('6-week bins', fontsize=8)
            plt.ylabel('Clinical cases', fontsize=8)
            legend = plt.legend(loc=1, fontsize=8)
            
            '''
            init_font_size = 8
            for i,label in enumerate(legend.get_texts()):
                if i > 2:
                    label.set_fontsize(max(init_font_size - i, 5))
            '''

            plt.title('Clinical cases timeseries', fontsize = 8, fontweight = 'bold', color = 'black')
            plt.gca().tick_params(axis='x', labelsize=8)
            plt.gca().tick_params(axis='y', labelsize=8)
            plt.tight_layout()
            output_plot_file_path = os.path.join(self.root_sweep_dir, cc_traces_plots_dir, cc_traces_base_file_name + cluster_id + '.png')
            plt.savefig(output_plot_file_path, dpi = 300, format='png')
            plt.close()
            
    
    def plot_weighted_cc_per_hfca(self, weighted_ccs_model_agg_by_hfca, ccs_model_agg_by_hfca_cluster_id):
        
        
        clusters_processed = 0
        for hfca_id, weighted_ccs_combos in weighted_ccs_model_agg_by_hfca.iteritems():
        
            weighted_ccs_by_bin = {}
            for i in range(0, cc_num_fold_bins):
                weighted_ccs_by_bin[i] = []
                
            
            for weighted_ccs_combo in weighted_ccs_combos:
                sum_weighted_ccs = cc_num_fold_bins * [0]

                for weighted_ccs in weighted_ccs_combo:
                    sum_weighted_ccs = np.add(sum_weighted_ccs, weighted_ccs)
                    
                for i in range(0, cc_num_fold_bins):
                    weighted_ccs_by_bin[i].append(sum_weighted_ccs[i])
        
        
            per_bottom = []
            per_top = []
            per_median = []
            
            for i in range(0, cc_num_fold_bins):
                weighted_ccs_by_bin_idx = weighted_ccs_by_bin[i]
                per_bottom.append( satp(weighted_ccs_by_bin_idx, 2.5) )
                per_top.append( satp(weighted_ccs_by_bin_idx, 97.5) )
                per_median.append( satp(weighted_ccs_by_bin_idx, 50) )
                
            '''
            debug_p('length of weighted ccs_combos array ' + str(len(weighted_ccs_combos)))
            '''
            debug_p('length of bin 0 in weighted_ccs_by_bin ' + str(len(weighted_ccs_by_bin[0])))
           
            
            for cluster_id in hfca_id_2_cluster_ids(hfca_id):
                
                fig = plt.figure(cluster_id, figsize=(9.2, 4), dpi=100, facecolor='white')
                gs = gridspec.GridSpec(1, 4)
                     
                ax = plt.subplot(gs[0:4])

                x_smooth = np.linspace(0, cc_num_fold_bins-1,60)
                
                per_bottom_smooth = spline(range(0, cc_num_fold_bins),per_bottom,x_smooth)
                per_top_smooth = spline(range(0, cc_num_fold_bins),per_top,x_smooth)
                per_median_smooth = spline(range(0, cc_num_fold_bins),per_median,x_smooth)
                
                ax.plot(x_smooth, per_bottom_smooth, alpha=1, linewidth=0.5, color = 'black', linestyle=':', label = '2.5 percentile HS weighted: prevalence space samples', marker = None)
                ax.plot(x_smooth, per_top_smooth, alpha=1, linewidth=0.5, color = 'black', linestyle=':', label = '97.5 percentile HS weighted: prevalence space samples', marker = None)
                ax.plot(x_smooth, per_median_smooth, alpha=1, linewidth=2.0, color = 'magenta', linestyle='-', label = 'median HS weighted: prevalence space samples', marker = None)
                ax.fill_between(x_smooth, per_bottom_smooth, per_top_smooth, facecolor='gray', alpha=0.5, interpolate=True)
                
                cluster_cat = get_cluster_category(cluster_id)
                
                opt_group_key = self.best_fits[cluster_id]['group_key']
                
                opt_sim_key_cc = self.best_fits[cluster_id]['cc_penalty']['sim_key']
                cc_trace_opt_cc = self.calib_data[cluster_cat][opt_group_key][opt_sim_key_cc]
                
                opt_sim_key_prev = self.best_fits[cluster_id]['mse']['sim_key']
                cc_trace_opt_prev = self.calib_data[cluster_cat][opt_group_key][opt_sim_key_prev]
                
                opt_sim_key_fit = self.best_fits[cluster_id]['fit']['sim_key']
                cc_trace_opt_fit = self.calib_data[cluster_cat][opt_group_key][opt_sim_key_fit]
            
                ccs_model_agg_cc, ccs_ref_agg = get_cc_model_ref_traces(cc_trace_opt_cc, cluster_id)
                ccs_model_agg_prev, ccs_ref_agg = get_cc_model_ref_traces(cc_trace_opt_prev, cluster_id)
                ccs_model_agg_fit, ccs_ref_agg = get_cc_model_ref_traces(cc_trace_opt_fit, cluster_id)
                
                facility = hfca_id_2_facility(hfca_id)
                ax.plot(range(0, len(ccs_model_agg_cc)), ccs_model_agg_cc, alpha=1, linewidth=1, color = 'blue', label = 'Best fit: clinical cases', marker = 's')
                ax.plot(range(0, len(ccs_model_agg_prev)), ccs_model_agg_prev, alpha=1, linewidth=1, color = 'magenta', label = 'Best fit: prevalence', marker = 'o')
                ax.plot(range(0, len(ccs_model_agg_fit)), ccs_model_agg_fit, alpha=1, linewidth=1, color = 'black', label = 'Best fit: prevalence + clinical cases', marker = '*')
                ax.plot(range(0, len(ccs_ref_agg)), ccs_ref_agg, alpha=1, linewidth=2.0, linestyle = '-', color = 'red', label = 'Observed in ' + facility, marker = None)    
                
                for i,sample_ccs in enumerate(ccs_model_agg_by_hfca_cluster_id[hfca_id][cluster_id]['unweighted']):
                      
                    if i == 0:
                        ax.plot(range(0, cc_num_fold_bins), sample_ccs[2], alpha=0.5, linewidth=0.5, color = 'magenta', label = 'Opt 5-percentile samples for cluster ' + cluster_id, marker = None)
                        #ax.plot(range(0, cc_num_fold_bins), sample_ccs[2])
                    else:
                        ax.plot(range(0, cc_num_fold_bins), sample_ccs[2], alpha=0.5, linewidth=0.5, color = 'magenta', marker = None)
        
                plt.xlabel('6-week bins', fontsize=8)
                plt.ylabel('Clinical cases', fontsize=8)
                legend = plt.legend(loc=1, fontsize=8)
            
                
                plt.xlim(0,8)
                plt.title('Clinical cases timeseries', fontsize = 8, fontweight = 'bold', color = 'black')
                plt.gca().tick_params(axis='x', labelsize=8)
                plt.gca().tick_params(axis='y', labelsize=8)
                plt.tight_layout()
                output_plot_file_path = os.path.join(self.root_sweep_dir, weighted_cc_traces_plots_dir, weighted_cc_traces_base_file_name + cluster_id + '.png')
                plt.savefig(output_plot_file_path, dpi = 300, format='png')
                plt.close()
                
                clusters_processed = clusters_processed + 1
                
                debug_p('Processed weighting and plotting clinical cases for ' + str(clusters_processed) + ' clusters') 
            
    
    def get_marker(self, sim_key, count_traces):
        
        if not sim_key in self.fit_entries_2_markers:
            marker = markers[count_traces % len(markers)]
            self.fit_entries_2_markers[sim_key] = marker
        else:
            marker = self.fit_entries_2_markers[sim_key]
            
        return marker
    
    
    def get_colorbar_ticks(self, min_z, max_z, z):
        
        blevels = []
        ticks = [min_z]
        
        current = 0.1
        previous = 0.0
        index = 0
        
        self.generate_blevel_sequence(blevels, current, previous, index)
    
        for blevel in blevels:
            tick = min_z + blevel*min_z
            if tick <= max_z:
                ticks.append(tick)
        
        # not efficient; but there are many parts of this that aren't;
        # need to refactor at some point if speed is an issue
        # avoid premature optimization
        
        ticks.append(satp(z, sample_size_percentile))
        
                
        return sorted(ticks)
            
        
    def generate_blevel_sequence(self, blevels, current, previous, index):
        if index == 20:
            return current
        else:
            next = current + previous
            blevels.append(next)
            return self.generate_blevel_sequence(blevels, next, current, index + 1)
        
    
    def custom_cmap(self, numcolors=13, name='custom_cmap', mincol='blue', midcol='white', maxcol='black'):
        
        from matplotlib.colors import LinearSegmentedColormap 

        cmap = LinearSegmentedColormap.from_list(name=name, colors =[mincol, midcol, maxcol], N=numcolors)
        return cmap