#!/usr/bin/env python2.7

# Beta test-suite (see ticket #65)
# Code for calculating and plotting the bit comparison tests 
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
# Adhill - 040516

import netCDF4
import numpy as np
import matplotlib.pyplot as p
import os
import sys
import subprocess as sub
#from prettytable import PrettyTable

#import local module
import get_data

def monc_calc(testcase, test_run, model_dict) :
    
    Bit_compare_threshold = 0.0

    # list of 3d fields to compare
    fields_3d = ['u', 'v', 'w', 'q_vapour', 'q_cloud_liquid_mass']

    # Work out the known good output (kgo) difference
    #
    # construct filename using model_dict information
    # KGO first
    new_file_test = True
    kgo_file_test = True
    kgo_fn_1pe =  model_dict['kgo_dir']+test_run+model_dict['no_procs'][0]+model_dict['files'][0]+'.nc'
    print kgo_fn_1pe
    if (not os.path.exists(kgo_fn_1pe)) :
        kgo_file_test = False
    kgo_fn_32pe = model_dict['kgo_dir']+test_run+model_dict['no_procs'][1]+model_dict['files'][0]+'.nc'
    print kgo_fn_32pe
    if (not os.path.exists(kgo_fn_32pe)) :
        kgo_file_test = False
    # New ouput 
    new_fn_1pe =  model_dict['new_dir']+test_run+model_dict['no_procs'][0]+model_dict['files'][0]+'.nc'
    if (not os.path.exists(new_fn_1pe)) :
        new_file_test = False
    new_fn_32pe = model_dict['new_dir']+test_run+model_dict['no_procs'][1]+model_dict['files'][0]+'.nc'
    print new_fn_32pe
    if (not os.path.exists(new_fn_32pe)) :
        new_file_test = False
    # check files exist
    if not new_file_test :
        with open('bit_compare_results.txt', "a") as f:
            f.write(test_run+' bit compare **FAIL** - file missing, continue to next test'+"\n")
        return True, kgo_file_test  
    else : # files exist, so carry on with bit comparison
        # open both KGO files (if they exist)
        # script will still do bit compare test on new files even if
        # the kgo does not exist
        if kgo_file_test :
            kgo_1pe_open=netCDF4.Dataset(kgo_fn_1pe,'r') 
            kgo_32pe_open=netCDF4.Dataset(kgo_fn_32pe,'r')
 
        # open both NEW test files
        new_1pe_open=netCDF4.Dataset(new_fn_1pe,'r') 
        new_32pe_open=netCDF4.Dataset(new_fn_32pe,'r')
     
        decomp_success = True
        # Work out the difference the new output on different decompositions 
        for field in fields_3d :
            
            if kgo_file_test :
                # Extract field from the KGO file
                kgo_w_1pe = kgo_1pe_open.variables[field][:,:,:,:]
                kgo_w_32pe = kgo_32pe_open.variables[field][:,:,:,:]

            # Extract field from the New output
            new_w_1pe = new_1pe_open.variables[field][:,:,:,:]
            new_w_32pe = new_32pe_open.variables[field][:,:,:,:]
        
            # work out the difference between the processor decompositions
            # KGO
            if kgo_file_test :
                kgo_diff_pe = kgo_w_1pe - kgo_w_32pe
                if (np.max(abs(kgo_diff_pe)) > Bit_compare_threshold) :
                    decomp_success = False
            # NEW 
            new_diff_pe = new_w_1pe - new_w_32pe 
            
            if (np.max(abs(new_diff_pe)) > Bit_compare_threshold) :
                if kgo_file_test :          
                    if decomp_success :
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' **FAIL** (KGO passed)'+"\n")
                            f.write('Field ='+field+' - Max new decomp diff ='+str(np.max(abs(new_diff_pe)))+"\n")
                        decomp_success = False
                    else :
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' **FAIL** (KGO Failed)'+"\n")
                            f.write('Field ='+field+' - Max KGO decomp diff ='+str(np.max(abs(kgo_diff_pe)))+"\n")
                            f.write('Field ='+field+' - Max new decomp diff ='+str(np.max(abs(new_diff_pe)))+"\n")
                else :
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' **FAIL** (NO KGO files to test)'+"\n")
                        f.write('Field ='+field+' - Max new decomp diff ='+str(np.max(abs(new_diff_pe)))+"\n")
                    decomp_success = False                            
            else :
                if kgo_file_test : 
                    if decomp_success :
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS (KGO passed)'+"\n")
                    else :
                        with open('bit_compare_results.txt', "a") as f:
                            f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS (KGO failed)'+"\n")
                else :
                    with open('bit_compare_results.txt', "a") as f:
                        f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS (NO KGO files to test)'+"\n")

            if kgo_file_test :
                diff_kgo_new_1pe = new_w_1pe - kgo_w_1pe                                                                                                    
                diff_kgo_new_32pe = new_w_32pe - kgo_w_32pe                                                                                                 
                                                                                                                                                        
                if (np.max(abs(diff_kgo_new_32pe)) > Bit_compare_threshold) :                                                                               
                    decomp_success = False                                                                                                                  
                    print "KGO problem 32pe", decomp_success, field                                                                                         
                elif (np.max(abs(diff_kgo_new_1pe)) > Bit_compare_threshold ) :                                                                             
                    decomp_success = False                                                                                                                  
                    print "KGO problem 1pe", decomp_success, field                                                                                          
                                                                                                                                                        
                if (np.max(abs(diff_kgo_new_1pe)) > Bit_compare_threshold ) :                                                                               
                    with open('bit_compare_results.txt', "a") as f:                                                                                         
                        f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs kgo '+model_dict['no_procs'][0]+'pe '+test_run+' **FAIL**'+"\n")  
                        f.write('Field ='+field+' - Max 1pe diff ='+str(np.max(abs(diff_kgo_new_1pe)))+"\n")                                                 
                    decomp_success = False                                                                                                                  
                else :                                                                                                                                      
                    with open('bit_compare_results.txt', "a") as f:                                                                                         
                        f.write('Field ='+field+' - New '+model_dict['no_procs'][0]+'pe vs kgo '+model_dict['no_procs'][0]+'pe '+test_run+' SUCCESS'+"\n")  
                                                                                                                                                        
                    if (np.max(abs(diff_kgo_new_32pe)) > Bit_compare_threshold ) :                                                                              
                        with open('bit_compare_results.txt', "a") as f:                                                                                         
                            f.write('Field ='+field+' - New '+model_dict['no_procs'][1]+'pe vs kgo '+model_dict['no_procs'][1]+'pe '+test_run+' **FAIL**'+"\n") 
                            f.write('Field ='+field+' - Max 32pe diff ='+str(np.max(abs(diff_kgo_new_32pe)))+"\n")                                              
                        decomp_success = False                                                                                                                  
                    else :                                                                                                                                      
                        with open('bit_compare_results.txt', "a") as f:                                                                                         
                            f.write('Field ='+field+' - New '+model_dict['no_procs'][1]+'pe vs kgo '+model_dict['no_procs'][1]+'pe '+test_run+' SUCCESS'+"\n")  
                with open('bit_compare_results.txt', "a") as f:                                                                                             
                    f.write('..........................................................................................'+"\n")
            else:
                with open('bit_compare_results.txt', "a") as f:
                    f.write('Field ='+field+' - New '+model_dict['no_procs'][1]+'pe vs kgo '+model_dict['no_procs'][1]+'pe '+test_run+' NO KGO DATA'+"\n") 
                    f.write('..........................................................................................'+"\n")
        if kgo_file_test:            
            kgo_1pe_open.close()
            kgo_32pe_open.close()

        new_1pe_open.close()
        new_32pe_open.close()    
        print 'end', decomp_success, kgo_file_test  
        return decomp_success, kgo_file_test  
      
def lem_calc(testcase, test_run, model_dict) :

    Bit_compare_threshold = 0.0

    for idx, file_name in enumerate(model_dict['files']) :
        # construct filename using model_dict information
        model_fn_1pe =  model_dict['kgo_dir']+test_run+model_dict['no_procs'][0]+file_name+'.nc'
        model_fn_32pe = model_dict['kgo_dir']+test_run+model_dict['no_procs'][1]+file_name+'.nc'
            
        # open both files
        f_1pe_open=netCDF4.Dataset(model_fn_1pe,'r') 
        f_32pe_open=netCDF4.Dataset(model_fn_32pe,'r')

        w_1pe = f_1pe_open.variables['W'][:,:,:]
        w_32pe = f_32pe_open.variables['W'][:,:,:]
            
        diff_pe = w_1pe - w_32pe

        if (np.max(abs(diff_pe)) > Bit_compare_threshold) :
            print 'LEM '+test_run+' bit compare test also **FAIL**'
            return False
        else :
            if idx == len(model_dict['files'])-1 :
                print 'LEM '+test_run+' bit compare test SUCCESS'
                return True

def monc_plot(testing, testcase, test_run, fig_title, monc_dict):
    
    #
    # plot the bit comparison results first
    #
    if testing == 'bit_compare':
        model_1Pe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][0]
        model_MultiPe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][1]
        label_prefix_1Pe = 'MONC 1 pe '
        label_prefix_32Pe = 'MONC 32 pe '
    elif testing == 'kgo_compare_1pe':
        model_1Pe_fn = monc_dict['kgo_dir']+test_run+monc_dict['no_procs'][0]
        model_MultiPe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][0]
        label_prefix_1Pe = 'KGO 1 pe '
        label_prefix_32Pe = 'NEW 1 pe '
    elif testing == 'kgo_compare_32pe':    
        model_1Pe_fn = monc_dict['kgo_dir']+test_run+monc_dict['no_procs'][1]
        label_prefix_1Pe = 'KGO 32 pe '
        label_prefix_32Pe = 'NEW 32 pe '
        model_MultiPe_fn = monc_dict['new_dir']+test_run+monc_dict['no_procs'][1]
        

    # this dictionary is used to output data from get_profile 
    #it is overwritten on each loop
    monc_1Pe_data = { 'zn': [],'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                      'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],
                      'vapour_mmr_mean': [], 'liquid_mmr_mean': []}
    print model_1Pe_fn
    get_data.profile_diagnostics(model_1Pe_fn, monc_dict['files'], monc_1Pe_data)

    # this dictionary is used to output data from get_profile 
    #it is overwritten on each loop
    monc_MultiPe_data = { 'zn': [],'ww_mean':[],'uu_mean':[], 'vv_mean':[], 
                          'theta_mean':[], 'u_mean':[], 'v_mean': [], 'w_mean': [],
                          'vapour_mmr_mean': [], 'liquid_mmr_mean': []}
    print model_MultiPe_fn
    get_data.profile_diagnostics(model_MultiPe_fn, monc_dict['files'], monc_MultiPe_data)

    z_plot =  monc_MultiPe_data['zn'][0][:]

    # Set up the velocity figure for the case 
    #
    velocity_fig = p.figure(figsize=(15,10))
    velocity_fig.suptitle(fig_title, fontsize=16 )
    ax1 = velocity_fig.add_subplot(231)
    ax2 = velocity_fig.add_subplot(232)
    ax3 = velocity_fig.add_subplot(233)
    ax4 = velocity_fig.add_subplot(234)
    ax5 = velocity_fig.add_subplot(235)
    ax6 = velocity_fig.add_subplot(236)

    for file_index, file_numbers in enumerate(monc_dict['files']) :
        print file_index
        ax1.plot(monc_1Pe_data['u_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index],label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
        ax1.plot(monc_MultiPe_data['u_mean'][0][file_index,:], z_plot[:]
                 , '--',color=monc_dict['time_color'][file_index],label=label_prefix_32Pe+monc_dict['t_labels'][file_index] )
            
        ax2.plot(monc_1Pe_data['v_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
        ax2.plot(monc_MultiPe_data['v_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
            
        ax3.plot(monc_1Pe_data['w_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
        ax3.plot(monc_MultiPe_data['w_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
        
        ax4.plot(monc_1Pe_data['uu_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
        ax4.plot(monc_MultiPe_data['uu_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
        
        ax5.plot(monc_1Pe_data['vv_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
        ax5.plot(monc_MultiPe_data['vv_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
        
        ax6.plot(monc_1Pe_data['ww_mean'][0][file_index,:], z_plot[:],
                 color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
        ax6.plot(monc_MultiPe_data['ww_mean'][0][file_index,:], z_plot[:], '--',
                 color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
    
    ax1.set_title('Mean U')
    ax2.set_title('Mean V')
    ax3.set_title('Mean W')
    ax4.set_title('Mean UU')
    ax5.set_title('Mean VV')
    ax6.set_title('Mean WW')
    ax1.set_ylabel('height (m)')
    ax4.set_ylabel('height (m)')
    ax1.set_xlabel('m s$^{-1}$')
    ax2.set_xlabel('m s$^{-1}$')
    ax3.set_xlabel('m s$^{-1}$')
    ax4.set_xlabel('m$^{2}$ s$^{-2}$')
    ax5.set_xlabel('m$^{2}$ s$^{-2}$')
    ax6.set_xlabel('m$^{2}$ s$^{-2}$')
        
    lgd = ax3.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
    
    p.savefig(testing+'_'+test_run+'velocity_monc_bit_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
    p.close(velocity_fig)

    if testcase == 'bubble' :
        theta_fig = p.figure(figsize=(8,5))
        theta_fig.suptitle(fig_title, fontsize=16 )
        ax1 = theta_fig.add_subplot(111)
        for file_index, file_numbers in enumerate(monc_dict['files']) :
            ax1.plot(monc_1Pe_data['theta_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
            ax1.plot(monc_MultiPe_data['theta_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
        ax1.set_title('Mean potential temp')
        ax1.set_ylabel('height (m)')
        ax1.set_xlabel('potential T (K)')
        
        lgd = ax1.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
        p.savefig(test_run+'theta_monc_bit_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
        p.close(theta_fig)
        

    if testcase == 'stratus' or testcase == 'shallow_convection' or testcase == 'rce':
        q_theta_fig = p.figure(figsize=(15,5))
        q_theta_fig.suptitle(fig_title, fontsize=16 )
        ax1 = q_theta_fig.add_subplot(131)
        ax2 = q_theta_fig.add_subplot(132)
        ax3 = q_theta_fig.add_subplot(133)
        for file_index, file_numbers in enumerate(monc_dict['files']) :
            print file_index
            ax1.plot(monc_1Pe_data['theta_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
            ax1.plot(monc_MultiPe_data['theta_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
            ax2.plot(monc_1Pe_data['vapour_mmr_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
            ax2.plot(monc_MultiPe_data['vapour_mmr_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
            ax3.plot(monc_1Pe_data['liquid_mmr_mean'][0][file_index,:], z_plot[:],
                     color=monc_dict['time_color'][file_index], label=label_prefix_1Pe+monc_dict['t_labels'][file_index])
            ax3.plot(monc_MultiPe_data['liquid_mmr_mean'][0][file_index,:], z_plot[:], '--',
                     color=monc_dict['time_color'][file_index], label=label_prefix_32Pe+monc_dict['t_labels'][file_index])
                    
        ax1.set_title('Mean potential temp')
        ax1.set_ylabel('height (m)')
        ax1.set_xlabel('potential T (K)')
        ax2.set_title('Vapour mass mixing ratio')
        ax2.set_xlabel('qv (kg kg$^{-1}$')
        ax3.set_title('Liquid mass mixing ratio')
        ax3.set_xlabel('ql (kg kg$^{-1}$')

        lgd = ax3.legend(labelspacing=0.2, loc='center left', bbox_to_anchor=(1, 0.5))
        p.savefig(test_run+'q_theta_monc_bit_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight')
        p.close(q_theta_fig)
                     
     

    
