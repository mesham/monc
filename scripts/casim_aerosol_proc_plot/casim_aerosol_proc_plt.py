#!/usr/bin/env python2.7

# Beta test-suite (see ticket #65) 
# This has been developed to automate the component testing of MONC
# The attached script is an example of funtionality required for the development of a 
# Rose suite to automate the testing of the MONC. The script below is specific to MONSOON
# The script below is used to analyse the output from monc_test_harness.py and compare with 
# the LEM, where appropriate.
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
# Adhill - 040516

import sys
import os
import getopt

#import local module
import bit_compare
import dictionary_setup
import list_setup
import test_harness_setup

suite='casim_aerosol_processing'
print 'casim_socrates_processing'
procs = [72, 36]
compare_trunk_pkg = False
testcases = [ 'stratus_hilletal', 'lem_bomex', 'rce_deep' ]
        
# test_cfg={ 
#     'stratus_hilletal':[ 'ScFull_2M_fullproc_iopt3',  'ScFull_2M_fullproc_iopt5',  'ScFull_2M_NdFix',
#                          'ScFull_2M_noproc_iopt3',  'ScFull_2M_noproc_iopt5',  'ScFull_2M_passiveproc_iopt3',
#                          'ScFull_2M_passiveproc_iopt5' ],
#     'lem_bomex':['CuNoDamp_2M_fullproc_iopt3',  'CuNoDamp_2M_fullproc_iopt5',  'CuNoDamp_2M_NdFix',
#                  'CuNoDamp_2M_noproc_iopt3',  'CuNoDamp_2M_noproc_iopt5',  'CuNoDamp_2M_passiveproc_iopt3',
#                  'CuNoDamp_2M_passiveproc_iopt5'],
#     'rce_deep':['RCENoDamp_2M_fullproc_iopt3',  'RCENoDamp_2M_fullproc_iopt5',  
#                 'RCENoDamp_2M_noproc_iopt3',  'RCENoDamp_2M_noproc_iopt5',  'RCENoDamp_2M_passiveproc_iopt3',
#                 'RCENoDamp_2M_passiveproc_iopt5']
# }
test_cfg_arg={ 
    'stratus_hilletal':['ScFull_2M_fullproc_iopt3','ScFull_2M_noproc_iopt3', 'ScFull_2M_passiveproc_iopt3']
    'lem_bomex':['CuNoDamp_2M_fullproc_iopt3', 'CuNoDamp_2M_noproc_iopt3', 'CuNoDamp_2M_passiveproc_iopt3'],
    'rce_deep':['RCENoDamp_2M_fullproc_iopt3','RCENoDamp_2M_noproc_iopt3', 'RCENoDamp_2M_passiveproc_iopt3',]
}

if compare_trunk_pkg :
    monc_kgo_datadir = '/data/d04/fra23/monc_branches/monc_trunk_r8166/test_harness/MO-XC40_casim_aerosol_processing_cray_safe/'
    new_monc_datadir = '/data/d04/fra23/monc_branches/r8166_monc_w_casim_vn0p5/test_harness/MO-XC40_casim_aerosol_processing_cray_safe/'
    #monc_casim_gordon2020 = '/data/d04/fra23/monc_branches/r8725_demistify_aero_test_processing/test_harness/MO-XC40_casim_aerosol_processing_cray_safe/'
else :
    monc_kgo_datadir = '/data/d04/fra23/monc_branches/r8166_monc_w_casim_vn0p5/test_harness/MO-XC40_casim_aerosol_processing_cray_safe/'
    new_monc_datadir = '/data/d04/fra23/monc_branches/r8725_demistify_aero_test_processing/test_harness/MO-XC40_casim_aerosol_processing_cray_safe/'

# read in lists and dictionaries, which are
# required to setup and run the test_harness
dict_list = test_harness_setup.dict_list_read(suite)
procs = dict_list[0]
test_case = dict_list[1]

with open('bit_compare_results.txt', "w") as f:
    f.write('****************Bit comparison results from MONC-MAIN component testing********************'+"\n")
    f.write('Known good MONC output = '+monc_casim_kgo+"\n")
    f.write('New MONC output = '+monc_casim_test+"\n")
    
for case in test_case :
    # set up the dictionaries for each case, so diags names can be linked to dir and file nameare linked to 
    if case == 'stratus_hilletal':
        monc_dict = { 'kgo_dir': kgo_datadir,
                      'new_dir': new_datadir,
                      'files':['diagnostics_7260.nc', 'diagnostics_7260.nc', 'diagnostics_7260.nc',
                              'diagnostics_7260.nc'],
                      't_labels':['30 mins', '60 mins', '90 mins', '120 min'],
                      'time_color':[ 'blue', 'red', 'y', 'c' ],
                      'no_procs':['72_']} 
    elif case == 'lem_bomex':
        monc_dict = {'kgo_dir': kgo_datadir,
                     'new_dir': new_datadir,
                     'files':['diagnostics_21600.nc', 'diagnostics_21600.nc', 'diagnostics_21600.nc',
                              'diagnostics_21600.nc'],
                     't_labels':['210 mins', '240 mins', '270 mins', '300 mins', '330 mins', 
                                 '360 mins'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y' ],
                     'no_procs':['72_']}
    elif case == 'rce_deep':
        monc_dict = { 'kgo_dir': kgo_datadir,
                      'new_dir': new_datadir,
                      'files':['diagnostics_86400.nc', 'diagnostics_86400.nc', 'diagnostics_86400.nc'],
                     't_labels':['22 hours', '23 hours', '24 hours'],
                     'time_color':[ 'm', 'g', 'k'],
                     'no_procs':['72_']}   
    else :
        sys.exit('testcase is not defined. Please check testcase name, EXITING')

    with open('bit_compare_results.txt', "a") as f:
        f.write('********************************************************************************************'+"\n")
        f.write('******************************************'+case+'******************************************'+"\n")

    # do bit comparison check on MONC and maybe LEM
    for test_idx, test_run in enumerate(test_list) :
        print test_run
        with open('bit_compare_results.txt', "a") as f:
            f.write('=========================================================================================='+"\n")
        bit_compare_success = bit_compare.monc_calc(case, test_run, monc_dict)
        
        if not bit_compare_success[0] :  
            print 'plot bit comparison from monc to assess issue'
            bit_compare.monc_plot('bit_compare', case, test_run, figure_title_list[test_idx], monc_dict) 

            if bit_compare_success[1] :
                bit_compare.monc_plot('kgo_compare_1pe', case, test_run, figure_title_list[test_idx], monc_dict)  
                bit_compare.monc_plot('kgo_compare_32pe', case, test_run, figure_title_list[test_idx], monc_dict)
                


