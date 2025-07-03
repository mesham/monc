#!/usr/bin/env python2.7

# Beta test-suite (see ticket #65) 
# This has been developed to automate the component testing of MONC
# The attached script is an example of funtionality required for the development of a 
# Rose suite to automate the testing of the MONC. The script below is specific to MONSOON
# The script below is used to analyse the output from monc_test_suite.py and compare with 
# the LEM, where appropriate.
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
# Adhill - 040516

#import local module
import bit_compare
import dictionary_setup
import list_setup
import sys
import os

print 'monc_control_testsuite called'
       
test_case = [ 'stratus', 'stratus_diurnal', 'shallow_convection', 'rce']
optimisation_level = 3
checkpoint_test = True
compiler = 'gnu'
kgo_dir = 'latest_kgo/r5815_vn0.3casim'
branch_path= os.environ.get("BRANCH-PATH",
                            "please provide the full path for the test branch, e.g. /some/default/path/r<my_branch_name>")
if len(sys.argv) > 1: 
    branch_path = sys.argv[1]
else :
    sys.exit('EXITING - This script needs the full path for your test branch at command line, e.g. monc_control_testsuite_analysis.py /some/default/path/r<my_branch_name>')

print 'branch to test = ', branch_path

########## Main part of script to perform bit comparison test on MONC ###################################
if optimisation_level == 1 :
    #casim_kgo_datadir = 
    sys.exit('EXITING - optimisation level 1 data not available')
else :
    casim_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/'+compiler+'/opt3/casim/'+kgo_dir+'/'

new_casim_datadir = branch_path
    
with open('bit_compare_results.txt', "w") as f:
    f.write('****************Bit comparison results from*************************'+"\n")
    f.write('Known good MONC output = '+casim_kgo_datadir+"\n")
    f.write('New MONC output = '+new_casim_datadir+"\n")

    
for case in test_case :
    with open('bit_compare_results.txt', "a") as f:
        f.write('********************************************************************************************'+"\n")
        f.write('******************************************'+case+'******************************************'+"\n")
    
    test_list = list_setup.test_names(case)
    figure_title_list = list_setup.figure_names(case)
    
    casim_dict = dictionary_setup.casim(case, casim_kgo_datadir, new_casim_datadir)

    # do bit comparison check on MONC and maybe LEM
    for test_idx, test_run in enumerate(test_list) :
        print test_run
        with open('bit_compare_results.txt', "a") as f:
            f.write('=========================================================================================='+"\n")
        bit_compare_success = bit_compare.monc_calc(case, test_run, casim_dict)
        print bit_compare_success
        
        if not bit_compare_success[0] :  
            print 'plot bit comparison from monc to assess issue'
            bit_compare.monc_plot('bit_compare', case, test_run, figure_title_list[test_idx], casim_dict) 
            
            if bit_compare_success[1]:
                bit_compare.monc_plot('kgo_compare_1pe', case, test_run, figure_title_list[test_idx], casim_dict) 
                bit_compare.monc_plot('kgo_compare_32pe', case, test_run, figure_title_list[test_idx], casim_dict)
                


