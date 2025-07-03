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
import monc_vs_lem
import bit_compare
import dictionary_setup
import list_setup
import sys
import os

print 'monc_control_testsuite called'
       
test_case = [ 'bubble',  'drybl', 'stratus', 'shallow_convection']
optimisation_level = 3
monc_bit_compare = True
monc_lem_compare = False
checkpoint_test = True
compiler = 'cray'
#kgo_dir = 'latest_kgo/r4677_vn0.9_part3_r4849'
if compiler == 'cray' :
    kgo_dir = 'latest_kgo/r6006/'
elif compiler == 'gnu' :
    kgo_dir = 'latest_kgo/r6006/'

branch_path= os.environ.get("BRANCH-PATH",
                            "please provide the full path for the test branch, e.g. /some/default/path/r<my_branch_name>")
if len(sys.argv) > 1: 
    branch_path = sys.argv[1]
else :
    sys.exit('EXITING - This script needs the full path for your test branch at command line, e.g. monc_control_testsuite_analysis.py /some/default/path/r<my_branch_name>')

print 'branch to test = ', branch_path

########## Main part of script to perform bit comparison test on MONC ###################################
if monc_bit_compare :
    if optimisation_level == 1 :
        monc_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/'+compiler+'/opt1/monc/'
    else :
        monc_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/'+compiler+'/opt3/monc/'+kgo_dir+'/'
    if checkpoint_test :
        #new_monc_datadir = branch_path+'test_suite/monc_main_test_suite/'+compiler+'/ncfiles/'
        new_monc_datadir = branch_path
    else :
        new_monc_datadir = branch_path+'test_suite/'+compiler+'/ncfiles/'
    
    with open('bit_compare_results.txt', "w") as f:
        f.write('****************Bit comparison results from*************************'+"\n")
        f.write('Known good MONC output = '+monc_kgo_datadir+"\n")
        f.write('New MONC output = '+new_monc_datadir+"\n")

    
    for case in test_case :
        with open('bit_compare_results.txt', "a") as f:
            f.write('********************************************************************************************'+"\n")
            f.write('******************************************'+case+'******************************************'+"\n")

        test_list = list_setup.test_names(case)
        figure_title_list = list_setup.figure_names(case)

        monc_dict = dictionary_setup.monc(case, monc_kgo_datadir, new_monc_datadir)

        # do bit comparison check on MONC and maybe LEM
        for test_idx, test_run in enumerate(test_list) :
            print test_run
            with open('bit_compare_results.txt', "a") as f:
                f.write('=========================================================================================='+"\n")
            bit_compare_success = bit_compare.monc_calc(case, test_run, monc_dict)
        
            if not bit_compare_success[0] :  

                print 'plot bit comparison from monc to assess issue'
                bit_compare.monc_plot('bit_compare', case, test_run, figure_title_list[test_idx], monc_dict) 

                if bit_compare_success[1]:
                    bit_compare.monc_plot('kgo_compare_1pe', case, test_run, figure_title_list[test_idx], monc_dict) 
                    bit_compare.monc_plot('kgo_compare_32pe', case, test_run, figure_title_list[test_idx], monc_dict)
                

########################################################################################################
#
#
######### Qualitatively compare LEM and MONC output from component testing #############################
if monc_lem_compare :
    if optimisation_level == 1 :
        monc_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/'+compiler+'/opt1/monc/'
    else :
        monc_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/'+compiler+'/opt3/monc/'+kgo_dir+'/'
    new_monc_datadir = branch_path
    
    lem_kgo_datadir = '/projects/monc/fra23/LEM_MONC_comparison/cray/opt1/lem/'
    
    for case in test_case :
    
        test_list = list_setup.test_names(case)
        figure_title_list = list_setup.figure_names(case)

        monc_dict = dictionary_setup.monc(case, monc_kgo_datadir, new_monc_datadir)
        lem_dict = dictionary_setup.lem(case, lem_kgo_datadir)
        
        monc_vs_lem.plotting(case, test_list, figure_title_list, monc_dict, lem_dict)
