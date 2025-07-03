#!/usr/bin/env python2.7

import subprocess as sub
import sys

test_suite_script_dir= os.environ.get("HARNESS_SCRIPT_DIRECTORY",
                                      "please provide the full path for the test suite script directory, e.g. /some/default")
if len(sys.argv) > 1: test_suite_script_dir = sys.argv[1]

print test_suite_script_dir

testcases = ['bubble','stratus', 'shallow_convection','dryBl']
tc_out_time = ['900', '36000', '40000', '64802']

simulation_configs={'bubble':['ColdTvdNoSmag', 'ColdTvdNoSmag1d', 'WarmTvdNoSmag', 'WarmTvdNoSmag1d',
                                  'ColdPwNoSmag', 'ColdPwNoSmag1d','WarmPwNoSmag', 'WarmPwNoSmag1d',
                                  'ColdTvdNoLbc', 'ColdTvdNoLbc1d', 'WarmTvdNoLbc', 'WarmTvdNoLbc1d',
                                  'ColdPwNoLbc', 'ColdPwNoLbc1d','WarmPwNoLbc', 'WarmPwNoLbc1d'],
                    'stratus':['ScFull', 'ScFull1d','ScFullnoDmp', 'ScFullnoDmp1d'], 
                    'shallow_convection':['CuFull1d', 'CuFullnoDmp1d','CuFull', 'CuFullnoDmp'],
                    'dryBl':['DryBlFull', 'DryBlFull1d']
                         }
procs = ['1', '32']

# the control data is on MONSOON
control_output_dir = '/projects/monc/adhill/monc_test_harness_data/'
multi_pe_output_dir = 'test_suite/ncfiles/'

for tc_idx, tc in enumerate(testcases) :
    for cfg in simulation_configs[tc] :
        file_1pe = control_output_dir+cfg+'_'+procs[0]+'_'+tc_out_time[tc_idx]+'.nc'
        file_32pe = multi_pe_output_dir+cfg+'_'+procs[1]+'_'+tc_out_time[tc_idx]+'.nc'
        
        print file_1pe, file_32pe

        p=sub.call([test_suite_script_dir'/run_compare.sh', file_1pe, file_32pe, cfg])
        
