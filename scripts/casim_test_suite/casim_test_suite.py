#!/usr/bin/env python2.7

# vn1.0 test-suite 
# This has been developed to automate the component testing of MONC. This is done
# by running a suite of tests that incrementally increase the complexity of 
# of the simulation by switching on components within the case. The component 
# testing is invoked by setting suite = 'component' and this will run the 
# configurations in the test_suite directory
# 
# As well as the component testing, this script can also run the standard cases
# so that a user can check that code change has not broken any of the standard cases.
# This is invoked by setting suite = 'standard'. This will use the configurations in 
# the testcases directory, but will copy the configurations and io configurations to
# the test_suite directory, so that the testcases remain stand alone.
#
# This script is an example of funtionality required for the development of a 
# Rose suite to automate the testing of the MONC. The script below is specific to MONSOON
# This script does not adhere to coding standards for MONC but is added to the trunk so 
# that the testing is traceable.
#
# Adhill - 120617

import shutil
import subprocess as sub
import os
import sys

def replaceInFile(file,x,y,fileout=None):
    """ replace occurences of x with y in a file """
    if fileout is None: fileout=file
    fd=open(file,'r')
    text=fd.read()
    fd.close()
    fd=open(fileout, 'w')
    text=text.replace(x,y)
    fd.write(text)

build = True
compiler = ['gnu']
system = 'MO-XC40'
qsub_script_template = 'test_suite/monc_casim_test_suite/submonc_template.pbs'

# default io config directory is the same irrespective of whether checkpoint is used or not
ioserver_cfg_dir = 'io/io_cfg_files/'

procs = [72,36] 
multifile_dgs = False
testcases = [ 'stratus', 'shallow_convection', 'rce' ]

test_cfg={ 
    'stratus':[ 'ScNoSubDamp_2M_Ndfix', 'ScNoDamp_2M_Ndfix', 'ScFull_2M_Ndfix', 
                'ScNoSubDamp_Socrates_2M_Ndfix', 'ScNoDamp_Socrates_2M_Ndfix', 'ScFull_Socrates_2M_Ndfix', 
                'ScNoDamp_2M_Ndfix_diurnal' ],
    'shallow_convection':['CuNoSubDamp_2M_Ndfix','CuNoDamp_2M_Ndfix', 'CuNoDamp_Socrates_2M_Ndfix', 'CuFull_2M_Ndfix'],
    'rce':['RCENoDampSocrates_2M_Ndfix','RCENoDampForce_2M_Ndfix', 'RCENoDamp_2M_Ndfix', 'RCE_2M_Ndfix', 
           'RCESocrates_2M_Ndfix', 'RCENoDampNoUVforce_2M_Ndfix']
    }

for comp in compiler :
    compiler_dir = 'test_suite/monc_casim_test_suite/'+comp
        
    if not os.path.exists(compiler_dir):
        os.makedirs(compiler_dir)
   
    if build :
        p=sub.call(['monc_build_cray.sh', comp, system])

## This script makes the following directories, if they do not exist
##
# this directory contains the configs that are run
# the files are built using the base configs in bubble, stratus
# directories. These will all be the same except the xml file
    tmp_config_dir =compiler_dir+'/tmp_config_dir/'
    tmp_ioserver_cfg_dir = compiler_dir+'/tmp_ioserver_cfg_dir/'
# Submission scripts are saved in this directory
    tmp_qsub_dir = compiler_dir+'/tmp_qsub_dir/'
# Netcdf output is saved in this directory
    nc_output_dir = compiler_dir+'/ncfiles/'
# Checkpoint output is save in this directory
    checkpoint_output_dir = compiler_dir+'/chk_point_dump/'
# the standard output is stored here. This is required for 
# debug and can be used for control of the test harness
    out_dir = compiler_dir+'/monc_stdout/'

# test for directory structure
    if not os.path.exists(tmp_config_dir):
        os.makedirs(tmp_config_dir)

    if not os.path.exists(tmp_ioserver_cfg_dir):
        os.makedirs(tmp_ioserver_cfg_dir)

    if not os.path.exists(tmp_qsub_dir):
        os.makedirs(tmp_qsub_dir)

    if not os.path.exists(nc_output_dir):
        os.makedirs(nc_output_dir)
    
    if not os.path.exists(checkpoint_output_dir):
        os.makedirs(checkpoint_output_dir)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    for tc in testcases :
        # set up the base directory with the default cfg files
        base_cfg_dir = 'test_suite/monc_casim_test_suite/'+tc+'/'
        ioserver_cfg_file=ioserver_cfg_dir+'data_write_1file.xml' 

        for cfg in test_cfg[tc]:
            
            for pe in procs : 
                # 1. set up config file and write to the temporary location
                config_file = base_cfg_dir+cfg
                new_cfg_file = tmp_config_dir+cfg+'_'+str(pe)
                # 2. set up the netcdf output file name and add to the ioserver_cfg file
                output_data = nc_output_dir+cfg+'_'+str(pe)
                ioserver_cfg_write_file = tmp_ioserver_cfg_dir+cfg+'_'+str(pe)+'.xml'
                # modify config file with the new specific ioserver cfg file path
                # following line adds path of the ioserver_config to the temporary config
                # this is only needed for component testing
                replaceInFile(config_file, 'add_ioconfig_file_path', ioserver_cfg_write_file, new_cfg_file)
                # 3. setup the checkpoint file name in the config file
                checkpoint_file=checkpoint_output_dir+cfg+'_'+str(pe)+'.nc'
                replaceInFile(new_cfg_file, 'add_checkpoint_file', checkpoint_file, new_cfg_file)
                #
                # all diagnostics written to one file
                replaceInFile(ioserver_cfg_file, 'diagnostic_files/diagnostics_ts.nc', output_data+'_diagnostics.nc', ioserver_cfg_write_file)
                # set up qsub submission file
                new_scriptfile = tmp_qsub_dir+'submonc_'+cfg+'_'+str(pe) 
                stdout_file = out_dir+'out_'+cfg+'_'+str(pe)
                replaceInFile(qsub_script_template, 'np', str(pe), new_scriptfile)
                if pe > 36 :
                     replaceInFile(new_scriptfile, 'select=1', 'select=2', new_scriptfile)
                replaceInFile(new_scriptfile, 'add_scriptname', new_scriptfile, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_testcase_mcf', new_cfg_file, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_stdout_dirname', out_dir, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_checkpoint_dirname', checkpoint_output_dir, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_testsuite_jobname', cfg+'_'+str(pe)+'_' , new_scriptfile)
                replaceInFile(new_scriptfile, 'out2', stdout_file , new_scriptfile)
                           
                print 'submitting '+new_scriptfile+' to xcs-c'
                p=sub.call(['qsub', new_scriptfile ])
                




        
        
        
        

