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
# Adhill - 310717

import subprocess as sub
import os
import sys
import getopt
from distutils.util import strtobool

def replaceInFile(file,x,y,fileout=None):
    """ replace occurences of x with y in a file """
    if fileout is None: fileout=file
    fd=open(file,'r')
    text=fd.read()
    fd.close()
    fd=open(fileout, 'w')
    text=text.replace(x,y)
    fd.write(text)

# default script settings
print 'defaults settings'
print 'component testing with cray compiler on MO system'
suite = 'component'
build = True
compiler = ['cray']
system = 'MO-XC40'

print suite, system, build, compiler
if len(sys.argv) > 1 :
    options, args = getopt.getopt(sys.argv[1:],"ht:c:s:b:",["test_suite=", "build=", "compiler=", "system="])

    for opt, arg in options:
        if opt == '-h' :
            print 'Usage of monc_test_suite.py'
            print 'monc_test_suite.py -t <standard or component> -b <True or False> -c <cray or gnu> -s <MO-XC40 or ARCHER>'
            print '-t or --test_suite is the suite type, i.e. component or standard'
            print '-b or --build is whether to build MONC, True or False'
            print '-c or --compiler is the compiler, i.e. cray or gnu'
            print '-s or --system  is the system to run the test suite on, MO-XC40 or ARCHER'
            print 'defaults settings ='
            print 'component, with build = True, cray compiler on MO system'
            sys.exit()
        if opt in ("-t", "--test_suite"):
            suite = arg
        if opt in ("-s", "--system"):
            system = arg   
        if opt in ("-b", "--build"):
            build = strtobool(arg)
        if opt in ('-c', '--compiler'):
            compiler = [arg]

qsub_script_template = 'test_suite/monc_main_test_suite/submonc_template.pbs'

# default io config directory, use the standard monc io configs
ioserver_cfg_dir = 'io/io_cfg_files/'

if suite == 'component' :
    procs = [36, 2] 
    multifile_dgs = False
    testcases = [ 'bubble' , 'dryBl', 'stratus', 'shallow_convection' ]
    #testcases = [ 'shallow_convection' ]
    

    test_cfg={ 'bubble':[ 'ColdNoSmagGalAdv', 'ColdNoSmagGalMomAdv',
                          'ColdPw','ColdPwNoSmag','ColdPwNoSmagGal', 
                          'ColdTvd','ColdTvdNoSmag','ColdTvdNoSmagGal',
                          'WarmNoSmagGalAdv', 'WarmNoSmagGalMomAdv',
                          'WarmPw', 'WarmPwNoSmag', 'WarmPwNoSmagGal',  
                          'WarmTvd', 'WarmTvdNoSmag', 'WarmTvdNoSmagGal'],
               'dryBl':['DryBlFull', 'DryBlNoSmagGal', 'DryBlNoGal'],
               'stratus':['ScFull', 'ScNoDamp', 'ScNoSrfRadSubDamp', 'ScNoRadGal', 
                          'ScNoGalSrfRadSubDamp', 'ScNoSrfSubDamp', 'ScNoSubDamp', 
                          'ScNoRad', 'ScFixFluxNoSubDamp'], 
               'shallow_convection':[ 'CuFull', 'CuNoDamp' , 'CuNoGalSrfForceSubDamp', 
                                      'CuNoSrfForceSubDamp', 'CuNoSrfSubDamp', 'CuNoSubDamp'] 
           }
elif suite == 'standard':
    procs = [2, 36]
    multifile_dgs = False
    testcases = ['tank_experiments', 'drybl', 'stratus', 'shallow_convection']
    test_cfg={'tank_experiments':['cold_bubble', 'warm_bubble', 'straka', 'up_down_bubble',
                                  'bubble_simple'],
              'drybl':['drybl'],
              'stratus':['fire_sc', 'fire_sc_casim'],
              'shallow_convection':['bomex', 'bomex_casim']
          }                                       

for comp in compiler :

    compiler_dir = 'test_suite/monc_main_test_suite/'+comp
        
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
        if suite == 'component' :
            base_cfg_dir = 'test_suite/monc_main_test_suite/'+tc+'/'
        else :
            base_cfg_dir = 'testcases/'+tc+'/'
      
        ioserver_cfg_file=ioserver_cfg_dir+'data_write_1file.xml' 

        for cfg in test_cfg[tc]:
            
            for pe in procs : 
                if suite == 'component' :
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
                    if multifile_dgs :
                        
                        # with multifile_dgs, netcdf output goes to a 3 files, one for profiles, one for scalars and
                        # one for 3-D output
                        replaceInFile(ioserver_cfg_file, 'profile_ts.nc', output_data+'_profile_ts.nc', ioserver_cfg_write_file)
                        replaceInFile(ioserver_cfg_write_file, 'scalar_ts.nc', output_data+'_scalar_ts.nc', ioserver_cfg_write_file)
                        replaceInFile(ioserver_cfg_write_file, '3dfields_ts.nc', output_data+'_3dfields_ts.nc', ioserver_cfg_write_file)
                    else :
                        # all diagnostics written to one file
                        replaceInFile(ioserver_cfg_file, 'diagnostic_files/diagnostics_ts.nc', output_data+'_diagnostics.nc', ioserver_cfg_write_file)
                else :

                    # 1. set up config file and write to the temporary location
                    config_file = base_cfg_dir+cfg+'.mcf'
                    new_cfg_file = tmp_config_dir+cfg+'_'+str(pe)+'.mcf'
                    # 2. replace the dump file name with a unique name for the config
                    chkpoint_fname = cfg+'_dump.nc'
                    new_chkpoint_fname = checkpoint_output_dir+cfg+'_'+str(pe)+'.nc'
                    replaceInFile(config_file, chkpoint_fname, new_chkpoint_fname, new_cfg_file)

                    # 3. set up the netcdf output file name and add to the ioserver_cfg file
                    output_data = nc_output_dir+cfg+'_'+str(pe)
                    ioserver_cfg_write_file = tmp_ioserver_cfg_dir+cfg+'_'+str(pe)+'.xml'

                    # modify config file with the new specific ioserver cfg file path
                    # Standard config should include an existing ioserver_config that points 
                    # to io/io_cfg_files/
                    replaceInFile(new_cfg_file, ioserver_cfg_file, ioserver_cfg_write_file, new_cfg_file)
                    # all diagnostics written to one file in the standard case
                    replaceInFile(ioserver_cfg_file, 'diagnostic_files/diagnostics_ts.nc', output_data+'_diagnostics.nc', ioserver_cfg_write_file)
               
            # set up qsub submission file
                new_scriptfile = tmp_qsub_dir+'submonc_'+cfg+'_'+str(pe) 
                stdout_file = out_dir+'out_'+cfg+'_'+str(pe)
                replaceInFile(qsub_script_template, 'np', str(pe), new_scriptfile)
                replaceInFile(new_scriptfile, 'add_scriptname', new_scriptfile, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_testcase_mcf', new_cfg_file, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_stdout_dirname', out_dir, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_checkpoint_dirname', checkpoint_output_dir, new_scriptfile)
                replaceInFile(new_scriptfile, 'add_testsuite_jobname', cfg+'_'+str(pe)+'_' , new_scriptfile)
                replaceInFile(new_scriptfile, 'out2', stdout_file , new_scriptfile)
                                
                #if suite == 'standard':
                #    replaceInFile(new_scriptfile, "output_file", new_chkpoint_fname, new_scriptfile)
            
                print 'submitting '+new_scriptfile+' to xcs-c'
                p=sub.call(['qsub', new_scriptfile ])
                




        
        
        
        

