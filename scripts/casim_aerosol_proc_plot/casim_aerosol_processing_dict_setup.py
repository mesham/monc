#!/usr/bin/env python2.7

import sys

def monc_casim(case, kgo_datadir, new_datadir) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus':
        monc_dict = {'files':['diagnostics_7260.nc', 'diagnostics_7260.nc', 'diagnostics_7260.nc',
                              'diagnostics_7260.nc'],
                     't_labels':['30 mins', '60 mins', '90 mins', '120 min'],
                     'time_color':[ 'blue', 'red', 'y', 'c' ],
                     'no_procs':['36_','72_']} 
    elif case == 'shallow_convection':
        monc_dict = {'files':['diagnostics_21600.nc', 'diagnostics_21600.nc', 'diagnostics_21600.nc',
                              'diagnostics_21600.nc'],
                     't_labels':['210 mins', '240 mins', '270 mins', '300 mins', '330 mins', 
                                 '360 mins'],
                     'time_color':[ 'm', 'g', 'k',  'blue', 'red', 'y' ],
                     'no_procs':['36_','72_']}
    elif case == 'rce':
        monc_dict = {'files':['diagnostics_86400.nc', 'diagnostics_86400.nc', 'diagnostics_86400.nc'],
                     't_labels':['22 hours', '23 hours', '24 hours'],
                     'time_color':[ 'm', 'g', 'k'],
                     'no_procs':['36_','72_']}   
    else :
        sys.exit('testcase is not defined. Please check testcase name, EXITING')
    return monc_dict

