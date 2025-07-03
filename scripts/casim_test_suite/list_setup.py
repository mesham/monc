#!/usr/bin/env python2.7

import sys

def test_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus':
        test_list = ['ScFull_2M_Ndfix_', 'ScNoSubDamp_2M_Ndfix_', 'ScNoDamp_2M_Ndfix_',  
                     'ScFull_Socrates_2M_Ndfix_', 'ScNoSubDamp_Socrates_2M_Ndfix_', 'ScNoDamp_Socrates_2M_Ndfix_'] 
    elif case == 'stratus_diurnal':
        test_list = ['ScNoDamp_2M_Ndfix_diurnal_']     
    elif case == 'shallow_convection' :
        test_list = ['CuFull_2M_Ndfix_', 'CuNoDamp_2M_Ndfix_', 'CuNoDamp_Socrates_2M_Ndfix_', 'CuNoSubDamp_2M_Ndfix_']
    elif case == 'rce' :
        test_list = ['RCE_2M_Ndfix_', 'RCENoDamp_2M_Ndfix_', 'RCENoDampForce_2M_Ndfix_',
                     'RCENoDampSocrates_2M_Ndfix_', 'RCENoDampNoUVforce_2M_Ndfix_', 
                     'RCESocrates_2M_Ndfix_']
    else :
        sys.exit('testcase is not defined in test_names. Please check testcase name, EXITING')
    return test_list

def figure_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'stratus':
        fig_title_list = [ 'Stratus - Full, casim, Nd = 50',
                           'Stratus - No subsidence, damping, casim, Nd = 50', 
                           'Stratus - No damping, casim, Nd = 50', 
                           'Stratus - Full, casim, socreate, Nd = 50',
                           'Stratus - No subsidence, damping, casim, socrates, Nd = 50', 
                           'Stratus - No damping, casim, socrates, Nd = 50']
    elif case == 'stratus_diurnal' :
        fig_title_list = ['Diurnal Stratus - No damping, casim, Nd = 50'] 
    elif case == 'shallow_convection':
        fig_title_list = [ 'Shallow convection - Full, casim, Nd = 50',
                           'Shallow convection - No damping, casim, Nd = 50', 
                           'Shallow convection - No damping, casim, socrates, Nd = 50',
                           'Shallow convection - No subsidence, damping, casim, Nd = 50']
    elif case == 'rce' :
        fig_title_list = [ 'RCE - Full, prescribed cooling, casim, Nd = 50',
                           'RCE - No Damp , casim, Nd = 50', 
                           'RCE - No Damp, no force, casim, Nd = 50',
                           'RCE - No Damp, socrates, casim, Nd = 50',
                           'RCE - No Damp, no wind relax, socrates, casim, Nd = 50', 
                           'RCE - Full, socrates, casim, Nd = 50']
    else :
        sys.exit('testcase is not defined in figure_names. Please check testcase name, EXITING')
    return fig_title_list
