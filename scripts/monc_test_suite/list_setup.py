#!/usr/bin/env python2.7

import sys

def test_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'bubble':
        test_list = ['WarmPw_', 'WarmPwNoSmag_', 
                     'WarmTvd_', 'WarmTvdNoSmag_',
                     'ColdPw_', 'ColdPwNoSmag_',   
                     'ColdTvd_', 'ColdTvdNoSmag_',
                     'WarmPwNoSmagGal_', 'WarmTvdNoSmagGal_', 
                     'ColdPwNoSmagGal_', 'ColdTvdNoSmagGal_',
                     'ColdNoSmagGalAdv_', 'ColdNoSmagGalMomAdv_',
                     'WarmNoSmagGalAdv_', 'WarmNoSmagGalMomAdv_' ] 

    elif case == 'stratus':
        test_list = ['ScFull_' , 'ScNoGalSrfRadSubDamp_', 'ScNoSrfSubDamp_',
                     'ScNoDamp_', 'ScNoSrfRadSubDamp_', 'ScNoSubDamp_', 'ScNoRadGal_'] 
        
    elif case == 'drybl':
        test_list = ['DryBlFull_', 
                     'DryBlNoSmagGal_', 
                     'DryBlNoGal_']

    elif case == 'shallow_convection' :
        test_list = ['CuFull_', 'CuNoDamp_', 'CuNoGalSrfForceSubDamp_', 
                     'CuNoSrfForceSubDamp_', 'CuNoSrfSubDamp_', 'CuNoSubDamp_']
        #test_list = ['CuFull_', 'CuNoDamp_']
    else :
        sys.exit('testcase is not defined in test_names. Please check testcase name, EXITING')
    return test_list

def figure_names(case) :
    # call module that sets up the lists and dictionaries for the case
    if case == 'bubble':
        fig_title_list = [ 'warm bubble - P-W mom advection, Galilean Trans & Smag mixing on',
                           'warm bubble - P-W mom advection, Galilean Trans', 
                           'warm bubble - TVD mom advection, Galilean Trans & Smag mixing on',
                           'warm bubble - TVD mom advection, Galilean Trans', 
                           'cold bubble - P-W mom advection, Galilean Trans & Smag mixing on',
                           'cold bubble - P-W mom advection, Galilean Trans',
                           'cold bubble - TVD mom advection, Galilean Trans & Smag mixing on',
                           'cold bubble - TVD mom advection, Galilean Trans',
                           'warm bubble - pressure, P-W momentum, TVD scalars',
                           'warm bubble - pressure, TVD momentum, TVD scalars',
                           'cold bubble - pressure, P-W momentum, TVD scalars',
                           'cold bubble - pressure, TVD momentum, TVD scalars',
                           'cold bubble - pressure, NO advection',
                           'cold bubble - pressure, TVD scalar advection, NO momentum advection',
                           'warm bubble - pressure, NO advection',
                           'warm bubble - pressure, TVD scalar advection, NO momentum advection',
                           'warm bubble - pressure, NO scalar advection, PW momentum advection',
                           'cold bubble - pressure, NO scalar advection, PW momentum advection',
                           'cold bubble - pressure, PW scalar and momentum advection',
                           'warm bubble - pressure, PW scalar and momentum advection']
    elif case == 'stratus':
        fig_title_list = [ 'Stratus - Full',
                           'Stratus - No Galilean, surface flux, radiation param, subsidence, damping',
                           'Stratus - No surface flux, subsidence, damping', 
                           'Stratus - No damping',
                           'Stratus - No surface flux, radiation param, subsidence, damping',
                           'Stratus - No subsidence, damping', 
                           'Stratus - No radiation, galilean transformation']
            
    elif case == 'drybl':
        fig_title_list = [ 'Dry Boundary layer - Full',
                           'Dry Boundary layer - No Smag or Galilean trans',
                           'Dry Boundary layer - No Galilean']
    elif case == 'shallow_convection':
        fig_title_list = [ 'Shallow convection - Full',
                           'Shallow convection - No damping',
                           'Shallow convection - No Galilean, surface flux, forcing, subsidence, damping', 
                           'Shallow convection - No surface flux, forcing, subsidence, damping',
                           'Shallow convection - No surface flux, subsidence, damping',
                           'Shallow convection - No subsidence, damping']
    else :
        sys.exit('testcase is not defined in figure_names. Please check testcase name, EXITING')
    return fig_title_list
