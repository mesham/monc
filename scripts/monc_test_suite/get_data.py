#!/usr/bin/env python2.7

import netCDF4
import numpy as np
import sys

def profile_diagnostics(model_fn, file_numbers, model_data) :
    # first test for whether the model is the LEM
    lem_model=model_fn.count('lem')
    sc_q_incl=model_fn.count('Sc')
    cu_q_incl=model_fn.count('Cu')
    dry_bl=model_fn.count('Dry')
    if not lem_model :
        f_string = model_fn+file_numbers[0]+'.nc'
        print f_string
        f_open=netCDF4.Dataset(f_string,'r') 
        zn = f_open.variables['zn']
        model_data['zn'].append(zn[:])
        ww_bar = f_open.variables['ww_mean']
        uu_bar = f_open.variables['uu_mean']
        vv_bar = f_open.variables['vv_mean']
        model_data['ww_mean'].append(ww_bar[:,:])
        model_data['uu_mean'].append(uu_bar[:,:])
        model_data['vv_mean'].append(vv_bar[:,:])
        u_bar = f_open.variables['u_wind_mean']
        v_bar = f_open.variables['v_wind_mean']
        w_bar = f_open.variables['w_wind_mean']
        model_data['u_mean'].append(u_bar[:,:])
        model_data['v_mean'].append(v_bar[:,:])
        model_data['w_mean'].append(w_bar[:,:])
        if not dry_bl :
            theta = f_open.variables['theta_mean']
            model_data['theta_mean'].append(theta[:,:])
            if sc_q_incl or cu_q_incl :
                ql = f_open.variables['liquid_mmr_mean']
                model_data['liquid_mmr_mean'].append(ql[:,:])
                qv = f_open.variables['vapour_mmr_mean']
                model_data['vapour_mmr_mean'].append(qv[:,:])
        f_open.close()
    else :
        for idx, n in enumerate(file_numbers):
            f_string = model_fn+n+'.nc'
            print f_string
            f_open=netCDF4.Dataset(f_string,'r')
            if (idx == 0) :
                z = f_open.variables['Z']
                zn = f_open.variables['ZN']
                model_data['z'].append(z[:])
                model_data['z'].append(zn[:])

            ww_bar = f_open.variables['WW']
            uu_bar = f_open.variables['UU']
            vv_bar = f_open.variables['VV']            
            model_data['ww_mean'].append(ww_bar[:])
            model_data['uu_mean'].append(uu_bar[:])
            model_data['vv_mean'].append(vv_bar[:])
            u_bar = f_open.variables['UBRAV']
            v_bar = f_open.variables['VBRAV']
            w_bar = f_open.variables['WBRAV']
            model_data['u_mean'].append(u_bar[:])
            model_data['v_mean'].append(v_bar[:])
            model_data['w_mean'].append(w_bar[:])
            if not dry_bl :
                theta = f_open.variables['THBAV']
                model_data['theta_mean'].append(theta[:])
                if q_incl :
                    ql = f_open.variables['QBAV02']
                    model_data['liquid_mmr_mean'].append(ql[:])
                    qv = f_open.variables['QBAV01']
                    model_data['vapour_mmr_mean'].append(qv[:])
            f_open.close()

