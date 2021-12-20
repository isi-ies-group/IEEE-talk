from math import sqrt
import matplotlib.pyplot as plt

import pvlib

import cpvlib

A_ref = 10 # la eficiencia está ajustada para el FF con A_ref=10

# Canadian Solar CS1U-410MS - PVSyst
eff_pv = 20
A_pv = 2.061 * (20.5/eff_pv)  # m2

corr_pv = A_ref / A_pv
A_pv *= corr_pv

pv_mod_params = {
    "alpha_sc": 4.8e-3,  # coef. temp. Isc
    "gamma_ref": 0.967,  # "Datos básicos"
    "mu_gamma": -0.00042,  # "Parámetros modelo" [1/K]
    "I_L_ref": 9.7 *sqrt(corr_pv), # Isc
    "I_o_ref": 0.03e-9,  # "Datos básicos"
    "R_sh_ref": 600,  # R paral ref "Parámetros modelo"
    "R_sh_0": 2500,  # R paral G=0 W/m2 "Parámetros modelo"
    "R_s": 0.291,  # R serie "Parámetros modelo"
    "cells_in_series": 81 *sqrt(corr_pv),
}

eff_cpv = 30
# Soitec CX-M500
A_cpv = 7.386 * (34.87/eff_cpv)  # m2

corr_cpv = A_ref / A_cpv
A_cpv *= corr_cpv
cpv_mod_params = {
    "alpha_sc": 0.00,
    "gamma_ref": 3.664,
    "mu_gamma": 0.003,
    "I_L_ref": 3.861 *1.274 *sqrt(corr_cpv),
    "I_o_ref": 0.005e-9,
    "R_sh_ref": 3461,
    "R_sh_0": 25000,
    "R_s": 0.61,
    "EgRef": 3.91,
    "cells_in_series": 240 *sqrt(corr_cpv),
    "irrad_ref":943,
    "temp_ref":64
}

UF_parameters_cpv = {
    "IscDNI_top": 1,
    "am_thld": 1,#1.7,
    "am_uf_m_low": 0,#0.1,
    "am_uf_m_high": 0,#-0.1,
    # "ta_thld": 25,
    # "ta_uf_m_low": 0.005,
    # "ta_uf_m_high": 0,
    # "weight_am": 0.55,
    # "weight_temp": 0.45,
}

cpv_mod_params.update(UF_parameters_cpv)

parameters_tracker = {
    "axis_tilt":10,
    "axis_azimuth":180,
    "max_angle":90, 
    "backtrack":True,
    "gcr":2/7
    }

def genera_pot_pv(location, solpos, data, tilt, diffuse_model, in_singleaxis_tracker=False):
    
    temp_mod_params = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['pvsyst']['freestanding']

    pv_sys = pvlib.pvsystem.PVSystem(
        surface_tilt=tilt,
        surface_azimuth=180,
        module_parameters=pv_mod_params,
        temperature_model_parameters=temp_mod_params,
        modules_per_string=1,
    )
    
    if in_singleaxis_tracker:
        
        tracking_info = pvlib.tracking.singleaxis(
                apparent_zenith=solpos['zenith'], 
                apparent_azimuth=solpos['azimuth'], **parameters_tracker)

        surface_tilt = tracking_info.surface_tilt
        surface_azimuth = tracking_info.surface_azimuth
        
        aoi = pvlib.tracking.singleaxis(
                apparent_zenith=solpos['zenith'], 
                apparent_azimuth=solpos['azimuth'], **parameters_tracker).aoi
    else:
        surface_tilt = pv_sys.surface_tilt
        surface_azimuth = pv_sys.surface_azimuth
        
        aoi = pv_sys.get_aoi(
            solar_zenith=solpos['zenith'],
            solar_azimuth=solpos['azimuth'],
        )

    pv_irr = pvlib.irradiance.get_total_irradiance(
        surface_tilt=surface_tilt,
        surface_azimuth=surface_azimuth,
        solar_zenith=solpos['zenith'],
        solar_azimuth=solpos['azimuth'],
        ghi=data['ghi'],
        dhi=data['dhi'],
        dni=data['dni'],
        dni_extra=pvlib.irradiance.get_extra_radiation(data.index),
        model=diffuse_model
        )
    
    # Modelo pérdidas angulares
    irradiance = pv_irr['poa_global']
    effective_irradiance = irradiance * pvlib.iam.martin_ruiz(aoi, a_r=0.16) # default value
    
    cell_temp = pv_sys.pvsyst_celltemp(
        poa_global=effective_irradiance,
        temp_air=data['temp_air'],
        wind_speed=data['wind_speed']
    )
    
    diode_parameters = pv_sys.calcparams_pvsyst(
    effective_irradiance=effective_irradiance,
    temp_cell=cell_temp,
    )
    
    power = pv_sys.singlediode(*diode_parameters)
    p_mp = power['p_mp']
    
    diode_parameters_25 = pv_sys.calcparams_pvsyst(
    effective_irradiance=effective_irradiance,
    temp_cell=25,
    )
    p_mp_25 = pv_sys.singlediode(*diode_parameters_25)['p_mp']
    
    # calcula Pmp STC
    Pdc_stc = pvlib.pvsystem.singlediode(*pvlib.pvsystem.PVSystem(
        module_parameters=pv_mod_params
        ).calcparams_pvsyst(
        effective_irradiance=1000,
        temp_cell=25))['p_mp']
    
    eff_a = Pdc_stc / (1000 * A_pv)
    
    return irradiance, p_mp, p_mp_25, Pdc_stc, eff_a

def genera_pot_cpv(location, solpos, data, tilt, eff_opt_cpv):
    
    cpv_mod_params_copy = cpv_mod_params.copy()
    cpv_mod_params_copy["I_L_ref"] *= eff_opt_cpv
    
    temp_mod_params = {"eta_m": 0.32, "u_c":29.0, "u_v":0.6} # Gerstmaier, Tobias et al «Validation of the PVSyst Performance Model for the Concentrix CPV Technology»

    cpv_sys = cpvlib.CPVSystem(
        module_parameters=cpv_mod_params_copy,
        temperature_model_parameters=temp_mod_params,
        modules_per_string=1,
    )
    
    #print(cpv_sys.module_parameters)
    irradiance = data['dni']
    
    cell_temp = cpv_sys.pvsyst_celltemp(
    poa_global=irradiance,
    temp_air=data['temp_air'],
    wind_speed=data['wind_speed']
    )

    diode_parameters = cpv_sys.calcparams_pvsyst(
        effective_irradiance=irradiance,
        temp_cell=cell_temp,
    )
    
    power = cpv_sys.singlediode(*diode_parameters)
    p_mp = power['p_mp']

    diode_parameters_25 = cpv_sys.calcparams_pvsyst(
    effective_irradiance=irradiance,
    temp_cell=25,
    )
    p_mp_25 = cpv_sys.singlediode(*diode_parameters_25)['p_mp']
    
    data['am'] = location.get_airmass(data.index).airmass_absolute

    uf_cpv = cpv_sys.get_am_util_factor(data['am'])
    
    p_mp_uf = p_mp * uf_cpv
    p_mp_uf_25 = p_mp_25 * uf_cpv
    
    # calcula Pmp STC
    Pdc_stc = pvlib.pvsystem.singlediode(*cpvlib.CPVSystem(
        module_parameters=cpv_mod_params_copy
        ).calcparams_pvsyst(
        effective_irradiance=1000,
        temp_cell=25))['p_mp']
    
    eff_a = Pdc_stc / (1000 * A_cpv)
    
    return irradiance, p_mp_uf, p_mp_uf_25, Pdc_stc, eff_a

def genera_pot_static_cpv(location, solpos, data, tilt, eff_opt_cpv, in_singleaxis_tracker=False):

    cpv_mod_params_copy = cpv_mod_params.copy()
    cpv_mod_params_copy["I_L_ref"] *= eff_opt_cpv

    cpv_temp_mod_params = {"eta_m": 0.32, "u_c":29.0, "u_v":0.6}

    static_cpv_sys = cpvlib.StaticCPVSystem(
            surface_tilt=tilt,
            surface_azimuth=180,
            module_parameters=cpv_mod_params_copy,
            temperature_model_parameters=cpv_temp_mod_params,
            modules_per_string=1,
            in_singleaxis_tracker=in_singleaxis_tracker,
            parameters_tracker=parameters_tracker
    )

    irradiance = static_cpv_sys.get_irradiance(
        solpos['zenith'], solpos['azimuth'], data['dni'])

    aoi = static_cpv_sys.get_aoi(
        solar_zenith=solpos['zenith'],
        solar_azimuth=solpos['azimuth'],
    )

    # theta_ref = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    # iam_ref = [1.000, 1.007, 0.998, 0.991, 0.971, 0.966, 0.938, 0.894, 0.830, 0.790, 0.740, 0.649, 0.387]
    
    theta_ref = [0, 30, 55, 75]
    iam_ref = [1.0, 1.0, 0.9, 0]
    
    cpv_effective_irradiance = irradiance * pvlib.iam.interp(aoi, theta_ref, iam_ref, method='linear')

    #pd.Series(iam_ref, theta_ref).plot()

    cell_temp = static_cpv_sys.pvsyst_celltemp(
        poa_global=cpv_effective_irradiance,
        temp_air=data['temp_air'],
        wind_speed=data['wind_speed']
    )

    diode_parameters = static_cpv_sys.calcparams_pvsyst(
        effective_irradiance=cpv_effective_irradiance,
        temp_cell=cell_temp,
    )

    power_no_uf = static_cpv_sys.singlediode(*diode_parameters)

    diode_parameters_25 = static_cpv_sys.calcparams_pvsyst(
    effective_irradiance=irradiance,
    temp_cell=25,
    )
    power_no_uf_25 = static_cpv_sys.singlediode(*diode_parameters_25)

    data['am'] = location.get_airmass(data.index).airmass_absolute

    uf_cpv = static_cpv_sys.get_am_util_factor(data['am'])

    p_mp_uf = power_no_uf['p_mp'] * uf_cpv
    p_mp_uf_25 = power_no_uf_25['p_mp'] * uf_cpv
    
    # calcula Pmp STC
    Pdc_stc = pvlib.pvsystem.singlediode(*cpvlib.StaticCPVSystem(
        module_parameters=cpv_mod_params_copy
        ).calcparams_pvsyst(
        effective_irradiance=1000,
        temp_cell=25))['p_mp']
    
    eff_a = Pdc_stc / (1000 * A_cpv)
    
    return irradiance, p_mp_uf, p_mp_uf_25, Pdc_stc, eff_a

def genera_pot_flatplate(location, solpos, data, diffuse_model, tilt, 
                          aoi_limit, eff_opt_pv, cpv_irradiance_spillage, 
                          type_irr_input, in_singleaxis_tracker=False):

    pv_mod_params_copy = pv_mod_params.copy()
    pv_mod_params_copy["I_L_ref"] *= eff_opt_pv
    
    pv_temp_mod_params = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['pvsyst']['freestanding']

    static_flatplate_sys = cpvlib.StaticFlatPlateSystem(
        surface_tilt=tilt,
        surface_azimuth=180,
        module_parameters=pv_mod_params_copy,
        temperature_model_parameters=pv_temp_mod_params,
        modules_per_string=1,
        in_singleaxis_tracker=in_singleaxis_tracker,
        parameters_tracker=parameters_tracker
    )
    
    if in_singleaxis_tracker:
        tracking_info = pvlib.tracking.singleaxis(
        apparent_zenith=solpos['zenith'], 
        apparent_azimuth=solpos['azimuth'], **static_flatplate_sys.parameters_tracker)

        surface_tilt = tracking_info.surface_tilt
        surface_azimuth = tracking_info.surface_azimuth
        
    else:
        surface_tilt = static_flatplate_sys.surface_tilt
        surface_azimuth = static_flatplate_sys.surface_azimuth
        
    # el objeto static_flatplate_sys ya contiene el flag de tracker y sabe devolver bien
    aoi = static_flatplate_sys.get_aoi(
        solar_zenith=solpos['zenith'],
        solar_azimuth=solpos['azimuth']
        )
            
    pv_irradiance = pvlib.irradiance.get_total_irradiance(
        surface_tilt=surface_tilt, surface_azimuth=surface_azimuth,
        solar_zenith=solpos['zenith'], solar_azimuth=solpos['azimuth'],
        dni=data['dni'], ghi=data['ghi'], dhi=data['dhi'],
        dni_extra=pvlib.irradiance.get_extra_radiation(data.index), model=diffuse_model
    )['poa_diffuse']
    
    if type_irr_input == 'diffuse':
        pv_effective_irradiance = (
            static_flatplate_sys.get_irradiance(
            solar_zenith=solpos['zenith'], solar_azimuth=solpos['azimuth'],
            aoi=aoi, aoi_limit=aoi_limit,
            dni=data['dni'], ghi=data['ghi'], dhi=data['dhi'], model=diffuse_model)
        ) * pvlib.iam.martin_ruiz(aoi, a_r=0.16) + cpv_irradiance_spillage
    
    elif type_irr_input == 'dni':            
        pv_effective_irradiance = (
            pvlib.irradiance.beam_component(
                surface_tilt=surface_tilt,
                surface_azimuth=surface_azimuth,
                solar_zenith=solpos['zenith'],
                solar_azimuth=solpos['azimuth'],
                dni=data['dni'])
            ) * pvlib.iam.martin_ruiz(aoi, a_r=0.16) + cpv_irradiance_spillage
    else:
        raise SystemError
    
    cell_temp = static_flatplate_sys.pvsyst_celltemp(
        poa_flatplate_static=pv_effective_irradiance,
        temp_air=data['temp_air'],
        wind_speed=data['wind_speed']
    )
    
    diode_parameters = static_flatplate_sys.calcparams_pvsyst(
        effective_irradiance=pv_effective_irradiance,
        temp_cell=cell_temp,
    )
    
    # plt.figure();cell_temp.hist(bins=50)
    # plt.figure();data['temp_air'].hist(bins=50)
    
    power = static_flatplate_sys.singlediode(*diode_parameters)

    diode_parameters_25 = static_flatplate_sys.calcparams_pvsyst(
        effective_irradiance=pv_effective_irradiance,
        temp_cell=25,
    )
    power_25 = static_flatplate_sys.singlediode(*diode_parameters_25)
    
    p_mp = power['p_mp']
    p_mp_25 = power_25['p_mp']
    
    # calcula Pmp STC
    Pdc_stc = pvlib.pvsystem.singlediode(*cpvlib.StaticFlatPlateSystem(
        module_parameters=pv_mod_params_copy
        ).calcparams_pvsyst(
        effective_irradiance=1000,
        temp_cell=25))['p_mp']
    
    eff_a = Pdc_stc / (1000 * A_pv)
    
    return pv_irradiance, p_mp, p_mp_25, Pdc_stc, eff_a