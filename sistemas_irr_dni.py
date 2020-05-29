from math import sqrt

import pvlib

import cpvlib

def genera_pot_cpv(location, solpos, data, tilt, eff_opt_cpv, spillage_factor):
    A_ref = 10

    # Soitec CX-M500

    A = 7.386  # m2

    corr = A_ref / A
    A *= corr
    cpv_mod_params = {
        "alpha_sc": 0.00,
        "gamma_ref": 3.664,
        "mu_gamma": 0.003,
        "I_L_ref": 3.861 *1.274*eff_opt_cpv *sqrt(corr),
        "I_o_ref": 0.005e-9,
        "R_sh_ref": 3461,
        "R_sh_0": 25000,
        "R_s": 0.61,
        "EgRef": 3.91,
        "cells_in_series": 240 *sqrt(corr),
        "irrad_ref":943,
        "temp_ref":64
    }

    UF_parameters_cpv = {
        "IscDNI_top": 1,
        "am_thld": 1.7,
        "am_uf_m_low": 0.1,
        "am_uf_m_high": -0.1,
        "ta_thld": 25,
        "ta_uf_m_low": 0.005,
        "ta_uf_m_high": 0,
        "weight_am": 0.55,
        "weight_temp": 0.45,
    }

    cpv_mod_params.update(UF_parameters_cpv)

    cpv_temp_mod_params = {"eta_m": 0.32, "u_c":29.0, "u_v":0.6}

    static_cpv_sys = cpvlib.StaticCPVSystem(
            surface_tilt=tilt,
            surface_azimuth=180,
            module_parameters=cpv_mod_params,
            temperature_model_parameters=cpv_temp_mod_params,
            modules_per_string=1,
    )

    cpv_irradiance = static_cpv_sys.get_irradiance(
        solpos['zenith'], solpos['azimuth'], data['dni'])

    aoi = static_cpv_sys.get_aoi(
        solar_zenith=solpos['zenith'],
        solar_azimuth=solpos['azimuth'],
    )

    cpv_irradiance_spillage = cpv_irradiance * spillage_factor

    theta_ref = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60]
    iam_ref = [1.000, 1.007, 0.998, 0.991, 0.971, 0.966, 0.938, 0.894, 0.830, 0.790, 0.740, 0.649, 0.387]

    cpv_effective_irradiance = cpv_irradiance * pvlib.iam.interp(aoi, theta_ref, iam_ref, method='linear')

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

    cpv_power_no_uf = static_cpv_sys.singlediode(*diode_parameters)

    data['am'] = location.get_airmass(data.index).airmass_absolute

    uf_cpv = static_cpv_sys.get_am_util_factor(data['am'])

    cpv_power_p_mp = cpv_power_no_uf['p_mp'] * uf_cpv
    
    return cpv_power_p_mp, aoi, cpv_irradiance_spillage

def genera_pot_pv(location, solpos, data, diffuse_model, tilt, aoi_limit, eff_opt_pv, cpv_irradiance_spillage):
    A_ref = 10

    # Canadian Solar CS1U-410MS - PVSyst
    A = 2.061  # m2

    corr = A_ref / A
    A *= corr

    pv_mod_params = {
        "alpha_sc": 4.8e-3,  # coef. temp. Isc
        "gamma_ref": 0.967,  # "Datos básicos"
        "mu_gamma": -0.00042,  # "Parámetros modelo" [1/K]
        "I_L_ref": 9.7 *eff_opt_pv *sqrt(corr), # Isc
        "I_o_ref": 0.03e-9,  # "Datos básicos"
        "R_sh_ref": 600,  # R paral ref "Parámetros modelo"
        "R_sh_0": 2500,  # R paral G=0 W/m2 "Parámetros modelo"
        "R_s": 0.291,  # R serie "Parámetros modelo"
        "cells_in_series": 81 *sqrt(corr),
    }

    pv_temp_mod_params = pvlib.temperature.TEMPERATURE_MODEL_PARAMETERS['pvsyst']['freestanding']

    static_flatplate_sys = cpvlib.StaticFlatPlateSystem(
        surface_tilt=tilt,
        surface_azimuth=180,
        module_parameters=pv_mod_params,
        temperature_model_parameters=pv_temp_mod_params,
        modules_per_string=1,
    )
    
    pv_irradiance = pvlib.irradiance.get_total_irradiance(
        static_flatplate_sys.surface_tilt, static_flatplate_sys.surface_azimuth,
        solar_zenith=solpos['zenith'], solar_azimuth=solpos['azimuth'],
        dni=data['dni'], ghi=data['ghi'], dhi=data['dhi']
    )['poa_diffuse']
    
    aoi = static_flatplate_sys.get_aoi(
        solar_zenith=solpos['zenith'],
        solar_azimuth=solpos['azimuth'],
    )
    
    pv_effective_irradiance = (
        static_flatplate_sys.get_irradiance(
        solar_zenith=solpos['zenith'], solar_azimuth=solpos['azimuth'],
        aoi=aoi, aoi_limit=aoi_limit,
        dni=data['dni'], ghi=data['ghi'], dhi=data['dhi'])
    ) * pvlib.iam.martin_ruiz(aoi, a_r=0.16) + cpv_irradiance_spillage
    
    pv_cell_temp = static_flatplate_sys.pvsyst_celltemp(
        poa_flatplate_static=pv_effective_irradiance,
        temp_air=data['temp_air'],
        wind_speed=data['wind_speed']
    )
    
    pv_diode_parameters = static_flatplate_sys.calcparams_pvsyst(
        effective_irradiance=pv_effective_irradiance,
        temp_cell=pv_cell_temp,
    )
    
    pv_power = static_flatplate_sys.singlediode(*pv_diode_parameters)

    return pv_power, aoi