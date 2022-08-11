import sys
import os
import copy as cp
os.environ["OMP_NUM_THREADS"] = "1"
import numpy as np
from scipy.interpolate import interp1d,CubicSpline
from petitRADTRANS import nat_cst as nc
from petitRADTRANS.retrieval import cloud_cond as fc
from .util import surf_to_meas, calc_MMW
"""
Models Module

This module contains a set of functions that generate the spectra used
in the petitRADTRANS retrieval. This includes setting up the
pressure-temperature structure, the chemistry, and the radiative
transfer to compute the emission or transmission spectrum.

All models must take the same set of inputs:

    pRT_object : petitRADTRANS.RadTrans
        This is the pRT object that is used to compute the spectrum
        It must be fully initialized prior to be used in the model function
    parameters : dict
        A dictionary of Parameter objects. The naming of the parameters
        must be consistent between the Priors and the model function you
        are using.
    PT_plot_mode : bool
        If this argument is True, the model function should return the pressure
        and temperature arrays before computing the flux.
    AMR : bool
        If this parameter is True, your model should allow for reshaping of the
        pressure and temperature arrays based on the position of the clouds or
        the location of the photosphere, increasing the resolution where required.
        For example, using the fixed_length_amr function defined below.
"""
# Global constants to reduce calculations and initializations.
PGLOBAL = np.logspace(-6,3,1000)

def emission_model_diseq(pRT_object,
                         parameters,
                         PT_plot_mode = False,
                         AMR = True):
    """
    Disequilibrium Chemistry Emission Model

    This model computes an emission spectrum based on disequilibrium carbon chemistry,
    equilibrium clouds and a spline temperature-pressure profile. (Molliere 2020).

    Args:
        pRT_object : object
            An instance of the pRT class, with optical properties as defined in the RunDefinition.
        parameters : dict
            Dictionary of required parameters:
                *  D_pl : Distance to the planet in [cm]
                *  log_g : Log of surface gravity
                *  R_pl : planet radius [cm]
                *  T_int : Interior temperature of the planet [K]
                *  T3 : Innermost temperature spline [K]
                *  T2 : Middle temperature spline [K]
                *  T1 : Outer temperature spline [K]
                *  alpha : power law index in tau = delta * press_cgs**alpha
                *  log_delta : proportionality factor in tau = delta * press_cgs**alpha
                *  sigma_lnorm : Width of cloud particle size distribution (log normal)
                *  log_pquench : Pressure at which CO, CH4 and H2O abundances become vertically constant
                *  Fe/H : Metallicity
                *  C/O : Carbon to oxygen ratio
                *  log_kzz : Vertical mixing parameter
                *  fsed : sedimentation parameter
                *  log_X_cb : Scaling factor for equilibrium cloud abundances.
        PT_plot_mode : bool
            Return only the pressure-temperature profile for plotting. Evaluate mode only.
        AMR :
            Adaptive mesh refinement. Use the high resolution pressure grid around the cloud base.

    Returns:
        wlen_model : np.array
            Wavlength array of computed model, not binned to data [um]
        spectrum_model : np.array
            Computed emission spectrum [W/m2/micron]
    """

    pglobal_check(pRT_object.press,
                    parameters['pressure_simple'].value,
                    parameters['pressure_scaling'].value)
    #for key, val in parameters.items():
    #    print(key,val.value)

    # Priors for these parameters are implemented here, as they depend on each other
    T3 = ((3./4.*parameters['T_int'].value**4.*(0.1+2./3.))**0.25)*(1-parameters['T3'].value)
    T2 = T3*(1.0-parameters['T2'].value)
    T1 = T2*(1.0-parameters['T1'].value)
    delta = ((10**(-3.0+5.0*parameters['log_delta'].value))*1e6)**(-parameters['alpha'].value)

    # Make the P-T profile
    temp_arr = np.array([T1,T2,T3])
    if AMR:
        p_use = PGLOBAL
    else:
        p_use = pRT_object.press/1e6
    temperatures = PT_ret_model(temp_arr, \
                            delta,
                            parameters['alpha'].value,
                            parameters['T_int'].value,
                            p_use,
                            parameters['Fe/H'].value,
                            parameters['C/O'].value,
                            conv=True)
    if PT_plot_mode:
        return p_use, temperatures

    # If in evaluation mode, and PTs are supposed to be plotted
    abundances, MMW, small_index = get_abundances(p_use,
                                                  temperatures,
                                                  pRT_object.line_species,
                                                  pRT_object.cloud_species,
                                                  parameters,
                                                  AMR =AMR)
    Kzz_use = (10**parameters['log_kzz'].value ) * np.ones_like(p_use)

    # Only include the high resolution pressure array near the cloud base.
    pressures = p_use
    if AMR:
        temperatures = temperatures[small_index]
        pressures = PGLOBAL[small_index]
        MMW = MMW[small_index]
        Kzz_use = Kzz_use[small_index]
        #pRT_object.setup_opa_structure(pressures)
    # Calculate the spectrum
    if pressures.shape[0] != pRT_object.press.shape[0]:
        return None,None
    pRT_object.press = pressures*1e6

    gravity = 2.0
    R_pl = 1.0
    if 'log_g' in parameters.keys() and 'mass' in parameters.keys():
        gravity = 10**parameters['log_g'].value
        R_pl = np.sqrt(nc.G*parameters['mass'].value/gravity)
    elif 'log_g' in parameters.keys():
        gravity= 10**parameters['log_g'].value
        R_pl = parameters['R_pl'].value
    elif 'mass' in parameters.keys():
        R_pl = parameters['R_pl'].value
        gravity = nc.G * parameters['mass'].value/R_pl**2
    else:
        print("Pick two of log_g, R_pl and mass priors!")
        sys.exit(5)
    pRT_object.calc_flux(temperatures,
                        abundances,
                        gravity,
                        MMW,
                        contribution = False,
                        fsed = parameters['fsed'].value,
                        Kzz = Kzz_use,
                        sigma_lnorm = parameters['sigma_lnorm'].value)
    # Getting the model into correct units (W/m2/micron)
    wlen_model = nc.c/pRT_object.freq/1e-4
    wlen = nc.c/pRT_object.freq
    f_lambda = pRT_object.flux*nc.c/wlen**2.
    # convert to flux per m^2 (from flux per cm^2) cancels with step below
    #f_lambda = f_lambda * 1e4
    # convert to flux per micron (from flux per cm) cancels with step above
    #f_lambda = f_lambda * 1e-4
    # convert from ergs to Joule
    f_lambda = f_lambda * 1e-7
    spectrum_model = surf_to_meas(f_lambda,
                                  R_pl,
                                  parameters['D_pl'].value)    #print(wlen_model,spectrum_model)
    return wlen_model, spectrum_model

def guillot_free_emission(pRT_object, \
                            parameters, \
                            PT_plot_mode = False,
                            AMR = False):
    """
    Free Chemistry Emission Model

    This model computes an emission spectrum based on free retrieval chemistry,
    free Ackermann-Marley clouds and a Guillot temperature-pressure profile. (Molliere 2018).

    Args:
        pRT_object : object
            An instance of the pRT class, with optical properties as defined in the RunDefinition.
        parameters : dict
            Dictionary of required parameters:
                *  D_pl : Distance to the planet in [cm]
                *  log_g : Log of surface gravity
                *  R_pl : planet radius [cm]
                *  T_int : Interior temperature of the planet [K]
                *  T_equ : Equilibrium temperature of the planet
                *  gamma : Guillot gamma parameter
                *  log_kappa_IR : The log of the ratio between the infrared and optical opacities
                *  sigma_lnorm : Width of cloud particle size distribution (log normal)
                *  log_kzz : Vertical mixing parameter
                *  fsed : sedimentation parameter
                *  species : Log abundances for each species in rd.line_list
                   (species stands in for the actual name)
                *  log_X_cb : Log cloud abundances.
                *  Pbase : log of cloud base pressure for each species.
        PT_plot_mode : bool
            Return only the pressure-temperature profile for plotting. Evaluate mode only.
        AMR :
            Adaptive mesh refinement. Use the high resolution pressure grid around the cloud base.

    Returns:
        wlen_model : np.array
            Wavlength array of computed model, not binned to data [um]
        spectrum_model : np.array
            Computed emission spectrum [W/m2/micron]
    """

    #for key, val in parameters.items():
    #    print(key,val.value)

    pglobal_check(pRT_object.press,
                  parameters['pressure_simple'].value,
                  parameters['pressure_scaling'].value)
    if AMR:
        p_use = PGLOBAL
    else:
        p_use = pRT_object.press/1e6
    temperatures = nc.guillot_global(p_use, \
                                10**parameters['log_kappa_IR'].value,
                                parameters['gamma'].value, \
                                10**parameters['log_g'].value, \
                                parameters['T_int'].value, \
                                parameters['T_equ'].value)
    Pbases = {}
    # TODO - identify species rather than hard coding
    Pbases['Fe(c)'] = fc.simple_cdf_Fe_free(p_use, temperatures,
                                10**parameters['log_X_cb_Fe(c)'].value)
    Pbases['MgSiO3(c)'] = fc.simple_cdf_MgSiO3_free(p_use, temperatures,
                                10**parameters['log_X_cb_MgSiO3(c)'].value)

    if AMR:
        p_clouds = np.array(list(Pbases.values()))
        pressures,small_index = fixed_length_amr(p_clouds,
                                                 p_use,
                                                 parameters['pressure_scaling'].value,
                                                 parameters['pressure_width'].value)
        pRT_object.press = pressures * 1e6
        temperatures = temperatures[small_index]
    else:
        pressures = pRT_object.press/1e6

    # If in evaluation mode, and PTs are supposed to be plotted
    if PT_plot_mode:
        return pressures, temperatures
    abundances = {}
    msum = 0.0
    for species in pRT_object.line_species:
        abundances[species] = 10**parameters[species.split("_R_")[0]].value * np.ones_like(pressures)
        msum += 10**parameters[species.split("_R_")[0]].value
    abundances['H2'] = 0.766 * (1.0-msum) * np.ones_like(pressures)
    abundances['He'] = 0.234 * (1.0-msum) * np.ones_like(pressures)

    MMW = calc_MMW(abundances)
    for cloud in pRT_object.cloud_species:
        cname = cloud.split('_')[0]
        pbase = Pbases[cname]
        abundances[cname] = np.zeros_like(temperatures)
        msum += 10**parameters['log_X_cb_'+cname].value
        try:
            abundances[cname][pressures < pbase] = \
                            10**parameters['log_X_cb_'+cname].value *\
                            ((pressures[pressures <= pbase]/pbase)**parameters['fsed'].value)
        except:
            return None,None

    # This is bad for some reason...
    #if msum > 1.0:
    #    return None, None
    gravity = 2.0
    R_pl = 1.0
    if 'log_g' in parameters.keys() and 'mass' in parameters.keys():
        gravity = 10**parameters['log_g'].value
        R_pl = np.sqrt(nc.G*parameters['mass'].value/gravity)
    elif 'log_g' in parameters.keys():
        gravity= 10**parameters['log_g'].value
        R_pl = parameters['R_pl'].value
    elif 'mass' in parameters.keys():
        R_pl = parameters['R_pl'].value
        gravity = nc.G * parameters['mass'].value/R_pl**2
    else:
        print("Pick two of log_g, R_pl and mass priors!")
        sys.exit(5)
    pRT_object.calc_flux(temperatures, \
                     abundances, \
                     gravity, \
                     MMW, \
                     contribution = False,
                     fsed = parameters['fsed'].value,
                     Kzz = 10**parameters['log_kzz'].value * np.ones_like(pressures),
                     sigma_lnorm = parameters['sigma_lnorm'].value)

    wlen_model = nc.c/pRT_object.freq/1e-4
    wlen = nc.c/pRT_object.freq
    f_lambda = pRT_object.flux*nc.c/wlen**2.
    # convert to flux per m^2 (from flux per cm^2) cancels with step below
    #f_lambda = f_lambda * 1e4
    # convert to flux per micron (from flux per cm) cancels with step above
    #f_lambda = f_lambda * 1e-4
    # convert from ergs to Joule
    f_lambda = f_lambda * 1e-7
    spectrum_model = surf_to_meas(f_lambda,
                                  R_pl,
                                  parameters['D_pl'].value)
    return wlen_model, spectrum_model

def guillot_eqchem_transmission(pRT_object, \
                                    parameters, \
                                    PT_plot_mode = False,
                                    AMR = False):
    """
    Equilibrium Chemistry Transmission Model, Guillot Profile

    This model computes a transmission spectrum based on equilibrium chemistry
    and a Guillot temperature-pressure profile.

    Args:
        pRT_object : object
            An instance of the pRT class, with optical properties as defined in the RunDefinition.
        parameters : dict
            Dictionary of required parameters:
                *  Rstar : Radius of the host star [cm]
                *  log_g : Log of surface gravity
                *  R_pl : planet radius [cm]
                *  T_int : Interior temperature of the planet [K]
                *  T_equ : Equilibrium temperature of the planet
                *  gamma : Guillot gamma parameter
                *  log_kappa_IR : The log of the ratio between the infrared and optical opacities
                *  Fe/H : Metallicity
                *  C/O : Carbon to oxygen ratio
                *  Pcloud : optional, cloud base pressure of a grey cloud deck.
        PT_plot_mode : bool
            Return only the pressure-temperature profile for plotting. Evaluate mode only.
        AMR :
            Adaptive mesh refinement. Use the high resolution pressure grid around the cloud base.

    Returns:
        wlen_model : np.array
            Wavlength array of computed model, not binned to data [um]
        spectrum_model : np.array
            Computed transmission spectrum R_pl**2/Rstar**2
    """

    try:
        from petitRADTRANS import poor_mans_nonequ_chem as pm
    except ImportError:
        print("Could not import poor_mans_nonequ_chemistry. Exiting.")
        sys.exit(2)
    # Make the P-T profile
    pressures = pRT_object.press/1e6
    temperatures = nc.guillot_global(pressures, \
                                    10**parameters['log_kappa_IR'].value, \
                                    parameters['gamma'].value, \
                                    10**parameters['log_g'].value, \
                                    parameters['Tint'].value, \
                                    parameters['Tequ'].value)

    # If in evaluation mode, and PTs are supposed to be plotted
    if PT_plot_mode:
        return pressures, temperatures

    # Make the abundance profile
    COs = parameters['C/O'].value * np.ones_like(pressures)
    FeHs = parameters['[Fe/H]'].value * np.ones_like(pressures)

    abundances_interp = pm.interpol_abundances(COs, \
                                               FeHs, \
                                               temperatures, \
                                               pressures)

    abundances = {}
    for species in pRT_object.line_species:
        abundances[species] = abundances_interp[ \
                    species.split('_R_')[0].replace('_all_iso', '').replace('C2H2','C2H2,acetylene')]
    abundances['H2'] = abundances_interp['H2']
    abundances['He'] = abundances_interp['He']

    MMW = abundances_interp['MMW']

    # Calculate the spectrum
    if 'Pcloud' in parameters.keys():
        pRT_object.calc_transm(temperatures, \
                               abundances, \
                               10**parameters['log_g'].value, \
                               MMW, \
                               R_pl=parameters['R_pl'].value, \
                               P0_bar=0.01,
                               Pcloud = parameters['Pcloud'].value)
    else:
        pRT_object.calc_transm(temperatures, \
                               abundances, \
                               10**parameters['log_g'].value, \
                               MMW, \
                               R_pl=parameters['R_pl'].value, \
                               P0_bar=0.01)

    wlen_model = nc.c/pRT_object.freq/1e-4
    spectrum_model = (pRT_object.transm_rad/parameters['Rstar'].value)**2.

    return wlen_model, spectrum_model

def isothermal_eqchem_transmission(pRT_object, \
                                    parameters, \
                                    PT_plot_mode = False,
                                    AMR = False):
    """
    Equilibrium Chemistry Transmission Model, Isothermal Profile

    This model computes a transmission spectrum based on equilibrium chemistry
    and a Guillot temperature-pressure profile.

    Args:
        pRT_object : object
            An instance of the pRT class, with optical properties as defined in the RunDefinition.
        parameters : dict
            Dictionary of required parameters:
                *  Rstar : Radius of the host star [cm]
                *  log_g : Log of surface gravity
                *  R_pl : planet radius [cm]
                *  T_int : Interior temperature of the planet [K]
                *  T_equ : Equilibrium temperature of the planet
                *  Fe/H : Metallicity
                *  C/O : Carbon to oxygen ratio
                *  Pcloud : optional, cloud base pressure of a grey cloud deck.
        PT_plot_mode : bool
            Return only the pressure-temperature profile for plotting. Evaluate mode only.
        AMR :
            Adaptive mesh refinement. Use the high resolution pressure grid around the cloud base.

    Returns:
        wlen_model : np.array
            Wavlength array of computed model, not binned to data [um]
        spectrum_model : np.array
            Computed transmission spectrum R_pl**2/Rstar**2
    """

    try:
        from petitRADTRANS import poor_mans_nonequ_chem as pm
    except ImportError:
        print("Could not import poor_mans_nonequ_chemistry. Exiting.")
        sys.exit(2)
    # Make the P-T profile
    pressures = pRT_object.press/1e6
    temperatures = parameters['Temp'].value * np.ones_like(pressures)

    # If in evaluation mode, and PTs are supposed to be plotted
    if PT_plot_mode:
        return pressures, temperatures

    # Make the abundance profile
    COs = parameters['C/O'].value * np.ones_like(pressures)
    FeHs = parameters['[Fe/H]'].value * np.ones_like(pressures)

    abundances_interp = pm.interpol_abundances(COs, \
                                               FeHs, \
                                               temperatures, \
                                               pressures)

    abundances = {}
    for species in pRT_object.line_species:
        abundances[species] = abundances_interp[ \
                    species.split('_R_')[0].replace('_all_iso', '').replace('C2H2','C2H2,acetylene')]
    abundances['H2'] = abundances_interp['H2']
    abundances['He'] = abundances_interp['He']

    MMW = abundances_interp['MMW']
    pcloud = None
    if 'log_Pcloud' in parameters.keys():
        pcloud = 10**parameters['log_Pcloud'].value
    elif 'Pcloud' in parameters.keys():
        pcloud = parameters['log_Pcloud'].value
    if pcloud is not None:
        # P0_bar is important for low gravity transmission
        # spectrum. 100 is standard, 0.01 is good for small,
        # low gravity objects
        pRT_object.calc_transm(temperatures, \
                               abundances, \
                               10**parameters['log_g'].value, \
                               MMW, \
                               R_pl=parameters['R_pl'].value, \
                               P0_bar=0.01,
                               Pcloud = pcloud)
    else:
        pRT_object.calc_transm(temperatures, \
                               abundances, \
                               10**parameters['log_g'].value, \
                               MMW, \
                               R_pl=parameters['R_pl'].value, \
                               P0_bar=0.01)
    wlen_model = nc.c/pRT_object.freq/1e-4
    spectrum_model = (pRT_object.transm_rad/parameters['Rstar'].value)**2.
    return wlen_model, spectrum_model


def isothermal_free_transmission(pRT_object, \
                                parameters, \
                                PT_plot_mode = False,
                                AMR = False):
    """
    Free Chemistry Transmission Model, Guillot Profile

    This model computes a transmission spectrum based on free retrieval chemistry
    and an isothermal temperature-pressure profile.

    Args:
        pRT_object : object
            An instance of the pRT class, with optical properties as defined in the RunDefinition.
        parameters : dict
            Dictionary of required parameters:
                *  Rstar : Radius of the host star [cm]
                *  log_g : Log of surface gravity
                *  R_pl : planet radius [cm]
                *  Temp : Isothermal temperature [K]
                *  species : Abundances for each species used in the retrieval
                *  Pcloud : optional, cloud base pressure of a grey cloud deck.
        PT_plot_mode : bool
            Return only the pressure-temperature profile for plotting. Evaluate mode only.
        AMR :
            Adaptive mesh refinement. Use the high resolution pressure grid around the cloud base.

    Returns:
        wlen_model : np.array
            Wavlength array of computed model, not binned to data [um]
        spectrum_model : np.array
            Computed transmission spectrum R_pl**2/Rstar**2
    """

    # Make the P-T profile
    pressures = pRT_object.press/1e6
    temperatures = parameters['Temp'].value * np.ones_like(pressures)
    #for key,value in parameters.items():
    #    print(key,value.value)
    # If in evaluation mode, and PTs are supposed to be plotted
    if PT_plot_mode:
        return pressures, temperatures
    abundances = {}
    msum = 0.0
    for species in pRT_object.line_species:
        spec = species.split('_R_')[0]
        abundances[species] = 10**parameters[spec].value * np.ones_like(pressures)
        msum += 10**parameters[spec].value
    if msum > 1.0:
        return None, None
    abundances['H2'] = 0.766 * (1.0-msum) * np.ones_like(pressures)
    abundances['He'] = 0.234 * (1.0-msum) * np.ones_like(pressures)

    #MMW = abundances_interp['MMW']
    MMW = calc_MMW(abundances)

    # Calculate the spectrum
    pcloud = 100.0
    if 'Pcloud' in parameters.keys():
        pcloud = parameters['Pcloud'].value
    elif 'log_Pcloud' in parameters.keys():
        pcloud = 10**parameters['log_Pcloud'].value
    # Calculate the spectrum

    if pcloud is not None:
        # P0_bar is important for low gravity transmission
        # spectrum. 100 is standard, 0.01 is good for small,
        # low gravity objects
        pRT_object.calc_transm(temperatures, \
                               abundances, \
                               10**parameters['log_g'].value, \
                               MMW, \
                               R_pl=parameters['R_pl'].value, \
                               P0_bar=0.01,
                               Pcloud = pcloud)
    else:
        pRT_object.calc_transm(temperatures, \
                               abundances, \
                               10**parameters['log_g'].value, \
                               MMW, \
                               R_pl=parameters['R_pl'].value, \
                               P0_bar=0.01)

    wlen_model = nc.c/pRT_object.freq/1e-4
    spectrum_model = (pRT_object.transm_rad/parameters['Rstar'].value)**2.
    return wlen_model, spectrum_model


##################
# Helper Functions
##################

### Global Guillot P-T formula with kappa/grav replaced by delta
def PT_ret_model(T3, delta, alpha, tint, press, FeH, CO, conv = True):
    '''
    Self-luminous retrieval P-T model.

    Arsg:
        T3 : np.array([t1, t2, t3])
            temperature points to be added on top
            radiative Eddington structure (above tau = 0.1).
            Use spline interpolation, t1 < t2 < t3 < tconnect as prior.
        delta :
            proportionality factor in tau = delta * press_cgs**alpha
        alpha:
            power law index in tau = delta * press_cgs**alpha
            For the tau model: use proximity to kappa_rosseland photosphere
            as prior.
        tint:
            internal temperature of the Eddington model
        press:
            input pressure profile in bar
        conv:
            enforce convective adiabat yes/no
        CO:
            C/O for the nabla_ad interpolation
        FeH:
            metallicity for the nabla_ad interpolation
    '''

    try:
        from petitRADTRANS import poor_mans_nonequ_chem as pm
    except ImportError:
        print("Could not import poor_mans_nonequ_chemistry. Exiting.")
        sys.exit(2)
    # Go grom bar to cgs
    press_cgs = press*1e6

    # Calculate the optical depth
    tau = delta*press_cgs**alpha

    # This is the eddington temperature
    tedd = (3./4.*tint**4.*(2./3.+tau))**0.25

    ab = pm.interpol_abundances(CO*np.ones_like(tedd), \
            FeH*np.ones_like(tedd), \
            tedd, \
            press)

    nabla_ad = ab['nabla_ad']

    # Enforce convective adiabat
    if conv:
        # Calculate the current, radiative temperature gradient
        nab_rad = np.diff(np.log(tedd))/np.diff(np.log(press_cgs))
        # Extend to array of same length as pressure structure
        nabla_rad = np.ones_like(tedd)
        nabla_rad[0] = nab_rad[0]
        nabla_rad[-1] = nab_rad[-1]
        nabla_rad[1:-1] = (nab_rad[1:]+nab_rad[:-1])/2.

        # Where is the atmosphere convectively unstable?
        conv_index = nabla_rad > nabla_ad

        # TODO: Check remains convective and convergence
        for i in range(10):
            if i == 0:
                t_take = cp.copy(tedd)
            else:
                t_take = cp.copy(tfinal)

            ab = pm.interpol_abundances(CO*np.ones_like(t_take), \
                FeH*np.ones_like(t_take), \
                t_take, \
                press)

            nabla_ad = ab['nabla_ad']

            # Calculate the average nabla_ad between the layers
            nabla_ad_mean = nabla_ad
            nabla_ad_mean[1:] = (nabla_ad[1:]+nabla_ad[:-1])/2.
            # What are the increments in temperature due to convection
            tnew = nabla_ad_mean[conv_index]*np.mean(np.diff(np.log(press_cgs)))
            # What is the last radiative temperature?
            tstart = np.log(t_take[~conv_index][-1])
            # Integrate and translate to temperature from log(temperature)
            tnew = np.exp(np.cumsum(tnew)+tstart)

            # Add upper radiative and
            # lower conective part into one single array
            tfinal = cp.copy(t_take)
            tfinal[conv_index] = tnew

            if np.max(np.abs(t_take-tfinal)/t_take) < 0.01:
                #print('n_ad', 1./(1.-nabla_ad[conv_index]))
                break

    else:
        tfinal = tedd

    # Add the three temperature-point P-T description above tau = 0.1
    def press_tau(tau):
        # Returns the pressure at a given tau, in cgs
        return (tau/delta)**(1./alpha)

    # Where is the uppermost pressure of the Eddington radiative structure?
    p_bot_spline = press_tau(0.1)

    for i_intp in range(2):

        if i_intp == 0:

            # Create the pressure coordinates for the spline support nodes at low pressure
            support_points_low = np.logspace(np.log10(press_cgs[0]), \
                             np.log10(p_bot_spline), \
                             4)

            # Create the pressure coordinates for the spline support nodes at high pressure,
            # the corresponding temperatures for these nodes will be taken from the
            # radiative+convective solution
            support_points_high = 1e1**np.arange(np.log10(p_bot_spline),
                                                 np.log10(press_cgs[-1]),
                                                 np.diff(np.log10(support_points_low))[0])

            # Combine into one support node array, don't add the p_bot_spline point twice.
            support_points = np.zeros(len(support_points_low)+len(support_points_high)-1)
            support_points[:4] = support_points_low
            support_points[4:] = support_points_high[1:]

        else:

            # Create the pressure coordinates for the spline support nodes at low pressure
            support_points_low = np.logspace(np.log10(press_cgs[0]), \
                             np.log10(p_bot_spline), \
                             7)

            # Create the pressure coordinates for the spline support nodes at high pressure,
            # the corresponding temperatures for these nodes will be taken from the
            # radiative+convective solution
            support_points_high = np.logspace(np.log10(p_bot_spline), np.log10(press_cgs[-1]), 7)

            # Combine into one support node array, don't add the p_bot_spline point twice.
            support_points = np.zeros(len(support_points_low)+len(support_points_high)-1)
            support_points[:7] = support_points_low
            support_points[7:] = support_points_high[1:]

        # Define the temperature values at the node points.
        t_support = np.zeros_like(support_points)

        if i_intp == 0:
            tfintp = interp1d(press_cgs, tfinal,kind='cubic')
            # The temperature at p_bot_spline (from the radiative-convectice solution)
            t_support[int(len(support_points_low))-1] = tfintp(p_bot_spline)
            # The temperature at pressures below p_bot_spline (free parameters)
            t_support[:(int(len(support_points_low))-1)] = T3
            # t_support[:3] = tfintp(support_points_low)
            # The temperature at pressures above p_bot_spline
            # (from the radiative-convectice solution)
            t_support[int(len(support_points_low)):] = \
                tfintp(support_points[(int(len(support_points_low))):])

        else:
            tfintp1 = interp1d(press_cgs, tret,kind='cubic')
            t_support[:(int(len(support_points_low))-1)] = \
                tfintp1(support_points[:(int(len(support_points_low))-1)])

            tfintp = interp1d(press_cgs, tfinal)
            # The temperature at p_bot_spline (from the radiative-convectice solution)
            t_support[int(len(support_points_low))-1] = tfintp(p_bot_spline)
            #print('diff', t_connect_calc - tfintp(p_bot_spline))
            t_support[int(len(support_points_low)):] = \
                tfintp(support_points[(int(len(support_points_low))):])

        # Make the temperature spline interpolation to be returned to the user
        cs = CubicSpline(np.log10(support_points), t_support)
        tret = cs(np.log10(press_cgs))

    tret[tret<0.0] = 10.0
    # Return the temperature, the pressure at tau = 1,
    # and the temperature at the connection point.
    # The last two are needed for the priors on the P-T profile.
    return tret#, press_tau(1.)/1e6, tfintp(p_bot_spline)

def _make_half_pressure_better(P_clouds, press):
    """
    deprecated
    """
    # Obsolete, replaced with fixed_length_amr
    press_plus_index = np.zeros((press.shape[0],2))
    press_plus_index[:,0] = press
    press_plus_index[:,1] = range(len(press))

    press_small = press_plus_index[::24, :]
    press_plus_index = press_plus_index[::2,:]

    indexes_small = press_small[:,0] > 0.
    indexes       = press_plus_index[:,0] > 0.

    for cname,P_cloud in P_clouds.items():
        indexes_small = indexes_small & \
            ((np.log10(press_small[:,0]/P_cloud) > 0.05) | \
            (np.log10(press_small[:,0]/P_cloud) < -0.3))
        indexes = indexes & \
            ((np.log10(press_plus_index[:,0]/P_cloud) > 0.05) | \
            (np.log10(press_plus_index[:,0]/P_cloud) < -0.3))

    press_cut = press_plus_index[~indexes, :]
    press_small_cut = press_small[indexes_small, :]

    press_out = np.zeros((len(press_cut)+len(press_small_cut))*2).reshape((len(press_cut)+len(press_small_cut)), 2)
    press_out[:len(press_small_cut), :] = press_small_cut
    press_out[len(press_small_cut):, :] = press_cut

    press_out = np.sort(press_out, axis = 0)
    return press_out[:,0],  press_out[:, 1].astype('int')

def fixed_length_amr(P_clouds, press, scaling = 10, width = 3):
    """
    This function takes in the cloud base pressures for each cloud,
    and returns an array of pressures with a high resolution mesh
    in the region where the cloud is located.

    Args:
        P_clouds : numpy.ndarray
            The cloud base pressures in bar
        press : np.ndarray
            The high resolution pressure array.
        scaling : int
            The factor by which the low resolution pressure array is scaled
        width : int
            The number of low resolution bins to be replaced for each cloud layer.
    """

    # P_clouds is array of pressures
    # press should be ~len scaling*100
    # guarantees total length will be press.shape[0] + P_clouds.shape[0]*width*(scaling - 1)
    # scaling is how many hi-res points per normal point. Must be int
    # width is the number of low res points to replace. Must be int.

    press_plus_index = np.zeros((press.shape[0],2))
    press_plus_index[:,0] = np.logspace(np.log10(np.min(press)),np.log10(np.max(press)),press.shape[0])
    press_plus_index[:,1] = range(len(press_plus_index[:,0]))
    # Set up arrays for indexing
    press_small= press_plus_index[::scaling,:]
    # Make some lists to store the replacement indices
    c_list = []
    for i,P_cloud in enumerate(P_clouds):
        # Find out where the clouds are in the high res grid
        idx = (np.abs(press_plus_index[:,0] - P_cloud)).argmin()
        # constant length list of indices around that point
        inds = np.linspace(int(idx-(width*scaling/2.0)),int(idx+(width*scaling/2.0))-1,int(scaling*width),dtype=int)
        c_list.append(inds)
    # We need to return a list that's always the same length
    # So we need to check for duplicates
    total_inds = []
    for j in range(len(c_list)):
        # At first, just copy in the list
        if j == 0:
            if np.any(c_list[j] < 0):
                total_inds.extend(np.linspace(0,len(c_list[j]),len(c_list[j]),dtype=int))
                continue
            if np.any(c_list[j] > press.shape[0]):
                total_inds.extend(np.linspace(press.shape[0]-len(c_list[j]),press.shape[0],dtype=int))
                continue
            total_inds.extend(c_list[j])
            continue
        # Check if the next set of indices is lower than the current minimum
        # if so, we want to add scaling*width indices below the current minimum
        if min(c_list[j])<=min(np.array(total_inds)):
            start = c_list[j][-1]
            sl = len(total_inds)
            ind = 0
            done = 0
            while len(total_inds) < sl+int(scaling*width):
                if start-ind < 0:
                    start = start + len(c_list[j])
                    ind = done
                if np.in1d(start-ind,np.array(total_inds)).any():
                    ind += 1
                    continue
                else:
                    total_inds.append(int(start-ind))
                    ind+=1
                    done += 1
        # Check if the smallest new index is larger than the current max
        # if so, we can just add the indexes
        # I can probably replace all this with total_inds.extend(c_list[j])
        elif max(c_list[j])>=max(np.array(total_inds)):
            start = c_list[j][0]
            sl = len(total_inds)
            ind = 0
            done = 0
            while len(total_inds) < sl+int(scaling*width):
                if (start+ind) >= (len(press_plus_index)-1):
                    start = start - len(c_list[j])
                    ind = done
                    continue
                if np.in1d(start+ind,np.array(total_inds)).any():
                    ind += 1
                    continue
                else:
                    total_inds.append(int(start+ind))
                    ind+=1
                    done += 1
        else:
            # This loop takes care of cases where we're between existing entries
            # it adds indices until duplicates are found, then keeps incrementing
            # until there is a free index to add.
            start = c_list[j][0]
            sl = len(total_inds)
            ind = 0
            done = 0
            while len(total_inds) < sl+int(scaling*width):
                if (start+ind) >= (len(press_plus_index)-1):
                    start = start - len(c_list[j])
                    ind = done
                    continue
                if np.in1d(start+ind,np.array(total_inds)).any():
                    ind += 1
                    continue
                else:
                    total_inds.append(int(start+ind))
                    ind+=1
                    done += 1
    #total_inds = np.array(sorted(total_inds,reverse=False))
    # Stack the low res and high res grids, sort it, and take the unique values
    try:
        press_out = np.vstack((press_small,press_plus_index[total_inds]))
    except:
        print("AMR returned incorrect length")
        return PGLOBAL, np.array([0])
    press_out = np.sort(press_out, axis = 0)
    p_out,ind = np.unique(press_out[:,0],return_index = True)
    return p_out,  press_out[ind, 1].astype('int')


def get_abundances(pressures, temperatures, line_species, cloud_species, parameters, AMR = False):
    """
    This function takes in the C/O ratio, metallicity, and quench pressures and uses them
    to compute the gas phase and equilibrium condensate abundances.
    Clouds are currently hard coded into the function.

    Args:
        pressures : numpy.ndarray
            A log spaced pressure array. If AMR is on it should be the full high resolution grid.
        temperatures : numpy.ndarray
            A temperature array with the same shape as pressures
        line_species : List(str)
            A list of gas species that will contribute to the line-by-line opacity of the pRT atmosphere.
        cloud_species : List(str)
            A list of condensate species that will contribute to the cloud opacity of the pRT atmosphere.
        parameters : dict
            A dictionary of model parameters, in particular it must contain the names C/O, Fe/H and
            log_pquench. Additionally the cloud parameters log_X_cb_Fe(c) and MgSiO3(c) must be present.
        AMR : bool
            Turn the adaptive mesh grid on or off. See fixed_length_amr for implementation.

    Returns:
        abundances : dict
            Mass fraction abundances of all atmospheric species
        MMW : numpy.ndarray
            Array of the mean molecular weights in each pressure bin
        small_index : numpy.ndarray
            The indices of the high resolution grid to use to define the adaptive grid.
    """

    try:
        from petitRADTRANS import poor_mans_nonequ_chem as pm
    except ImportError:
        print("Could not import poor_mans_nonequ_chemistry. Exiting.")
        sys.exit(2)
    # Make the abundance profile
    COs = parameters['C/O'].value * np.ones_like(pressures)
    FeHs = parameters['Fe/H'].value * np.ones_like(pressures)

    # Prior check all input params
    XFe = fc.return_XFe(parameters['Fe/H'].value, parameters['C/O'].value)
    XMgSiO3 = fc.return_XMgSiO3(parameters['Fe/H'].value, parameters['C/O'].value)

    clouds = {}
    clouds['Fe(c)'] = 10**parameters['log_X_cb_Fe(c)'].value*XFe
    clouds['MgSiO3(c)'] = 10**parameters['log_X_cb_MgSiO3(c)'].value*XMgSiO3

    pquench_C = None
    if 'log_pquench' in parameters.keys():
        pquench_C = 10**parameters['log_pquench'].value
    abundances_interp = pm.interpol_abundances(COs, \
                                               FeHs, \
                                               temperatures, \
                                               pressures,
                                               Pquench_carbon = pquench_C)
    MMW = abundances_interp['MMW']

    Pbases = {}
    Pbases['Fe(c)'] = fc.simple_cdf_Fe(pressures, temperatures,
                                parameters['Fe/H'].value, parameters['C/O'].value, np.mean(MMW))
    Pbases['MgSiO3(c)'] = fc.simple_cdf_MgSiO3(pressures, temperatures,
                                parameters['Fe/H'].value, parameters['C/O'].value, np.mean(MMW))
    #Pbases['KCL(c)'] = fc.simple_cdf_KCL(pressures, temperatures,
    #                            parameters['Fe/H'].value, parameters['C/O'].value, np.mean(MMW))
    #Pbases['Na2S(c)'] = fc.simple_cdf_Na2S(pressures, temperatures,
    #                            parameters['Fe/H'].value, parameters['C/O'].value, np.mean(MMW))
    fseds = {}
    abundances = {}
    # Clouds
    p_clouds = []
    for key, val in Pbases.items():
        p_clouds.append(val)
    p_clouds = np.array(p_clouds)

    if AMR:
        press_use, small_index = fixed_length_amr(p_clouds,
                                                  pressures,
                                                  parameters['pressure_scaling'].value,
                                                  parameters['pressure_width'].value)
    else :
        #TODO: Test
        press_use = pressures
        small_index = np.linspace(press_use[0],press_use[-1],press_use.shape[0],dtype = int)

    for cloud in cp.copy(cloud_species):
        cname = cloud.split('_')[0]
        if 'log_X_cb_'+cname not in parameters.keys():
            continue
        abundances[cname] = np.zeros_like(temperatures)
        abundances[cname][pressures < Pbases[cname]] = \
                        clouds[cname] *\
                        (pressures[pressures <= Pbases[cname]]/\
                        Pbases[cname])**parameters['fsed'].value
        fseds[cname] = parameters['fsed'].value

    if AMR:
        abundances['Fe(c)'] = abundances['Fe(c)'][small_index]
        abundances['MgSiO3(c)'] = abundances['MgSiO3(c)'][small_index]
        for species in line_species:
            abundances[species] = abundances_interp[species.split('_')[0]][small_index]
        abundances['H2'] = abundances_interp['H2'][small_index]
        abundances['He'] = abundances_interp['He'][small_index]

    else:
        for species in line_species:
            abundances[species] = abundances_interp[species.split('_')[0]]
        abundances['H2'] = abundances_interp['H2']
        abundances['He'] = abundances_interp['He']
    if 'FeH' in abundances.keys():
        abundances['FeH'] = abundances['FeH']/2.
    return abundances,MMW,small_index

def pglobal_check(press,shape,scaling):
    """
    Check to ensure that the global pressure array has the correct length.
    Updates PGLOBAL.

    Args:
        press : numpy.ndarray
            Pressure array from a pRT_object. Used to set the min and max values of PGLOBAL
        shape : int
            the shape of the pressure array if no AMR is used
        scaling :
            The factor by which the pressure array resolution should be scaled.
    """

    global PGLOBAL
    if PGLOBAL.shape[0] != int(scaling*shape):
        PGLOBAL = np.logspace(np.log10(press[0]),
                              np.log10(press[-1]),
                              int(scaling*shape))

def set_resolution(lines,abundances,resolution):
    """
    deprecated
    """
    # Set correct key names in abundances for pRT, with set resolution
    # Only needed for free chemistry retrieval
    #print(lines)
    #print(abundances)
    if resolution is None:
        return abundances
    for line in lines:
        abundances[line] = abundances[line.split("_R_"+str(resolution))[0]]
        del abundances[line.split("_R_"+str(resolution))]
    return abundances
