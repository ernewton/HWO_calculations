# sbin/sbin.py

import numpy as np
import pandas as pd

from sbin import parameters




def suppression_factor_simple(a_values):
    a_inner_true = 10  # AU (suppression 100%)
    a_outer_true = 200  # AU (suppression 0%)
    results = (np.log10(a_values) - np.log10(a_inner_true)) / (np.log10(a_outer_true) - np.log10(a_inner_true))
    return np.clip(results, a_min=0, a_max=1)

def suppression_factor(a_values):
    """
    Return the suppression factor S for a set of orbital semi-major axes.
    
    S is defined as a piece-wise linear function (linear in log space of 
    semi-major axis) by three break points as specified in Moe \& Kratter 
    (2021).

    Parameters
    ----------
    a_values : array‑like
        One‑dimensional collection of semi‑major axes (in astronomical units).
 
    Returns
    -------
    np.ndarray
        An array of the same shape as ``a_values`` containing the suppression
        factor for each input value, expressed as a fraction from [0,1] where
        0 means complete suppression and 1 means no suppression.
    """
    
    # Change to log space for semi-major axis
    log_a = np.log10(np.asarray(a_values)) # Semi-major axes for which to calculate S
    
    # Define the function
    log_points = np.array([np.log10(1), np.log10(10), np.log10(200)]) # break points in log(au)
    S_points = np.array([0, 0.15, 1])  # suppression value at break points
    
    # Apply the formula
    S_val = np.interp(log_a, log_points, S_points)
    
    return np.clip(S_val, a_min=0, a_max=1)



def suppression_simulation(planets_cat, separations=None, join_on='KOI', prad_col='koi_prad'):
    """
    Simulate the suppression of planet formation by stellar binaries. Return 
    the original (single-star host) planet population from ``planets_cat'',
    suppressed as if all the host stars are binaries with semi-major axes 
    drawn from the distribution defined by ``suppression_cat''.

    Parameters
    ----------
    planets_cat : pandas.DataFrame
        Table of single‑star host planets (must contain the radius column 
        given by ``prad_col`` and the join column given by ``join_on``).
    separations : np.array, sort of optional
        Array-like list of binary star separations from which to draw, in 
        units of astronomical units.
    join_on : str, optional
        Column name used to merge the two catalogs (default ``'KOI'``).
    prad_col : str, optional
        Column name for planet radius in ``planets_cat`` (default ``'Rp'``).

    Returns
    -------
    planet_radius : np.ndarray
        Radii of planets that survive the suppression process.
    planet_counts : pandas.DataFrame
        Number of surviving planets per KOI (columns ``'KOI'`` and 
        ``'n_planets'``).
    frac_super_earths :
        Fraction of surviving planets with radius < ``rad_valley``.
    frac_multiplanet : float
        Fraction of surviving KOIs that are multiplanet systems.
    """
    
    # ----------------------------------------------
    # Set-up
    # ----------------------------------------------

    # Number of UNIQUE stellar hosts
    d = {'KOI': planets_cat['KOI'].unique()}
    suppression_cat = pd.DataFrame(data=d)

    n_stars = len(suppression_cat)
    n_planets = len(planets_cat)

    # ----------------------------------------------
    # Assign a binary star to each planetary system
    # ----------------------------------------------   
    
    # Randomly draw a binary separation for each stellar host
    if separations is not None:
        random_separations = np.random.choice(separations, 
                                          size=n_stars, replace=True)
    else:
        raise()
        
    # Add the jitter to the sampled values
    error_std = 0.1 * random_separations
    random_error = np.random.normal(loc=0.0, scale=error_std)
    random_separations = random_separations + random_error

    # Suppress planet formation (per STAR)
    suppression_cat['my_factor'] = suppression_factor(random_separations)

    # Match STAR suppression to each PLANET
    # ``realization'' will hold the outcome of this simulation
    realization = planets_cat.merge(suppression_cat, on=join_on)
    realization['planet_exists'] = np.zeros(n_planets, dtype=bool) 

    # ----------------------------------------------
    # Suppress planets
    # ----------------------------------------------   

    # Condition 1: If Rp < 1.8, planet formation is not suppressed
    realization.loc[realization[prad_col] <= parameters.radius_valley, 'planet_exists'] = True

    # Condition 2: If Rp >= 1.8, planet formation probabilistically suppressed
    mask = (realization[prad_col] > parameters.radius_valley)
    random_vals = np.random.rand(n_planets)  # uniform random [0,1) to compare to my_factor
    realization.loc[mask, 'planet_exists'] = random_vals[mask] < realization['my_factor'][mask]

    # ----------------------------------------------
    # Determine the properties of the suppressed population
    # ----------------------------------------------   
    
    # Only the planets that still exist
    obs = realization.loc[realization['planet_exists'] == 1].copy()
    planet_radius = obs[prad_col]

    n_super_earths_after = float(len(planet_radius[planet_radius < parameters.radius_valley]))
    n_planets_after = float(len(planet_radius))
    
    planet_counts = obs.groupby('KOI').size().reset_index(name='n_planets')    
    mtps = float(len(planet_counts.loc[planet_counts['n_planets']>1]))
    stps = float(len(planet_counts.loc[planet_counts['n_planets']==1]))
    
    
    return(planet_radius, planet_counts, 
           n_super_earths_after/n_planets_after, 
           mtps/(mtps+stps))