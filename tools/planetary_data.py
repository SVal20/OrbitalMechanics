from numpy import array

Gm = 6.67408e-11
G = Gm * 10**-9

sun={
    'name' : 'Sun',
    'mass' : 1.989e30,
    'mu' : 1.32712e11,
    'radius' : 695700.0,
    'deorbit_altitud' : 100.0
}

atm = array([[63.096, 2.059e-4],[251.189,5.909e-11],[1000.0,3.561e-15]])
earth={
    'name' : 'Earth',
    'mass' : 5.972e24,
    'mu' : 398600,
    'radius' : 6378.0,
    'J2' : -1.082635854e-3,
    'atm_rot_vector': array([0.0,0.0,72.9211e-6]),
    'zs' :atm[:,0],
    'rhos' :atm[:,1]*10**8,
    'deorbit_altitud' : 100.0
}

venus={
    'name' : 'Venus',
    'mass' : 4.867e24,
    'mu' : 324900,
    'radius' : 108.2e6
}

moon={
    'name': 'Moon',
    'mass' : 7.34767309e22,
    'mu' : 7.34767309e22 * G,
    'radius' : 1737.1,
    'distance2earth':384400.0,
    'spice_file': 'C:/Users/santi/OneDrive/Documents/OrbitalMechanicswithPython/spice_data/de432s.bsp'
}