import numpy as np
from math import sqrt
from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import tools as t

tspan = 30.0*3600
dt = 100.0

if __name__ == '__main__':
    r_mag = pd.earth['radius']+400.0
    v_mag = sqrt(pd.earth['mu']/r_mag)

    r0 = np.array([r_mag,0,0])
    v0 = np.array([0,v_mag,0])
    
    r_mag0 = pd.earth['radius']+1000.0
    v_mag0 = sqrt(pd.earth['mu']/r_mag)*1.3

    r00 = np.array([r_mag0,0,0])
    v00 = np.array([0,v_mag0,0.3])

    op0 = OP(r0,v0,tspan,dt)
    op00 = OP(r00,v00,tspan,dt)
    
    op0.propagate_orbit()
    op00.propagate_orbit()

    t.plot_n_orbits([op0.rs,op00.rs],['Orbit 1','Orbit 2'],show_plot=True,title='Two orbits')