from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')

import numpy as np
import tools as t
import spice_tools as st
import spiceypy as spice
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP

X = 2.300823339086520E+08; Y =-3.259424153551674E+07; Z = 2.476260652060311E+06
VX=-1.016365681477186E+00; VY= 2.168007889055697E+01; VZ= 2.861694061872742E-01

state0 = np.array([X,Y,Z,VX,VY,VZ])

DATE_I = '2018-02-08'
DATE_F = '2020-03-15'

STEPS = 1000

BODIES = ['MERCURY', 'VENUS', 'EARTH', 'MARS BARYCENTER', 'ROADSTER']
FRAME = 'ECLIPJ2000'
OBSERVER = 'SOLAR SYSTEM BARYCENTER'

if __name__ == '__main__':

    spice.furnsh('spice_data/solar_system_kernel.mk')

    t0 = spice.utc2et(DATE_I)
    tf = spice.utc2et(DATE_F)

    tspan = tf-t0
    
    times = st.timecov2array([t0,tf],STEPS)

    pos = []

    for name in BODIES[:-1]:
        pos.append(np.array(spice.spkezr(name,times,FRAME,'NONE' ,OBSERVER)[0]))
 
    op0 = OP(state0, tspan, STEPS, coes=False, cb=pd.sun)
    
    pos.append(op0.rs)

    t.plot_n_orbits(pos, BODIES, show_plot=True, title= 'Tesla Roadster',cb=pd.sun)