from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')

from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

import planetary_data as pd
import tools as t
import numpy as np

tspan = 3600 * 24 * 150
dt = 100.0 
cb = pd.earth
date0 = '2020-03-28'

if __name__ == '__main__':
    perts = null_perts()

    perts['thrust'] = 0.327 # N
    perts['isp'] = 4300
    perts['thrust_direction'] = 1

    sc = {'escape_velocity': True}

    mass0 = 50.0 #kg

    state0 = [cb['radius']+800.0, 0.01, 5.0, 0.0, 0.0, 0.0]

    op = OP(state0, tspan, dt, perts= perts, mass0=mass0, sc=sc)

    op.plot_3d(show_plot=True,title='Low Thrust Escape')

    exit()

    