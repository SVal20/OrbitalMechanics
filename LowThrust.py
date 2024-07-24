from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')

from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

import planetary_data as pd

tspan = 3600 * 30
dt = 100.0

cb = pd.earth

if __name__ == '__main__':

    perts = null_perts()
    perts['thrust'] = 0.327
    perts['isp'] = 4300
    perts['thrust_direction'] = -1

    sc = {'min_alt': 200}
    mass0 = 50.0 #kg

    state0 = [cb['radius'] + 1500, 0.1, 10.0, 0.0, 0.0, 0.0]

    op = OP(state0,tspan, dt, perts=perts, mass0= mass0, sc=sc)
    op.plot_alts(show_plot=True, hours= True)
    op.plot_3d(show_plot= True)
    op.calculate_OEs()
    op.plot_OEs(show_plot=True, hours=True)
