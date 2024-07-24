import numpy as np
import matplotlib.pyplot as plt 
from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP

cb=pd.earth

plt.style.use('dark_background')

if __name__=='__main__':
    r_mag = cb['radius']+500.0
    v_mag = np.sqrt(cb['mu']/r_mag)

    r0 = [r_mag,0,0]
    v0 = [0,v_mag,0]

    tspan = 3600*24.0
    dt=100.0

    op=OP(r0,v0,tspan,dt)
    op.propagate_orbit()
    op.plot_3d(True,title="Two body problem")