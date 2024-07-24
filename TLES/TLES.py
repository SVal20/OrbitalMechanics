from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import tools as t

cb = pd.earth
tspan = 24.0 * 3600
dt = 100.0

if __name__ == '__main__':

    #ISS
    ISS = OP(t.TLES2OES('ISS.txt'), tspan, dt,deg=False)

    #Hubble
    Hubble = OP(t.TLES2OES('Hubble.txt'), tspan, dt, deg=False)

    t.plot_n_orbits([ISS.rs,Hubble.rs],['ISS','Hubble'],show_plot=True,title = 'ISS and Hubble')