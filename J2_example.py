from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import tools as t

tspan = 3600 * 24 * 20.0
dt = 100.0

cb = pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['J2'] = True

    op = OP(t.TLES2OES('TLES/ISS.txt'),tspan,dt,perts=perts)

    op.plot_3d(show_plot=True)