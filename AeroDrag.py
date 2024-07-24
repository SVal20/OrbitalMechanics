from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts
import tools as t

tspan = 24.0*3600
dt  = 10.0
cb=pd.earth

if __name__ == '__main__':
    perts = null_perts()
    perts['aero'] = True
    perts['Cd']=2.2
    perts['A']=(1e-3)**2/4.0

    mass0 = 10.0

    rp = 215+cb['radius']
    ra = 300+cb['radius']
    
    raan = 340.0
    i = 65.2
    aop = 58.0
    ta = 332.0

    a = (rp+ra)/2.0
    e = (ra-rp)/(ra+rp)

    state0 = [a,e,i,ta,aop,raan]

    op = OP(state0, tspan, dt, deg = True, coes=True, mass0 = mass0, perts = perts)
    op.plot_alts(show_plot=True, hours=True)
    op.plot_3d(show_plot=True)
    op.calculate_OEs()
    op.plot_OEs(show_plot=True, hours=True)
