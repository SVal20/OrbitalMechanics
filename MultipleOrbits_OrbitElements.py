from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')
import planetary_data as pd
from OrbitPropagator import OrbitPropagator as OP
import tools as t

tspan = 24.0*3600
dt  = 10.0
cb=pd.earth

if __name__ == '__main__':

    c0 = [80000,1.4,30,30,60,40]
    c00 = [(411+425)/2+cb['radius'], 0.0010883, 51.6394, 0.0, 344.1821, 303.6364]

    op0 = OP(c0,tspan,dt,False)
    op00 = OP(c00,tspan,dt)
    

    t.plot_n_orbits([op0.rs,op00.rs],['Hyperbolic trajectory','International Space Station'],show_plot=True,title="Example 4.7 Orbital mechanics & ISS")