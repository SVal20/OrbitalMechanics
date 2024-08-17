from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')

from OrbitPropagator import OrbitPropagator as OP
from OrbitPropagator import null_perts

import planetary_data as pd
import tools as t
import numpy as np

tspan = 3600 * 24 * 3.0
dt = 100.0

cb =pd.earth

date0 = '2024-04-03'

if __name__ == '__main__':

    r0_apogee, v0_apogee = t.OrbElems2rv([cb['radius']+26600, 0.74, 35.0, 180.0, 0.0, 0.0])
    r0_perigee, v0_perigee = t.OrbElems2rv([cb['radius']+26600,0.74,35.0,0.0,0.0,0.0])
    OrbElems_rotated = [cb['radius']+26600, 0.74,35.0, 0.0,45.0,0.0]

    Vapo_circular = (pd.earth['mu']/np.linalg.norm(r0_apogee))**0.5
    Vperi_circular = (pd.earth['mu']/np.linalg.norm(r0_perigee))**0.5

    vel0apo_normed = t.normed(v0_apogee)
    vel0peri_normed = t.normed(v0_perigee)

    escape_vel_apo = t.escape_velocity(np.linalg.norm(r0_apogee))
    escape_vel_peri = t.escape_velocity(np.linalg.norm(r0_perigee))

    v0_apogee_norm = np.linalg.norm(v0_apogee)
    v0_perigee_norm = np.linalg.norm(v0_perigee)

    v0apo_escape_vel = vel0apo_normed * escape_vel_apo
    v0peri_escape_vel = vel0peri_normed * escape_vel_peri

    print('Current vel at apogee:\t%.2f km/s' % v0_apogee_norm)
    print('Circular vel at apogee:\t%.2f km/s' % Vapo_circular)
    print('Escape vel at apogee:\t%.2f km/s' % escape_vel_apo)
    print('Delta V:\t%.2f km/s' % (escape_vel_apo - v0_apogee_norm))

    print()

    print('Current vel at perigee:\t%.2f km/s' % v0_perigee_norm)
    print('Circular vel at perigee:\t%.2f km/s' % Vperi_circular)
    print('Escape vel at perigee:\t%.2f km/s' % escape_vel_peri)
    print('Delta V:\t%.2f km/s' % (escape_vel_peri - v0_perigee_norm))
    
    init_state = r0_apogee.tolist() + v0_apogee.tolist()
    init_state_apogee =  r0_apogee.tolist() + v0apo_escape_vel.tolist()
    init_state_perigee = r0_perigee.tolist() + v0peri_escape_vel.tolist()
   
    
    op_original = OP(np.array(init_state), tspan, dt, coes = False)
    op_apogee = OP(np.array(init_state_apogee), tspan, dt, coes = False)
    op_perigee = OP(np.array(init_state_perigee), tspan, dt, coes = False)
    op_rotated = OP(OrbElems_rotated, tspan, dt)

    # op_perigee.calculate_OEs()
    # op_perigee.plot_OEs(show_plot=True, days= True)

    t.plot_n_orbits([op_original.rs, op_rotated.rs, op_apogee.rs, op_perigee.rs], 
                    ['Original', 'Rotated', 'At apogee', 'At perigee'], show_plot= True, title= 'Escape Trajectory')


    