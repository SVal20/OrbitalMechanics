from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\tools')

import tools as t
import spice_tools as st
import planetary_data as pd
import spiceypy as spice 
import numpy as np

STEPS = 10000
FRAME = 'ECLIPJ2000'
OBSERVER = 'SUN'

if __name__ == '__main__':

    spice.furnsh('spice_data/solar_system_kernel.mk')

    ids, names, timecoverage_sec, timecoverage_cal = st.get_spk_binary_data(
        'C:/Users/santi/OneDrive/Documents/OrbitalMechanicswithPython/spice_data/de432s.bsp', display=True)

    names = [f for f in names if 'BARYCENTER' in f]

    times = st.timecov2array(timecoverage_sec, STEPS)

    PosComponents = []

    for name in names:
        PosComponents.append(np.array(spice.spkezr(name,times,FRAME,'NONE' ,OBSERVER)[0]))

    t.plot_n_orbits(PosComponents,names, cb=pd.sun,show_plot=True, title='Solar system')