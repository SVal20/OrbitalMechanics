import spiceypy as spice
import numpy as np

def get_spk_binary_data(bsp_kernel_file, display = False):
    objects = spice.spkobj(bsp_kernel_file)

    ids, names, timecov_sec, timecov_cal = [],[],[],[]

    n=0
    
    if display:
        print('\nObjects in %s:' %bsp_kernel_file)

    for obj in objects:
        ids.append(obj)

        #time coverage seconds
        tcs = spice.wnfetd(spice.spkcov(bsp_kernel_file, ids[n]),0)
        #time coverage calendar
    
        tcc = [spice.timout(f, "YYYY MON DD HR:MN:SC.### (TDB) ::TDB") for f in tcs]

        timecov_sec.append(tcs)
        timecov_cal.append(tcc)

        try:
            names.append(spice.bodc2n(obj))
        except:
            names.append("Unknown name")

        if display:
            print('id: %i\t\tname: %s\t\ttc: %s --> %s' % (ids[-1],names[-1],tcc[0],tcc[1]))

        n+=1
    return ids, names, timecov_sec, timecov_cal

def timecov2array(tcs,steps):
    arr = np.zeros((steps,1))
    arr[:,0] = np.linspace(tcs[0], tcs[1], steps)
    return arr