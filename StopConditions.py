from sys import path

path.append('C:\\Users\\santi\\OneDrive\\Documents\\OrbitalMechanicswithPython\\OrbitPropagator')

from OrbitPropagator import OrbitPropagator as OP
import planetary_data as pd

t_span = 3600 * 12.0
dt = 100.0

cb = pd.earth

if __name__ == '__main__':
    