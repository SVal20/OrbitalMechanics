import numpy as np

mu_sun = 1.327*10**11
mu_mars = 42830

Rearth = 149.6 *10**6
Rmars = 227.9*10**6

rmars = 3396

#a) the minimum delta-v required to place a spacecraft in orbit with a period of 7 h

Vinf = np.sqrt(mu_sun/Rmars) * (1-(np.sqrt(2*Rearth/(Rearth+Rmars))))
print(Vinf)
T = 7*3600

a = ((T*np.sqrt(mu_mars))/(2*np.pi))**(2/3)
print(a)
e = ((2*mu_mars)/(a*Vinf**2))-1
print(e)
dVmin = Vinf*np.sqrt((1-e)/2)
print(dVmin)

#b) The periapsis radius
