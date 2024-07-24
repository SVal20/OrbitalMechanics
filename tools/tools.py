import numpy as np
import matplotlib.pyplot as plt 
import math as m

import planetary_data as pd
plt.style.use('dark_background')

deg2rad = np.pi/180.0
rad2deg = 180/np.pi

def normed(v):
    return np.array(v)/np.linalg.norm(v)
    
def plot_n_orbits(rs,labels,cb=pd.earth,show_plot= False, save_plot = False, title="3D Plot"):
        
        fig = plt.figure(figsize=(12,6))
        ax=fig.add_subplot(111,projection='3d')

        n = 0
        for r in rs:
            ax.plot(r[:,0],r[:,1],r[:,2],label=labels[n])
            ax.plot([r[0,0]],[r[0,1]],[r[0,2]])
            n+=1

        _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x=cb['radius']*np.cos(_u)*np.sin(_v)
        _y=cb['radius']*np.sin(_u)*np.sin(_v)
        _z=cb['radius']*np.cos(_v)
        #ax.plot_surface(_x,_y,_z,cmap='Blues')
        ax.plot_wireframe(_x, _y, _z, color = "w", linewidth = 0.7)

        l = cb['radius']*2
        x,y,z=[[0,0,0],[0,0,0],[0,0,0]]
        u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
        ax.quiver(x,y,z,u,v,w,color='k')

        max_val=np.max(np.abs(rs))
    
        ax.set_xlim([-max_val,max_val])
        ax.set_ylim([-max_val,max_val])
        ax.set_zlim([-max_val,max_val])

        ax.set_xlabel(['X (km)'])
        ax.set_ylabel(['Y (km)'])
        ax.set_zlabel(['Z (km)'])

        ax.set_title(title)

        plt.legend()

        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png',dpi=300)

#oes = Orbital elements
def OrbElems2rv(oes,a,deg,mu=pd.earth['mu']):
    
    if a:
        a, e, i, tanom, argOfPrge, rAscnANode = oes

        if deg:
            i*=deg2rad
            tanom *=deg2rad
            argOfPrge *=deg2rad
            rAscnANode *=deg2rad

        E = ecc_anomaly([tanom,e],'tae')

        r_norm = a*(1-e**2)/(1+e*np.cos(tanom))

        rperifocal = r_norm*np.array([m.cos(tanom),m.sin(tanom),0])
        vperifocal = m.sqrt(mu*a)/r_norm*np.array([-m.sin(E),m.cos(E)*m.sqrt(1-e**2),0])

        periF2ECI = np.transpose(ECI2PeriF(rAscnANode,i,argOfPrge))

        r = np.dot(periF2ECI, rperifocal)
        v = np.dot(periF2ECI, vperifocal)
        return r,v


    else:
        h, e, i, tanom, argOfPrge, rAscnANode = oes

        if deg:
            i*=deg2rad
            tanom *=deg2rad
            argOfPrge *=deg2rad
            rAscnANode *=deg2rad

        r = (h**2/mu) * (1/(1+e*np.cos(tanom)))
        rperifocal = r*np.array([np.cos(tanom),np.sin(tanom),0])

        vperifocal = (mu/h) * np.array([-np.sin(tanom),e+np.cos(tanom),0])
        

        PeriF2GEC = np.transpose(GEC2PeriF(rAscnANode,i,argOfPrge))
        
        rGEC = np.matmul(PeriF2GEC, rperifocal)
        vGEC = np.matmul(PeriF2GEC, vperifocal)
    
        return rGEC,vGEC

#GEC: Geocentric equatorial coordinates
def GEC2PeriF(rAscNde,i,ArgOfP):

    matrix_RAAN = [
                    [np.cos(rAscNde), np.sin(rAscNde), 0],
                    [-np.sin(rAscNde), np.cos(rAscNde), 0],
                    [0, 0, 1]
                  ]

    matrix_i = [
                [1, 0, 0],
                [0, np.cos(i), np.sin(i)],
                [0, -np.sin(i), np.cos(i)]
               ]
    
    matrix_AOP = [
                  [np.cos(ArgOfP),np.sin(ArgOfP),0],
                  [-np.sin(ArgOfP),np.cos(ArgOfP),0],
                  [0, 0, 1]
                 ]
    
    matrix_RAAN_np = np.array(matrix_RAAN)
    matrix_i_np = np.array(matrix_i)
    matrix_AOP_np = np.array(matrix_AOP)

    TempMatrix = np.matmul(matrix_RAAN_np, matrix_i_np)
    CosineMatrix = np.matmul(TempMatrix, matrix_AOP_np)
    return CosineMatrix
#ECI: Earth centered intertial
def ECI2PeriF(rAscNde,i,ArgOfP):
    row0 = [-m.sin(rAscNde)*m.cos(i)*m.sin(ArgOfP)+m.cos(rAscNde)*m.cos(ArgOfP), m.cos(rAscNde)*m.cos(i)*m.sin(ArgOfP)+m.sin(rAscNde)*m.cos(ArgOfP),m.sin(i)*m.sin(ArgOfP)]
    row1 = [-m.sin(rAscNde)*m.cos(i)*m.cos(ArgOfP)-m.cos(rAscNde)*m.sin(ArgOfP), m.cos(rAscNde)*m.cos(i)*m.cos(ArgOfP)-m.sin(rAscNde)*m.sin(ArgOfP),m.sin(i)*m.cos(ArgOfP)]
    row2 = [m.sin(rAscNde)*m.sin(i),-m.cos(rAscNde)*m.sin(i),m.cos(i)]
    return np.array([row0,row1,row2])  
  
def ecc_anomaly(arr,method,tol=1e-8):


    if method == 'newton':
        Me,e=arr
        if Me < np.pi/2.0: E0=Me+e/2.0
        else: E0=Me-e
        for n in range(200):
            ratio = (E0-e*np.sin(E0)-Me)/(1-e*np.cos(E0))
            if abs(ratio)<tol:
                if n==0: return E0
                else: return E1
            else:
                E1 = E0-ratio
                E0 = E1
        return False
    elif method == 'tae':
        ta,e=arr
        return 2*m.atan(m.sqrt((1-e)/(1+e))*m.tan(ta/2.0))
    else:
        print("INVALID")

def TLES2OES(tles_file,mu=pd.earth['mu']):
    with open(tles_file,'r') as f:
        lines = f.readlines()
    line0=lines[0].strip()
    line1=lines[1].strip().split()
    line2=lines[2].strip().split()
    
    epoch = line1[3]

    i = float(line2[2])*deg2rad
    
    #Right ascension of the ascending node
    RAAN = float(line2[3])*deg2rad
    e_string = line2[4]
    e=float('0.'+e_string)
    #Argument of perigee
    AoP = float(line2[5])*deg2rad
    #Mean anomaly
    Me = float(line2[6])*deg2rad
    #Mean motion
    n = float(line2[7]) # rev/day
    #Period
    T = 1/n*24*3600
    #Semi-major axis
    a = (T**2*mu/4.0/np.pi**2)**(1/3.0)
    #Eccentric Anomaly
    E = ecc_anomaly([Me,e],'newton')
    #True anomaly
    tanom = 2*np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2.0))

    return a, e, i, tanom, AoP, RAAN

def rv2OEs(state_vector,mu = pd.earth['mu'], deg=False, print = False):
    rs,vs = state_vector
    r = m.sqrt(np.dot(rs,rs))
    v = m.sqrt(np.dot(vs,vs))

    radial_vel = np.dot(rs,vs/r)
    h = np.cross(rs,vs)
    ang_momentum = m.sqrt(np.dot(h,h))

    hz = h[2]
    i = m.cosh(hz/ang_momentum)
    Node_line = np.cross([0,0,1],h) 

    N = m.sqrt(np.dot(Node_line,Node_line))

    if Node_line[1]>0:
        RAAN = m.cosh(Node_line[0]/N)
    else:
        RAAN = 360 - m.cosh(Node_line[0]/N)
                
    eccentricity = (1/mu)*((v**2-(mu/r))*rs-(r*v*vs))
    e = m.sqrt(np.dot(eccentricity,eccentricity))

    ez = eccentricity[2] 

    if  ez > 0:
        AOP = m.cosh(Node_line*eccentricity/N*e)
    else:
        AOP = 360 - m.cosh(Node_line*eccentricity/N*e)

    if radial_vel >= 0:
        TrueAnom = m.cosh(np.dot(eccentricity/e, rs/r))
    else:
        TrueAnom = 360 - m.cosh(np.dot(eccentricity/e, rs/r))


    if print:
        print (h, i, RAAN, e, AOP, TrueAnom)

def rv2OEs2(r,v, mu = pd.earth['mu'], degrees = False, print = False):
    
    r_norm = np.linalg.norm(r)
    v_norm = np.linalg.norm(v)

    h = np.cross(r,v)
    h_norm = np.linalg.norm(h)

    i = m.acos(h[2]/h_norm)

    e = ((v_norm**2-mu/r_norm)*r-np.dot(r,v)*v)/mu
    e_norm = np.linalg.norm(e)

    N = np.cross([0,0,1],h)
    N_norm = np.linalg.norm(N)

    RAAN = m.acos(N[0]/N_norm)
    if N[1]<0: RAAN=2*np.pi-RAAN

    AOP = m.acos(np.dot(N,e)/N_norm/e_norm)
    if e[2]<0: AOP=2*np.pi-AOP

    ta = m.acos(np.dot(e,r)/e_norm/r_norm)
    if np.dot(r,v)<0: ta=2*np.pi-ta

    a = r_norm*(1+e_norm*m.cos(ta))/(1-e_norm**2)

    if print:
        print('a', a)
        print('e', e_norm)
        print('i', i*rad2deg)
        print('RAAN', RAAN*rad2deg)
        print('AOP', AOP*rad2deg)
        print('TA', ta*rad2deg)

    if degrees: return [a,e_norm, i*rad2deg, ta*rad2deg, AOP*rad2deg, RAAN*rad2deg]
    else: return [a,e_norm,i,ta,AOP,RAAN]

def calc_atmospheric_density(z):
    rhos, zs = find_rho_z(z)
    if rhos[0] == 0: return 0,0

    Hi = -(zs[1]-zs[0])/m.log(rhos[1]/rhos[0])

    return rhos[0]*m.exp(-(z-zs[0])/Hi)

def find_rho_z(z,zs = pd.earth['zs'],rhos=pd.earth['rhos']):
    if not 1.0<z<1000.0:
        return [[0,0,0],[0,0,0]]
    
    for n in range(len(rhos)-1):
        if zs[n]<z<zs[n+1]:
            return [[rhos[n],rhos[n+1]],[zs[n],zs[n+1]]]

    return [[0,0,0],[0,0,0]]