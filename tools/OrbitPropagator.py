import numpy as np
import matplotlib.pyplot as plt 
from scipy.integrate import ode

import planetary_data as pd
import tools as t

def null_perts():
    return {
        'J2': False,
        'aero': False,
        'thrust': 0,
        'isp': 0
    }

class OrbitPropagator:
    def __init__(self, state0, tspan, dt, a=True, coes=True, deg= True, cb=pd.earth, perts = null_perts(), mass0 = 0, sc = {}):
        if coes :
            self.r0,self.v0 = t.OrbElems2rv(state0,a,deg,cb['mu'])
        else:
            self.r0 = state0[:3]
            self.v0 = state0[3:]

        self.tspan = tspan
        self.dt = dt
        self.cb = cb
        self.mass0 = mass0

        self.n_steps = int(np.ceil(self.tspan/self.dt))+1

        self.ys = np.zeros((self.n_steps+1,7))
        self.ts = np.zeros((self.n_steps+1,1))
        self.alts = np.zeros((self.n_steps+1))
         
        self.y0 = self.r0.tolist()+self.v0.tolist()+[self.mass0]
        self.ys[0] = np.array(self.y0)
        self.alts[0] = np.linalg.norm(self.r0)-self.cb['radius']
        self.step = 0
        
        
        self.solver = ode(self.diffy_q)
        self.solver.set_integrator('lsoda')
        self.solver.set_initial_value(self.y0,0)

        self.perts = perts

        self.stop_conditions_dictionary = sc
        self.stop_conditions_map = { 'max_alt': self.check_max_alt, 'min_alt': self.check_min_alt}
        self.stop_conditions_functions = [self.check_deorbit]

        for key in self.stop_conditions_dictionary.keys():
            self.stop_conditions_functions.append(self.stop_conditions_map[key])

        self.propagate_orbit()
    
    def check_deorbit(self):
        if self.alts[self.step] < self.cb['deorbit_altitud']:
            print ("Spacecraft deorbited after %.1f seconds" %self.ts[self.step])
            return False
        else:
            return True
    
    def check_max_alt(self):
        if self.alts[self.step] > self.stop_conditions_dictionary['max_alt']:
            print("Spacecraft reached max altitude after %.1f seconds." %self.ts[self.step])
            return False
        else:
            return True

    def check_min_alt(self):
        if self.alts[self.step] < self.stop_conditions_dictionary['min_alt']:
            print("Spacecraft reached min altitude after %.1f seconds." %self.ts[self.step])
            return False
        else:
            return True

    def check_stop_conditions(self):
        for sc in self.stop_conditions_functions:
            if not sc():
                return False
        return True
            
                
    def propagate_orbit(self):
        while self.solver.successful() and self.step<self.n_steps and self.check_stop_conditions():
            self.solver.integrate(self.solver.t+self.dt)
            self.step+=1
            
            self.ts[self.step] = self.solver.t
            self.ys[self.step] = self.solver.y
            self.alts [self.step] = np.linalg.norm(self.solver.y[:3])-self.cb['radius']
            #print(self.alts[self.step])
        
        self.ts = self.ts[:self.step]
        self.rs = self.ys[:self.step,:3]
        self.vs = self.ys[:self.step,3:6]
        self.masses = self.ys[:self.step, -1]
        self.alts = self.alts[:self.step]
        
    def calculate_OEs(self,degrees=True):

        print("Calculating OEs...")

        self.OEs = np.zeros((self.n_steps,6))

        for n in range(self.step):
            self.OEs[n,:] = t.rv2OEs2(self.rs[n,:],self.vs[n,:],mu=self.cb['mu'],degrees=degrees)

    def plot_OEs(self, hours= False, days=False, show_plot = False, save_plot = False, title = "OEs", figsize=(16,8)):
        print("Plotting Orbital Elements...")

        fig,axs = plt.subplots(nrows=2,ncols=3,figsize=figsize)

        fig.suptitle(title, fontsize=20)

        if hours:
            ts=self.ts/3600.0
            xlabel = "Time elapsed (hours)"
        elif self.days:
            ts=self.ts/3600.0/24.0
            xlabel = "Time elapsed (days)"
        else:
            ts=self.ts
            xlabel = "Time elapsed (seconds)"

        axs[0,0].plot(ts, self.OEs[:self.step,3])
        axs[0,0].set_title("True anomaly vs time")
        axs[0,0].grid(True)
        axs[0,0].set_ylabel("Angle (deg)")
        axs[0,0].set_xlabel(xlabel)

        axs[1,0].plot(ts, self.OEs[:self.step,0])
        axs[1,0].set_title("Semi-major axis vs time")
        axs[1,0].grid(True)
        axs[1,0].set_ylabel("Semi-major axis (km)")
        axs[1,0].set_xlabel(xlabel)
        
        axs[0,1].plot(ts, self.OEs[:self.step,1])
        axs[0,1].set_title("Eccentricity vs time")
        axs[0,1].grid(True)
        axs[0,1].set_ylabel("Eccentricity")
        axs[0,1].set_xlabel(xlabel)

        axs[0,2].plot(ts, self.OEs[:self.step,4])
        axs[0,2].set_title("Argument of periapse vs time")
        axs[0,2].grid(True)
        axs[0,2].set_ylabel("Argument of periapse (deg)")
        axs[0,2].set_xlabel(xlabel)

        axs[1,1].plot(ts, self.OEs[:self.step,2])
        axs[1,1].set_title("Inclination vs time")
        axs[1,1].grid(True)
        axs[1,1].set_ylabel("Inclination (deg)")
        axs[1,1].set_xlabel(xlabel)

        axs[1,2].plot(ts, self.OEs[:self.step,5])
        axs[1,2].set_title("Right ascension of ascending node vs time")
        axs[1,2].grid(True)
        axs[1,2].set_ylabel("Right ascension of ascending node (deg)")
        axs[1,2].set_xlabel(xlabel)

        if show_plot:
            plt.show()

        if save_plot:
            plt.savefig(title+".png",dpi=300)

    def diffy_q(self,t_,y):
        rx,ry,rz,vx,vy,vz,mass = y
        r = np.array([rx,ry,rz])
        v = np.array([vx,vy,vz])

        norm_r = np.linalg.norm(r)

        a = -r*self.cb['mu']/norm_r**3

        if self.perts['J2']:
            z2 = r[2]**2
            r2 = norm_r**2
            tx = r[0]/norm_r*(5*z2/r2-1)
            ty = r[1]/norm_r*(5*z2/r2-1)
            tz = r[2]/norm_r*(5*z2/r2-3)

            a_j2 = 1.5*self.cb['J2']*self.cb['mu']*self.cb['radius']**2/norm_r**4*np.array([tx,ty,tz])
            a+= a_j2

        if self.perts['aero']:
            z = norm_r - self.cb['radius']
            rho = t.calc_atmospheric_density(z)
 
            v_rel = v-np.cross(self.cb['atm_rot_vector'],r)
            drag = v_rel*0.5*rho*np.linalg.norm(v_rel)*self.perts['Cd']*self.perts['A']/self.mass0

            a+= drag
        
        if self.perts['thrust']:
            a+= self.perts['thrust_direction'] * t.normed(v) * self.perts['thrust']/mass/1000.0

            dmdt = -self.perts['thrust']/self.perts['isp']/9.81

        return [vx,vy,vz,a[0],a[1],a[2],dmdt]
    
    def plot_3d(self,show_plot= False, save_plot = False, title="3D Plot"):
        fig = plt.figure(figsize=(12,6))
        ax=fig.add_subplot(111,projection='3d')

        ax.plot(self.rs[:,0],self.rs[:,1],self.rs[:,2],'w',label='Trajectory')
        ax.plot([self.rs[0,0]],[self.rs[0,1]],[self.rs[0,2]],'wo', label = "Initial Position")

        _u,_v=np.mgrid[0:2*np.pi:20j,0:np.pi:10j]
        _x=self.cb['radius']*np.cos(_u)*np.sin(_v)
        _y=self.cb['radius']*np.sin(_u)*np.sin(_v)
        _z=self.cb['radius']*np.cos(_v)
        #ax.plot_surface(_x,_y,_z,cmap='Blues')
        ax.plot_wireframe(_x, _y, _z, color = "k", linewidth = 0.7)

        '''max_val = np.max(np.abs(self.rs))
        ecliptic_x = np.linspace(-max_val, max_val, 100)
        ecliptic_y = np.linspace(-max_val, max_val, 100)
        ecliptic_x, ecliptic_y = np.meshgrid(ecliptic_x, ecliptic_y)
        ecliptic_z = np.zeros_like(ecliptic_x)
        ax.plot_surface(ecliptic_x, ecliptic_y, ecliptic_z, color='c', alpha=0.2, label='Ecliptic Plane')'''
        
        l = self.cb['radius']*2
        x,y,z=[[0,0,0],[0,0,0],[0,0,0]]
        u,v,w = [[l,0,0],[0,l,0],[0,0,l]]
        ax.quiver(x,y,z,u,v,w,color='k')

        max_val=np.max(np.abs(self.rs))
    
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

    def plot_alts(self,show_plot = False, save_plot = False, hours = False, days = False, title = 'Radial distance vs time', figsize = (16,8), altitude = True, dpi = 500):
        if hours:
            ts=self.ts/3600.0
            xlabel = "Time elapsed (hours)"
        elif self.days:
            ts=self.ts/3600.0/24.0
            xlabel = "Time elapsed (days)"
        else:
            ts=self.ts
            xlabel = "Time elapsed (seconds)"

        plt.figure(figsize=figsize)
        plt.plot(ts,self.alts,'w')
        plt.grid(True)
        plt.xlabel('Time (%s)' % xlabel)
        plt.ylabel('Altitude (km)')
        plt.title(title)
        if show_plot:
            plt.show()
        if save_plot:
            plt.savefig(title+'.png', dpi=dpi)