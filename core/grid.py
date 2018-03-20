#!/usr/bin/env python
import os
os.environ['GTK_MODULES'] = os.environ.get('GTK_MODULES', '').replace('pantheon-filechooser-module:', '')
# remove warning corresponding to missing GTK_MODULES

import code 										# file containing the c++ code for weave.inlining
import numpy as np, matplotlib.pyplot as plt 		# basic modules needed for computation and ploting
from matplotlib.patches import Circle, Rectangle 	# display current/voltage contacs
import time 										# measure computation time and give estimate
from datetime import datetime 						# convert datetimes in readable format
try:
    slow = False
    from scipy import weave
    from scipy.weave import converters
except ImportError:
    try:
        import weave
        from weave import converters
        print "Anaconda version of 4.2 or above found!\nUsing seperate Weave-module, may cause errors."
    except ImportError:
        print "No Weave module detected. Only slow calculations possible. To use fast methods install weave (via anaconda, ...)."
        slow = True  

# some constants for easy usage
ELECTRON_CHARGE = 1.6022e-19                   		# elementary charge of an electron in C
ELECTRON_MASS = 9.109e-31
ELECTRON_DENSITY_COPPER = lambda B: 8.5e22
ELECTRON_DENSITY_TaAs = lambda B: 3.5e19       		# zero field density for electrons in TaAs at room temperature in (cm)^-3 taken from Arnolds paper
ELECTRON_RELAX_TIME_COPPER = 2.5e-14
ELECTRON_MOBIL_TaAs = lambda T: 1.2e4          		# electron mobility calculated from dHvA oscillations in (cm^2/Vs) taken from Arnolds paper

# locations and shortened form (2d ignores 3rd coordinate)
TOP_LEFT_MIDDLE = [0.0, 0.5, 1.0]
TOP_RIGHT_MIDDLE = [1.0, 0.5, 1.0]
TOP_FRONT_MIDDLE = [0.5, 0.0, 1.0]
TOP_BACK_MIDDLE = [0.5, 1.0, 1.0]
LEFT = TOP_LEFT_MIDDLE
RIGHT = TOP_RIGHT_MIDDLE
FRONT = TOP_FRONT_MIDDLE
BACK = TOP_BACK_MIDDLE


class Grid:
    """A simple grid class that stores the details and computes the solution of the generalized Poisson equation with no-flux Neumann boundary conditions.
    """
    
    def __init__(self, size=[1.0, 1.0], ny=100, B=0., theta=0., T=300., mulGrid=1):
        if len(size) > 3 or len(size) < 2:
            raise UserWarning("Wrong dimensions!")
        size = sorted(size, reverse=True)
        if len(size) > 1:
            self.lx, self.ly = size[0]*1e-0, size[1]*1e-0
            self.dim = 2
            self.dx = self.ly*1e-6/(ny-1)
            self.nx = int(np.ceil(ny*self.lx/self.ly))
            self.ny = ny
            self.nz = 1
            self.lz = 0.
        if len(size) == 3:
            if slow:
                raise UserWarning("No 3d calculations with slow method possible!")
            self.lz = size[2]*1e-0
            self.dim = 3
            self.nz = int(np.ceil(ny*self.lz/self.ly))
        self.u = np.zeros((self.nx, self.ny, self.nz), dtype=np.float64)
        self.i = np.zeros((self.nx, self.ny, self.nz), dtype=np.float64)

        t = np.cos(np.pi/self.nx) + np.cos(np.pi/self.ny) # TODO check if really optimal solution
        if self.dim == 3: t += np.cos(np.pi(self.nz))
        self.wopt = (8-np.sqrt(64-16*t*t))/t**2

        self.B = B
        self.theta = theta
        self.T = T
        self.mulGrid = mulGrid
        self.currents = []
        self.probes = {}
        self.batches = {}
        self.outputs = []
        
        print "Initiate {}d grid with {}x{}".format(self.dim, self.nx, self.ny) + ("x{}".format(self.nz) if self.dim==3 else "") + " grid points"
        
    def getConductivity(self):
        if self.sigtype == 0: # returns the homogeneous sig_tensor as sig_xx, sig_xy, sig_yy, sig_zz
            if "mu" in self.rhoParams:
                mu = self.rhoParams["mu"](self.T)*1e-4
            else:
                mu = self.rhoParams["q"]*ELECTRON_CHARGE*self.rhoParams["tau"]/(self.rhoParams["m"]*ELECTRON_MASS)
            sig0 = self.rhoParams["n"](self.B)*1e6*self.rhoParams["q"]*mu
            div = 1+(mu*self.B*np.cos(self.theta))**2
            return (sig0/div, sig0*mu*self.B/div, sig0/(div*self.rhoParams["a"]), sig0)
        elif self.sigtype == 1: # returns the homogeneous sig_tensor as sig_xx, sig_xy, sig_yy, sig_zz
            if "mu1" in self.rhoParams:
                mu1, mu2 = self.rhoParams["mu1"](self.T)*1e-4, self.rhoParams["mu2"](self.T)*1e-4
            else:
                mu1 = self.rhoParams["q1"]*ELECTRON_CHARGE*self.rhoParams["tau1"]/(self.rhoParams["m1"]*ELECTRON_MASS)
                mu2 = self.rhoParams["q2"]*ELECTRON_CHARGE*self.rhoParams["tau2"]/(self.rhoParams["m2"]*ELECTRON_MASS)
            r1, r2 = 1./(self.rhoParams["q1"]*ELECTRON_CHARGE*self.rhoParams["n1"](self.B)*1e6), 1./(self.rhoParams["q2"]*ELECTRON_CHARGE*self.rhoParams["n2"](self.B)*1e6)
            rho1, rho2 = r1/mu1, r2/mu2
            r = (r1*rho2**2 + r2*rho1**2 + r1*r2*(r1+r2)*(self.B*np.cos(self.theta))**2)/((rho1+rho2)**2 + (r1+r2)**2*(self.B*np.cos(self.theta))**2)
            rho = (rho1*rho2*(rho1+rho2) + (rho1*r2**2+rho2*r1**2)*(self.B*np.cos(self.theta))**2)/((rho1+rho2)**2 + (r1+r2)**2*(self.B*np.cos(self.theta))**2)
            return (rho/((self.B*np.cos(self.theta))**2*r**2+rho**2), self.B*r/((self.B*np.cos(self.theta))**2*r**2+rho**2), rho/((self.B*np.cos(self.theta))**2*r**2+rho**2)/self.rhoParams["a"], 1./rho)
        elif self.sigtype == 2: # returns the conductivity as a nx*ny*nz array and factor a and b
            rho = np.zeros((self.nx, self.ny, self.nz), dtype=np.float64)
            nSurface = int(np.ceil(self.ny*self.rhoParams["surface_perc"]))
            rho[:,:,:] = self.rhoParams["sigxx"]/self.rhoParams["gSurBulk"]
            rho[:nSurface,:,:] = self.rhoParams["sigxx"]
            rho[-nSurface:,:,:] = self.rhoParams["sigxx"]
            rho[:,:nSurface,:] = self.rhoParams["sigxx"]
            rho[:,-nSurface:,:] = self.rhoParams["sigxx"]
            rho[:,:,:nSurface] = self.rhoParams["sigxx"]
            rho[:,:,-nSurface:] = self.rhoParams["sigxx"]
            return (rho, self.rhoParams["a"], 1)
            # TODO add factor between sigxx and sigxy, for actual computation not relevant
        elif self.sigtype == 3: # return the conductivity as a nx*ny*nz array
            return (self.rhoParams["sxx"], self.rhoParams["sxy"], self.rhoParams["syy"], self.rhoParams["sxx"]*0 if self.dim==2 else self.rhoParams["szz"])
        elif self.sigtype == 4: # homogeneous
        	sxx, syy = self.rhoParams["sxx"], self.rhoParams["syy"]
        	sxy, szz = self.rhoParams.get("sxy", 0.), self.rhoParams.get("szz", 0.)
        	return (sxx, sxy, syy, szz)
        elif self.sigtype == 5: # weyl
        	pass

    def setConductivity(self, sigtype="single", rand=[None], params={"q": 1, "m": 1, "n": ELECTRON_DENSITY_COPPER, "tau": ELECTRON_RELAX_TIME_COPPER, "a": 1}):
        if sigtype == "single":
            if all(k in params for k in ("q", "m", "n", "tau")) or all(k in params for k in ("q", "n", "mu")):
                print "Use single carrier model"
                self.sigtype = 0    
            else:
                raise UserWarning("Wrong amount of parameters for single carrier transport")       
        elif sigtype == "two":
            if all(k in params for k in ("q1", "q2", "m1", "m2", "n1", "n2", "tau1", "tau2")) or all(k in params for k in ("q1", "q2", "n1", "n2", "mu1", "mu2")):
                print "Use two carrier model"
                self.sigtype = 1
            else:
                raise UserWarning("Wrong amount of parameters for two carrier transport")
        elif sigtype == "ti":
            if self.dim != 3:
                raise UserWarning("Settings not compatible with dimensionality of grid")
            if len(params) != 6:
                raise UserWarning("Wrong amount of parameters for topological insulator model")
            if params["type"] not in ["costum", "single", "two"]:
                raise UserWarning("Unknown type for topological insulator model")
            self.sigtype = 2
            print "Use topological insulator model"
        elif sigtype == "costum":
        	if self.dim == 2:
        		if (params["sxx"].shape != (self.nx+1, self.ny+1, self.nz) or
        		params["sxy"].shape != (self.nx+1, self.ny+1, self.nz) or
        		params["syy"].shape != (self.nx+1, self.ny+1, self.nz)):
        			raise UserWarning("Sigma tensor doesnt match grid size + 1")
        	elif self.dim == 3:
	        	if (params["sxx"].shape != (self.nx+1, self.ny+1, self.nz+1) or
	        		params["sxy"].shape != (self.nx+1, self.ny+1, self.nz+1) or
	        		params["syy"].shape != (self.nx+1, self.ny+1, self.nz+1) or
	        		params["szz"].shape != (self.nx+1, self.ny+1, self.nz+1)):
	        		raise UserWarning("Sigma doesnt match grid size + 1")
	        self.sigtype = 3
	        print "Use costum model"
        elif sigtype == "homogeneous":
        	if "sxx" not in params or "syy" not in params:
        		raise UserWarning("Sigxx or Sigyy missing!")
        	else:
        		if self.dim == 3 and "szz" not in params:
        			raise UserWarning("Sigzz missing!")
        	self.sigtype = 4
        	print "Use homogeneous model"
        elif sigtype == "weyl": #TODO
        	if len(params) < 4:
        		raise UserWarning("Not enough parameters for semiclassical weyl model")
        	self.sigtype = 5
        	print "Use semiclassical weyl model"
        else:
            raise UserWarning(sigtype + " not supported")
            
        if len(rand) == 2:
            print "gaussian randomness introduced to sigma with amp={} and w={}".format(rand[0],rand[1])
            #TODO changes type of computation and requires additional array
        
        if "a" not in params:
            params["a"] = 1.
        self.rhoParams = params
          
    def addCurrent(self, mode="single", amp=1., params=LEFT):
        if mode == "single":
            self.i[int((self.nx-1)*params[0]), int((self.ny-1)*params[1]), int((self.nz-1)*params[2])] += amp*1e-3/self.dx**2
            current = {"amp": amp, "mode": mode, "pos": params}
            self.currents += [current]
            print "Added point-like current " + ("source" if amp>0. else "sink") + " with {}mA".format(amp) # TODO check if nA 
        elif mode == "circle":
            x, y, z = np.linspace(0, 1, self.nx), np.linspace(0, 1, self.ny), np.linspace(0, 1, self.nz)
            x, y, z = np.meshgrid(x, y, z, indexing='ij')
            current = {"amp": amp, "mode": mode, "pos": params[:3], "radius": params[3]}

            r = np.sqrt(((x - params[0])*self.lx)**2 + ((y - params[1])*self.ly)**2 + ((z - params[2])*self.lz)**2)
            inside = r < params[3]
            self.i[inside] += amp*1e-3/(inside.sum()*self.dx**2) # can throw a divide by zero if the circle is smaller than one grid point
            self.currents += [current]
            print "Added circle-like current " + ("source" if amp>0. else "sink") + " with %.1fmA and radius %.1fum"%(amp, params[3])
        elif mode == "line": #only top or bottom
            x, y = np.linspace(0, self.lx, self.nx), np.linspace(0, self.ly, self.ny)
            x, y = np.meshgrid(x, y, indexing='ij')
            current = {"amp": amp, "mode": mode, "pos1": params[:3], "pos2": params[3:6], "width": params[6]}

            x1, x2 = params[0]*self.lx, params[3]*self.lx
            y1, y2 = params[1]*self.ly, params[4]*self.ly
            z1, z2 = params[2]*self.lz, params[5]*self.lz
            w = params[6]

            theta = np.arctan(np.inf*(y2-y1) if x1==x2 else ((y2-y1)/(x2-x1)))
            dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            rotMatrix = np.array([[np.cos(theta),  np.sin(theta)],
                                  [-np.sin(theta), np.cos(theta)]])
            xt, yt = x-x1, y-y1
            xt, yt = np.einsum('ji, mni -> jmn', rotMatrix, np.dstack([xt, yt]))
            inside = (xt <= dist) & (xt > 0) & (np.abs(yt) <= w/2)
            self.i[inside, int(z1)] = amp*1e-3/(inside.sum()*self.dx**2) # can throw a divide by zero if area is smaller than one grid point
            self.currents += [current]
            print "Added line-like current " + ("source" if amp>0. else "sink") + " with %.1fmA and width %.1fum"%(amp, params[6])
        else:
            raise UserWarning(type + " not supported for current probes")     
            
    def addProbes(self, type1="point", pos1=FRONT, type2="point", pos2=BACK, name="V_H"): # Default args to measure Hall Voltage
        if isinstance(pos1, list) or isinstance(pos2, list):
            pos1 = np.array(pos1)
            pos2 = np.array(pos2)
            
        self.probes[name] = [type1, pos1, type2, pos2]
        print "Added probes with id: " + name

    def getPot(self, type, pos):
        if type == "point":
            return self.u[int(pos[0]*(self.nx-1)), int(pos[1]*(self.ny-1)), int(pos[2]*(self.nz-1))]
        elif type == "line":
            x, y = np.linspace(0, self.lx, self.nx), np.linspace(0, self.ly, self.ny)
            x, y = np.meshgrid(x, y, indexing='ij')

            x1, x2 = pos[0]*self.lx, pos[3]*self.lx
            y1, y2 = pos[1]*self.ly, pos[4]*self.ly
            z1, z2 = pos[2]*self.lz, pos[5]*self.lz
            w = pos[6]

            theta = np.arctan(np.inf*(y2-y1) if x1==x2 else ((y2-y1)/(x2-x1)))
            dist = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            rotMatrix = np.array([[np.cos(theta),  np.sin(theta)],
                                  [-np.sin(theta), np.cos(theta)]])
            xt, yt = x-x1, y-y1
            xt, yt = np.einsum('ji, mni -> jmn', rotMatrix, np.dstack([xt, yt]))
            inside = (xt <= dist) & (xt > 0) & (np.abs(yt) <= w)
            return self.u[inside, int(z1)].sum()/inside.sum() # can throw a divide by zero if area is smaller than one grid point
        
    def addBatch(self, var="B", params={"start": -1., "end": 1., "n": 20}):
        """B, T, mu, ... require start, end and n whereas voltages require only mode and n """

        if var in ["B", "T", "theta"] and len(params) == 3: # TODO add more possible vars
            tmp = np.linspace(params["start"], params["end"], params["n"])
            self.batches[var] = tmp
            print "Added batch job for %s from %.1f to %.1f with %d intervals"%(var, tmp[0], tmp[-1], tmp.size)
        elif var in self.probes and len(params) == 2 and params["mode"] in ["left", "right", "symm"]:
            dist = self.probes[var][0] - self.probes[var][1]
            dx = np.linspace(0, dist, params["n"])
            if params["mode"] == "left":
                self.batches[var] = [dx, dx*0.]
            elif params["mode"] == "right":
                self.batches[var] = [dx*0., -dx]
            else:
                self.batches[var] == [dx/2., -dx/2.]
            print "Added batch job for %s"%(var)
        else: # TODO if current or else 
            pass      
      
    def plotContour2d(self, n=20, side=0, vector=True, show=True, save=False, var="V_H"):
        """Convert the numeric representation of the selected surface to slices of the array"""
        if side < 4:
            r1 = slice(self.nx)
        elif side == 4:
            r1 = 0
        else:
            r1 = self.nx-1
            
        if side in [0, 2, 4, 5]:
            r2 = slice(self.ny)
        elif side == 1:
            r2 = 0
        else:
            r2 = self.ny-1
            
        if side in [1, 3, 4, 5]:
            r3 = slice(self.nz)
        elif side == 0:
            r3 = 0
        else:
            r3 = self.nz-1
        
        # calculate the gradient of the potential
        y, x = np.mgrid[0:self.ly:self.ny*1j, 0:self.lx:self.nx*1j]
        dy, dx = np.gradient(self.u[r1,r2,r3].T, np.diff(y[:2, 0]*1e-6), np.diff(x[0, :2]*1e-6))
        tmp = self.getConductivity()
        sigxx, sigxy, sigyy = tmp[0], tmp[1], tmp[2]
        
        if self.sigtype < 4:
            self.u[(sigxx[1:,1:,0] == 0.)] = np.nan
        # j = sig*E, convert from electric field to current
        #ddx = dx*sigxx-dy*sigxy
        #ddy = dx*sigxy+dy*sigyy
        #dx = ddx
        #dy = ddy
        
        u = self.u.copy()
        #fig, ax = plt.figure(figsize=(32,18)), plt.subplot(gs[0])#, plt.subplot(gs[1])
        min, max = np.nanmin(u[r1,r2,r3]), np.nanmax(u[r1,r2,r3])
        print min, max
        if np.abs(min) > np.abs(max):
            min = np.abs(max)*min/np.abs(min)
        else:
            max = np.abs(min)*max/np.abs(max)
        V = np.linspace(min,max,n)
        fig, ax = plt.subplots(figsize=(32,18))
        cont = ax.contourf(x, y, u[r1,r2,r3].T, V, cmap='gist_earth', alpha=0.599) # gist_earth, a=0.6
        ax.set_yticklabels([]) # 70 -> 40
        ax.set_xticklabels([])

        #dz = np.sqrt(dx**2+dy**2)*np.abs(self.u[:,:,0].T)
        #print dz.shape, self.u.shape, dx.shape
        #cont = ax.pcolormesh(x, y, dz, shading='gouraud', cmap='RdBu', vmin=0., vmax=dz.max()/2000)
        #cont = ax.imshow(dz)
        #fig.colorbar(cont)
        #print dz.max(), dz.min()


        for c in cont.collections:
            c.set_edgecolor("face")
        if vector: # add the electric field as streamplot
            c = np.max(np.hypot(dx, dy))/12 # define max line strength
            ax.streamplot(x, y, dx, dy, linewidth=np.hypot(dx, dy)/c*0+2, color=u[r1,r2,r3].T, density=0.4, cmap='gist_earth')
        
        # Add filled circles with text for the input currents 
        curr_colors = {True: '#aa0000', False: '#0000aa'}
        curr_labels = {True: 'I+', False: 'I-'}
        curr_align  = {True: 'right', False: 'left'}
        for current in self.currents:
            if current["mode"] in ["single", "circle"]:
                pos = current["pos"][:2]
                ax.add_patch(Circle(current["pos"]*np.array([self.lx,self.ly, 0]), 0.05*self.lx/10*self.ly, color=curr_colors[current["amp"]>0]))
                ax.annotate(curr_labels[current["amp"]>0], xy=(0,0), xycoords='axes fraction', xytext=pos+(np.array([0.5,0.5])-pos)/np.array([20,20]), fontsize=70, fontweight='bold', ha=curr_align[pos[0]>0.5])
            elif current["mode"] == "line":
                pos, width, vec = current["pos1"]*np.array([self.lx, self.ly, 0]), current["width"]/self.lx, (np.array(current["pos1"])-np.array(current["pos2"]))*np.array([self.lx, self.ly, 0])
                height, angle = np.linalg.norm(vec), np.arccos(np.clip(np.dot(vec/np.linalg.norm(vec), np.array([1,0,0])), -1.0, 1.0))
                anno = pos[:2]+np.array([np.cos(angle), np.sin(angle)])*height/2
                ax.add_patch(Rectangle(pos, width*2, height, angle/180, color=curr_colors[current["amp"]>0]))
                ax.annotate(curr_labels[current["amp"]>0], xy=(0,0), xycoords='data', xytext=anno+(np.array([0.5,0.5])-anno)/40, fontsize=40, fontweight='bold', ha='right')
        plt.savefig("homogen.png", frameon=False, dpi=100, bbox_inches='tight')
        return
        # Add filled circles with text for the voltage probes
        for key in self.probes:
            type1, pos1, type2, pos2 = self.probes[key] # TODO make for loop here
            if type1 == "point":
                ax.add_artist(Circle(pos1*np.array([self.lx,self.ly, 0]), 0.05*self.lx/10*self.ly, color='#455A64', alpha=0.7))
                ax.annotate(key[:1], xy=(0,0), xycoords='axes fraction', xytext=pos1[:2]+(np.array([0.65,0.5])-pos1[:2])/np.array([25,8]), fontsize=30, fontweight='bold', ha=curr_align[pos1[0]>0.5])
            elif type1 == "line":
                pos, width, vec = pos1[:3]*np.array([self.lx, self.ly, 0]), pos1[6]/self.lx, (pos1[:3]-pos1[3:6])*np.array([self.lx, self.ly, 0])
                height, angle = np.linalg.norm(vec), np.arccos(np.clip(np.dot(vec/np.linalg.norm(vec), np.array([1,0,0])), -1.0, 1.0))
                anno = (pos1[:2]+pos1[3:5])/2
                ax.add_patch(Rectangle(pos-width/2, width, height, angle/180, color='#455A64', alpha=0.3))
                ax.annotate(key[:1], xy=(0,0), xycoords='axes fraction', xytext=anno+(np.array([0.5,0.5])-anno)/np.array([40, 40])+np.array([1*width/self.lx,0])*0, fontsize=30, fontweight='bold', ha=curr_align[anno[0]>0.5])
            if type2 == "point":
                ax.add_artist(Circle(pos2*np.array([self.lx,self.ly, 0]), 0.05*self.lx/10*self.ly, color='#455A64', alpha=0.7))
                ax.annotate(key, xy=(0,0), xycoords='axes fraction', xytext=pos2[:2]+(np.array([0.65,0.5])-pos2[:2])/np.array([25,8]), fontsize=30, fontweight='bold', ha=curr_align[pos2[0]>0.5])
            elif type2 == "line":
                pos, width, vec = pos2[:3]*np.array([self.lx, self.ly, 0]), pos2[6]/self.lx, (pos2[:3]-pos2[3:6])*np.array([self.lx, self.ly, 0])
                height, angle = np.linalg.norm(vec), np.arccos(np.clip(np.dot(vec/np.linalg.norm(vec), np.array([1,0,0])), -1.0, 1.0))
                anno = (pos2[:2]+pos2[3:5]*np.array([1,0]))/2
                ax.add_patch(Rectangle(pos-width/2, width, height, angle/180, color='#455A64', alpha=0.4))
                ax.annotate(key, xy=(0,0), xycoords='axes fraction', xytext=anno+(np.array([0.5,0.5])-anno)/np.array([40, 40])-np.array([0.1*width/self.lx,0]), fontsize=30, fontweight='bold', ha=curr_align[anno[0]>0.5])
        
        var_table = {"B": self.B, "T": self.T, "theta": self.theta}
        #ax.set(aspect=1, title=r'Streamplot of electric field with contours of equipotential lines of $\phi$'+r' for ${}={:.2f}$'.format(var, var_table[var]))
        ax.set_title(r'Electric potential $\phi$'+r' for $A=\sigma_{xx}/\sigma_{yy}={%.1f}$ with a %dx%d grid'%(self.rhoParams["a"], self.nx, self.ny), fontsize=48, y=1.02)
        ax.set_xlabel('L (mm)', fontsize=40)
        ax.set_ylabel('W (mm)', fontsize=40)
        for label in ax.xaxis.get_majorticklabels():
            label.set_fontsize(32)
        for label in ax.yaxis.get_majorticklabels():
            label.set_fontsize(32)
        ax.yaxis.get_majorticklabels()[1].set_visible(False)
        """plt.tick_params(
            axis='x',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            bottom='off',      # ticks along the bottom edge are off
            top='off',         # ticks along the top edge are off
            labelbottom='off') # labels along the bottom edge are off
        plt.tick_params(
            axis='y',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            left='off',      # ticks along the bottom edge are off
            right='off',         # ticks along the top edge are off
            labelleft='off') # labels along the bottom edge are off
        """
        #cbar = fig.colorbar(cont, orientation='horizontal' )
        #cbar.ax.invert_xaxis() 
        #ax2.plot(np.linspace(0, self.lx, self.nx), u[:,self.ny/2, 0])
        ax.set_aspect(1)
        if save:            
            plt.savefig("{:%m%d}_2dcontour_{}_{:.2f}.png".format(datetime.now(), var, var_table[var]), frameon=False, dpi=200, bbox_inches='tight')
        if show:
            plt.show()
        else:
            plt.clf()

    def plotLinInt(self):
        fig, ax = plt.subplots(figsize=(12,18))
        x, midSection = np.linspace(0,4,self.nx), self.u[:,self.ny/2,0]
        linInt = lambda l: (midSection[self.nx/2+10]-midSection[self.nx/2-10])/(x[self.nx/2+10]-x[self.nx/2-10])*(l-x[self.nx/2]) + midSection[self.nx/2]
        ax.plot(x, midSection*2e10, color="b", linewidth=3, label="Middle section")
        ax.plot(x, linInt(x)*2e10, color="r", linewidth=4, alpha=0.5, label="Linear Interpolation")
        ax.set_xlabel('Length (mm)', fontsize=27)
        ax.set_ylabel(r"Potential $\phi$ (mV)", fontsize=33)
        for label in ax.xaxis.get_majorticklabels():
            label.set_fontsize(20)
        for label in ax.yaxis.get_majorticklabels():
            label.set_fontsize(24)
        
        ax.yaxis.get_majorticklabels()[0].set_visible(False)
        ax.grid()
        ax.legend(loc=1, fontsize=29)
        plt.savefig("middle.png", frameon=False, dpi=200, bbox_inches='tight')
        plt.show()  
             
    def plotContour3d(self, n=10, vector=True, show=True, save=False): # TODO
        pass       
       
    def plotVar(self, x, y, var="V_H", label=[], show=True, save=False):
        """Either makes a simple plot of one variable or if multidim. saves to disk as array"""
        if len(y.shape) == 1:
            fig, ax = plt.subplots()
            ax.set_title(var)
            ax.set_xlabel(label[0])
            ax.plot(x,y)
            if save:
                plt.savefig("{:%m%d}_plot_{}_{}.png".format(datetime.now()), var, label[0])
            if show:
                plt.show()
            else:
                plt.clf()
        else:
            if save:
                np.save(np.array([label, x, y]), "{:%m%d}_{}".format(datetime.now(), var))        
            
    def addOutput(self, mode="contour2dV", params={"n": 50, "side": 0}, show=True, save=False):
        """Adds the desired outputs to the simulation"""
            
        if mode == "contour2d" or mode == "contour2dV":
            self.outputs += [[mode, params, show, save]]
            print "Added %s output"%(mode)
        elif mode == "contour3d" or mode == "contour3dV":
            if self.dim == 2:
                raise UserWarning("3D plot not available for 2d grid")
            self.outputs += [[mode, params, show, save]]
            print "Added %s output"%(mode)
        elif mode == "plot":
            if params["var"] in self.probes:
                self.outputs += [[mode, params, show, save]]
                print "Added plot of %s"%(params["var"])
        elif mode == "current": # TODO add current measurements
            pass            
            
    def slowJacobi(self, tol=1e-4, nmax=5000): # only 2d, homogeneous conductivity
        """Calculate the solution only via simple numpy loops, no optimization"""

        u, I = self.u, self.i
        sig = self.getConductivity()[0]
        sxx, syy = sig[0], sig[2]
        nx, ny, dx2 = self.nx, self.ny, self.dx*self.dx
        
        for t in range(nmax):
            tmp = u.copy()

            i,j,k = 0, 0, 0 # left front bottom corner neumann boundary conditions with 1st order approximation
            u[i,j,k] = (sxx*u[i+1,j,k] + syy*u[i,j+1,k] + I[i,j,k]*dx2)/(sxx + syy)

            i, k = 0, 0 # left bottom edge
            u[i,1:-1,k] = (sxx*u[i+1,1:-1,k] + syy*(u[i,2:,k] + u[i,:-2,k]) + I[i,1:-1,k]*dx2)/(sxx + 2*syy)

            i, j, k = 0, ny-1, 0 # left back bottom corner
            u[i,j,k] = (sxx*u[i+1,j,k] + syy*u[i,j-1,k] + I[i,j,k]*dx2)/(sxx + syy)

            j, k = 0, 0 # front bottom edge
            u[1:-1,j,k] = (sxx*(u[2:,j,k] + u[:-2,j,k]) + syy*u[1:-1,j+1,k] + I[1:-1,j,k]*dx2)/(2*sxx + syy)

            k = 0 # bottom face
            u[1:-1,1:-1,k] = (sxx*(u[2:,1:-1,k] + u[:-2,1:-1,k]) + syy*(u[1:-1,2:,k] + u[1:-1,:-2,k]) + I[1:-1,1:-1,k]*dx2)/(2*sxx + 2*syy)

            j, k = ny-1, 0 # back bottom edge
            u[1:-1,j,k] = (sxx*(u[2:,j,k] + u[:-2,j,k]) + syy*u[1:-1,j-1,k] + I[1:-1,j,k]*dx2)/(2*sxx + syy)

            i, j, k = nx-1, 0, 0 # bottom front right corner
            u[i,j,k] = (syy*u[i,j+1,k] + sxx*u[i-1,j,k] + I[i,j,k]*dx2)/(syy + sxx)

            i, k = nx-1, 0 # bottom right edge
            u[i,1:-1,k] = (syy*(u[i,2:,k] + syy*u[i,:-2,k]) + sxx*u[i-1,1:-1,k] + I[i,1:-1,k]*dx2)/(2*syy + sxx)

            i, j, k = nx-1, ny-1, 0 # bottom back right corner
            u[i,j,k] = (sxx*u[i-1,j,k] + syy*u[i,j-1,k] + I[i,j,k]*dx2)/(sxx + syy)         

            sums = (tmp**2).sum()
            err = ((u-tmp)**2).sum()
            if sums != 0.:
                err = np.sqrt(err/sums)
                if err < tol:
                	return t
        return err      
         
    def fastSOR(self, w=1.2, tol=1e-4, nmax=10000): # homogeneous conductivity
        """Uses inline c++ code to increase the performance around 100x fold, add pointwise rel. residue and compare with tol"""

        u, I = self.u, self.i
        tmp = self.getConductivity()
        sxx, sxy, syy, szz = tmp[0], tmp[1], tmp[2], tmp[3]
        nx, ny, nz, _dx2 = self.nx, self.ny, self.nz, self.dx*self.dx

        if self.sigtype in [0, 1]:
        	_code = code.sor_2d_hom_xy if self.dim == 2 else code.sor_3d_inh_xy
        elif self.sigtype in [2, 3]:
        	_code = code.sor_2d_inh_xy if self.dim == 2 else code.sor_3d_inh_xy
        elif self.sigtype in [4, 5]:
            _code = code.sor_2d_hom if self.dim == 2 else code.sor_3d_hom

        # compiler keyword only needed on windows with MSVC installed
        t = weave.inline(_code, ['u', 'sxx', 'sxy', 'syy', 'szz', 'I', '_dx2', 'nx', 'ny', 'nz', 'nmax', 'w', 'tol'],
                         type_converters = converters.blitz) #,compiler = 'gcc')
        return t
                
    def run(self, w=None, tol=1e-6, nmax=10000, verbose=True, method="fast"):
        """Solves the generalized Poisson equation for a specific grid and tolerance and handles the output/batch jobs"""

        if w is None:
            w = float(self.wopt)
        if verbose:
            print "Solving poisson equation with relaxation parameter w=%.1f and a tolerance of %.1e"%(w, tol)

        if np.abs(self.i.sum()) > 1e-6: # if sum of j is bigger than a few uA
            raise UserWarning("Compatibility condition for flux conversation is not fulfilled %.1e"%(self.i.sum()))
        
        plot_tmp = {}
        var_table = {0: "B", 1: "B", 2: "T", 3: "theta"}
        for plot in self.outputs:
            mode, params, show, save = plot
            if mode == "plot":
                plot_tmp[params["var"]] = []
        
        def recursive(pos=0): # handles the recursive nature of batch jobs
            if "B" in self.batches and pos < 1:
                bRange = self.batches["B"]
                print "Vary magnetic field"
                for b in bRange:
                    self.B = b
                    print "Set B to %.2f"%(b)
                    recursive(1)
                
            elif "T" in self.batches and pos < 2:
                tRange = self.batches["T"]
                print "Vary temperature"
                for t in tRange:
                    self.T = t
                    print "Set T to %.2f"%(t)
                    recursive(2)
            elif "theta" in self.batches and pos < 3:
                thetaRange = self.batches["theta"]
                print "Vary angle of magnetic field"
                for theta in thetaRange:
                    self.theta = theta
                    print "Set theta to %.2f"%(theta)
                    recursive(3)
            elif pos < 4:
                #self.reset() # can't reuse old values due to the highly uneven nature of the potential
                t1 = time.clock()
                if method == "slow" or slow:
                	steps = self.slowJacobi(tol, nmax)
                else:
                	steps = self.fastSOR(w, tol, nmax)
                dt = time.clock()-t1
                    
                if steps < 1:
                    if verbose:
                        print "Couldn't reach tolerance after %f seconds with a rest error of %.2e"%(dt, steps)
                else:
                    if verbose:
                        print "Reached tolerance after %d steps and %f seconds"%(steps, dt)
                    for plot in self.outputs:
                        mode, params, show, save = plot
                        if mode == "plot":
                            if params["var"] in self.probes:
                                type1, pos1, type2, pos2 = self.probes[params["var"]]
                                potA = self.getPot(type1, pos1)
                                potB = self.getPot(type2, pos2)
                                plot_tmp[params["var"]] += [potA-potB]
                        elif mode == "contour2d": # at end because changes u
                            self.plotContour2d(vector=False, show=show, save=save, var=var_table[pos], **params)
                        elif mode == "contour2dV": # at end because changes u
                            self.plotContour2d(vector=True, show=show, save=save, var=var_table[pos], **params)
                        else: # TODO
                            pass
            else: # TODO probes and stuff
                pass
        
        recursive() # start going through all batches
        plot_var = np.array((1))
        plot_label = []
        for var in ["B", "T", "theta"]: # TODO add more ....
            if var in self.batches:
                if var == "B":
                    tmp = self.batches["B"]
                    plot_label += ["B"]
                elif var == "T":
                    tmp = self.batches["T"]
                    plot_label += ["T"]
                elif var == "theta":
                    tmp = self.batches["theta"]
                    plot_label += ["theta"]
                else: # TODO
                    break
                plot_var = np.outer(plot_var, tmp)
        for plot in self.outputs:
            mode, params, show, save = plot
            if mode == "plot":
                self.plotVar(np.squeeze(plot_var), np.array(plot_tmp[params["var"]]), params["var"], label=plot_label, show=show, save=save)
                
    def reset(self):
        """Sets the field back to zero"""
        self.u = np.zeros((self.nx, self.ny, self.nz))