from grid import *
import numpy as np

if __name__ == "__main__":
	# 2d example with homogeneous conductivity and no B field, just contour plot
    #grid = Grid(size=[4.46, 3.32], ny=150, B=0.) # Initiate grid
    grid = Grid(size=[4., 2.], ny=100, B=0.) # Initiate grid
    sxx = np.full((grid.nx+1, grid.ny+1, grid.nz), 1e6)
    sxy = np.full((grid.nx+1, grid.ny+1, grid.nz), 0)
    syy = np.full((grid.nx+1, grid.ny+1, grid.nz), 1e6)
    #sxx[120:,90:,0]=0.
    #syy[120:,90:,0]=0. 
    x, y = np.linspace(0, grid.lx, grid.nx), np.linspace(0, grid.ly, grid.ny)
    x, y = np.meshgrid(x, y, indexing='ij')
    inside = ((-(x-grid.lx)*0.75*grid.lx/grid.ly+0.7*grid.ly) < y)
    sxx2 = sxx[1:,1:,0]
    syy2 = syy[1:,1:,0]
    sxx2[inside] = 0.
    syy2[inside] = 0.
    params = {"sxx": sxx, "sxy": sxy, "syy": syy}
    grid.setConductivity(sigtype="homogeneous", params={"sxx": 1e+6*1.0, "syy": 1e+6, "a": 1}) # use simple homogeneous conductivity
    #grid.setConductivity(sigtype="costum", params=params)
    grid.addCurrent(amp=1, params=[0.0, 0.5, 0.0]) # add a current source on the left
    grid.addCurrent(amp=-1, params=[1.0, 0.5, 0.0]) # add a current source on the left
    #grid.addCurrent(amp=-1, mode="line", params=[0.99, -0.01, 0.0, 0.99, 0.67, 0.0, 0.1]) # add a current sink on the right
    grid.addProbes(pos1=LEFT+np.array([0.1,0.0,0.0]), pos2=RIGHT-np.array([0.1,-0.2,0.0]), name="V1") # add resisitivy probes V_R
    grid.addProbes(pos1=LEFT+np.array([0.1,0.0,0.0]), pos2=RIGHT-np.array([0.1,0.0,0.0]), name="V2") # add resisitivy probes V_R
    grid.addProbes(pos1=LEFT+np.array([0.1,0.0,0.0]), pos2=RIGHT-np.array([0.1,0.2,0.0]), name="V3") # add resisitivy probes V_R
    grid.addProbes(pos1=LEFT+np.array([0.1,0.0,0.0]), type2="line", pos2=np.array([0.9,0.0,0.0,0.9,1.0,0.0,0.1]), name="V4") # add resisitivy probes V_R
    #plt.imshow(grid.i[:,:,0].T,interpolation="none")
    #plt.show()
    grid.addOutput(mode="contour2d", show=True, save=False, params={"n": 24}) # add contour + vector stream plot
    grid.run(w=1.9, tol=1e-7, nmax=50000) # run the simulation with your settings

    b = grid.getPot("point", LEFT+np.array([0.1,0.0,0.0])) - grid.getPot("point", RIGHT-np.array([0.1,-0.2,0.0]))
    c = grid.getPot("point", LEFT+np.array([0.1,0.0,0.0])) - grid.getPot("point", RIGHT-np.array([0.1,0,0.0]))
    d = grid.getPot("point", LEFT+np.array([0.1,0.0,0.0])) - grid.getPot("point", RIGHT-np.array([0.1,0.2,0.0]))
    e = grid.getPot("point", LEFT+np.array([0.1,0.0,0.0])) - grid.getPot("line", np.array([0.9,0.0,0.0,0.9,1.0,0.0,0.1]))
    print "%.2f\t%.3e\t%.3e\t%.3e\t%.3e"%(50.0,b,c,d,e)
    plt.show()
    #grid.plotLinInt()
    #print grid.u.max(), grid.u.min()


    # two carrier, vary B field from -3T to 3T, plot Hall Voltage and probe in the middle
    #grid = Grid(size=[1.0, 0.1], ny=50, B=0) # Initiate grid    
    #grid.setConductivity(sigtype="two", rand=[None], params={"q1": 1, "n1": ELECTRON_DENSITY_TaAs, "mu1": ELECTRON_MOBIL_TaAs,
    #                                                         "q2":-1, "n2": ELECTRON_DENSITY_COPPER, "mu2": lambda T: 1e-0}) # two type carriers
    #grid.addCurrent() # add a current source on the left
    #grid.addCurrent(amp=-1., params=RIGHT) # add a current sink on the right
    #grid.addProbes() # add default Hall voltage V_H probes
    #grid.addProbes(pos1=LEFT+np.array([0.1,0.0,0.0]), pos2=RIGHT-np.array([0.1,0.0,0.0]), name="V_R") # add resisitivy probes V_R
    #grid.addBatch(var="B", params={"start": -3.0, "end": 3.0, "n": 20}) # allow the simulation to vary the magnetic field
    #grid.addOutput(mode="plot", params={"var": "V_H"}) # add plot of V_H(B)
    #grid.addOutput(mode="plot", params={"var": "V_R"}) # add plot of V_H(B)
    #grid.run(tol=1e-4, nmax=10000) # run the simulation with your settings