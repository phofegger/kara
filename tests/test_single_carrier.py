import sys
print(sys.path)

from ..core import grid

if __name__ == "__main__":
	# 2d example with homogeneous conductivity and no B field, just contour plot
    grid = Grid(size=[2., 1.], ny=100, B=0.) # Initiate grid
    grid.setConductivity(sigtype="homogeneous", params={"sxx": 1, "syy": 1}) # use homogeneous conductivity
    grid.addCurrent(amp=1, params=[0.0,0.5,0.0]) # add a current source on the left
    grid.addCurrent(amp=-1, params=[1.,0.5,0.0]) # add a current sink on the right
    grid.addOutput(mode="contour2dV", show=True, save=False, params={"n": 30}) # add contour + vector stream plot
    grid.run(w=1.8, tol=1e-4, nmax=10000) # run the simulation with your settings
    plt.show()