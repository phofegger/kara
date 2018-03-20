Grid class documentation
method __init__
	size as [lx, ly, lz] --	length of grid in um (y axis used as short axis, if no lz defaults to 2d)
	ny                   -- number of grid points along the y-axis
	B                    -- magnetic field strength in Tesla (homogeneous over whole sample)
    theta                -- enclosed angle of B to the z-axis in radians
	T                    -- temperature in Kelvin
	mulGrid 			 -- number of multigrid calculations, use only for complex computations (EXPERIMENTAL)
        
method setConductivity
	sigtype              -- select type of conductivity, possible values: 'single', 'two', 'ti', 'costum', h'omogeneous', 'weyl'
	rand                 -- add randomness to the conductivity to account for irregularites in the conductivity [amp, sig], (EXPERIMENTAL)
	params               -- single carrier (either q,m,n,tau,a or q,n,mu,a):
                            	q - charge in e units
                            	m - effective mass in terms of electron mass m_e
                        		n - carrier density in (cm)^-3
                        		tau - relaxation time in s
                        		a - anisotropy factor sigxx/sigyy, would be same for sigxx/sigzz on 3d
                        	two carrier (either q,m,n,tau,a or q,n,mu,a):
                        		q1, q2 - charge in e units
                        		m1, m2 - effective mass in terms of electron mass m_e
                        		n1, n2 - carrier density in m^-3
                        		tau1, tau2 - relaxation time in s
                        		a - anisotropy factor sigxx/sigyy, would be same for sigxx/sigzz on 3d
							topological insulator ti (TODO add single or two):
                        		sigxx - surface conductivity in S/m^2
                        		sigxy - diagonal element of the conductivity tensor
                        		a - ratio of longitudinal to transversal conductivity
                        		gSurBulk - ratio of surface to bulk conductivity (fomula acc. to Evidence of Distributed Robust Surface Current Flow in 3D TI)
                        		surface_per - fraction of surface to bulk
                        		type - use single or two carrier for surface TODO
                   			costum:
                        		sig - conductivity tensor in the form of sig_xx, sig_xy, sig_yy, sig_zz as (nx+1)*(ny+1)*(nz+1) arrays
                        	homogeneous:
                        		sxx - conductivity in x-direction
                        		sxy - off-diagonal elements
                        		syy - conductivity in y-direction
                        		szz - conductivity in z-direction
                        	weyl:
                        		cw	- positive fit variable from topological E*B term
                        		A 	- fit term from parabolic contribution
                        		D 	- fit term from antilocation term
                        		s0	- residual conductivity

method addCurrent
	mode    			 -- single point, circle, line
	amp     			 -- amplitude in mA (pos for source, neg for sink)
	params               -- single point (no corners!):
                        		xpos - relative x-position on grid
                        		ypos - relative y-position on grid
                        		zpos - relative z-position on grid (only 3d)
                   			circle (in xy-plane, TODO other planes):
		                        xpos - relative x-position on grid
		                        ypos - relative y-position on grid
		                        zpos - relative z-position on grid (only 3d)
                        		r    - radius in um
		                    line (in xy-plane, TODO other planes):
		                        x1, x2 - relative x-positions on grid
		                        y1, y2 - relative y-positions on grid
		                        z1, z2 - relative z-positions on grid (only 3d)
		                        w - width of the connecting line in um
                        
method addProbes (to measure electric potential difference, optional):
	x1, x2               -- relative x-positions on grid
	y1, y2               -- relative y-positions on grid
	z1, z2               -- relative z-positions on grid (only 3d)
	name                 -- identification
        
method addBatch (to simulate a sweep of different variables, up to two dimensions, optional):
	var                  -- variable to change
	val                  -- range of values to use
        
method addOutput (required):
	mode                 -- contour2d, contour2dV, contour3d (TODO), plot (2d or 3d depending on batches), current measurements
	params               -- contour2d and contour2dV:
                        		n - number of contours
                        		side - which surface to plot (only 3d)
                   			contour3d and contour3dV:
                       			n - number of contours 
                   			plot:
                        		var - what variable to plot
                   			current(TODO):
                        		none
	show                 -- display the result (default is true)
	save                 -- true to save on disk and false to display directly (default is false)
        
method run (solve the diff. eq. with the SOR-algorithm for the settings used):
	w                    -- overrelaxation var to improve convergence of the algorithm (poss. values: 0. to 2., usually best results with 1.2 to 1.4)
	tol                  -- criterium at which convergence is reached
	method               -- either "slow" or "fast", default is fast, slow doesnt require a g++ compiler