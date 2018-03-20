To use this module either start by modifying one of the example tests or use the __init__.py file
or by importing the grid module and following those directions:
	initialize a grid with desired dimensions via    	grid = new Grid([lx, ly], ny=20)
	set the conductivity via 						 	grid.setConductivity(sigtype=?, .....)
	add currents and voltage probes  					grid.addCurrent() or grid.addProbe()
	if necessary add sweeping variables 				grid.addBatch()
	add desired output/s 								grid.addOutput
	start the simulation with run 						grid.run()

Usual methods of execution:
	if in freya folder write in command line: python -m core.__init__
	or in core folder: python -m __init__


To make the module more accessable you can add path to PYTHONPATH env variable