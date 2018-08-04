class SolverGenerator():
	def generate_plotter(self, solver_num, solver):
		plotters = solver["plotters"]
		for j, plotter in enumerate(plotters):
			print("Generating plotter[%d] = %s for solver" % (j, plotter["name"]))
		
			# stub, do something
	
	def __init__(self,spec,verbose):
		solvers = spec["solvers"]
		
		for i, solver in enumerate(solvers):
			print("Generating solver[%d] = %s..." % (i, solver["name"]))
			
			# stub, do something
			
			self.generate_plotter(i, solver)
