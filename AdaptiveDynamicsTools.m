(* ::Package:: *)

PairwiseInvasibilityPlot::usage = "Shows a pairwise invasibility plot (wraps around the PairwiseInvasibilityFunction).";
PairwiseInvasibility::usage = "Function to perform a numerical invasion analysis across resident and mutant trait values. Returns a table with invasion fitness for each pair of mutant-resident trait value.";
PlotSelectionGradient::usage = "Plot the selection gradient along trait values (wraps around the ComputeSelectionGradients function)";
ComputeSelectionGradients::usage = "Function to calculate the selection gradient for a range of resident trait values. Returns a table with the selection gradient for each resident trait value.";
ZeroCrossings::usage = "Returns a list of positions between which the list crosses zero.";
EvaluateSingularStrategies::usage = " to numerically find and evaluate singular strategies given some parameters. Returns a list of singular strategies with their convergence and evolutionary stability and equilibrium population densities.";


(* Show the pairwise invasibility plot (wraps around the PairwiseInvasibilityFunction) *)
PairwiseInvasibilityPlot[fitness_, demography_, ninit_, parameters_, xrange_, opts_:OptionsPattern[]] := 

	(* Shows a pairwise invasibility plot *)

	(* fitness = an expression for the invasion fitness function *)
	(* demography = the equilibrium condition(s) for the ecological dynamics *)
	(* ninit = initial condition(s) for the ecological equilibrium numerical search *)
	(* xrange = range of trait values explored and step size *)
	(* opts = graphical options to be passed to ListContourPlot e.g. the fitness threshold defining the zero-relative fitness isocline *)

	Block[{data},
	
		data = PairwiseInvasibility[fitness, demography, ninit, parameters, xrange];
		ListContourPlot[data, FrameLabel->{"Resident trait value", "Mutant trait value"}, PlotLabel->"Invasion fitness", PlotLegends->Automatic, opts]
		
	]


(* Function to perform a numerical invasion analysis across resident and mutant trait values *)
PairwiseInvasibility[fitness_, demography_, ninit_, parameters_, xrange_] := 

	(* Returns a table with invasion fitness for each pair of mutant-resident trait value *)

	(* fitness = an expression for the invasion fitness function *)
	(* demography = the equilibrium condition(s) for the ecological dynamics *)
	(* ninit = initial condition(s) for the ecological equilibrium numerical search *)
	(* xrange = range of trait values explored and step size *)

	Block[{xvalues, combinations, equilibria, numDemography, numFitness, fitnessValues},
	
		(* Trait values to test *)
		xvalues = Range[xrange[[1]], xrange[[2]], xrange[[3]]];
		 
		(* Mutant-resident trait combinations *)
		combinations = Tuples[xvalues, 2];
		 
		(* Parametrization *)
		numDemography = demography /. parameters;
		numFitness = fitness /. parameters;
		 
		(* Solve the ecological equilibrium for each resident trait value *)
		equilibria = Reap[Map[Sow[FindRoot[numDemography /. x->#, ninit]]&, xvalues]][[2]][[1]];
		
		(* Repeat each equilibrium to match the number of combinations *)
		equilibria = Map[ConstantArray[#, Length[xvalues]]&, equilibria] // Catenate;
		
		(* Compute invasion fitness across trait space *)
		fitnessValues = Reap[MapThread[Sow[numFitness /. #2 /. y->#1[[2]] /. x->#1[[1]]]&, {combinations, equilibria}]][[2]][[1]];

		(* Package trait and fitness values into a table *)
		Join[combinations, List /@ fitnessValues, 2]
		
	]


(* Plot the selection gradient along trait values (wraps around the ComputeSelectionGradients function) *)
PlotSelectionGradient[gradient_, demography_, ninit_, parameters_, xrange_, opts_:OptionsPattern[]] := 

	(* Shows a plot of the selection gradient *)

	(* gradient = an expression for the selection gradient *)
	(* demography = the equilibrium condition(s) for the ecological dynamics *)
	(* ninit = initial condition(s) for the ecological equilibrium numerical search *)
	(* xrange = range of trait values explored and step size *)
	(* opts = graphical options to be passed to ListContourPlot e.g. color *)

	Block[{xvalues, gradientValues, data},
	
		xvalues = Range[xrange[[1]], xrange[[2]], xrange[[3]]];
		gradientValues = ComputeSelectionGradients[gradient, demography, ninit, parameters, xrange];
		data = {xvalues, gradientValues} // Transpose;
		ListPlot[data, AxesLabel->{"Resident trait value", "Selection gradient"}, opts]
		
	]


(* Function to calculate the selection gradient for a range of resident trait values *)
ComputeSelectionGradients[gradient_, demography_, ninit_, parameters_, xrange_] := 

	(* Returns a table with the selection gradient for each resident trait value *)

	(* gradient = an expression for the selection gradient *)
	(* demography = the equilibrium condition(s) for the ecological dynamics *)
	(* ninit = initial condition(s) for the ecological equilibrium numerical search *)
	(* xrange = range of trait values explored and step size *)

	Block[{xvalues, combinations, equilibria, numDemography, numGradient, gradientValues},
	
		(* Trait values to test *)
		xvalues = Range[xrange[[1]], xrange[[2]], xrange[[3]]];
		 
		(* Parametrization *)
		numDemography = demography /. parameters;
		numGradient = gradient /. parameters;
		 
		(* Solve the ecological equilibrium for each resident trait value *)
		equilibria = Reap[Map[Sow[FindRoot[numDemography /. x->#, ninit]]&, xvalues]][[2]][[1]];
		
		(* Compute selection gradient across trait space *)
		gradientValues = Reap[MapThread[Sow[numGradient /. #2 /. x->#1]&, {xvalues, equilibria}]][[2]][[1]]
		
	]


(* Find where a list crosses zero *)
ZeroCrossings[l_List] := 

	(* Returns a list of positions between which the list crosses zero *)
	
	(* Taken from David Zhang *)
	(* https://mathematica.stackexchange.com/questions/10640/find-zero-crossing-in-a-list *)

	Block[{t, u, v},
		t = {Sign[l], Range[Length[l]]} // Transpose; (* List of -1, 0, 1 only *)
		u = Select[t, First[#] != 0 &];               (* Ignore zeros *)
		v = SplitBy[u, First];                        (* Group into runs of + and - values *)
		{Most[Max[#[[All, 2]]] & /@ v], Rest[Min[#[[All, 2]]] & /@ v]} // Transpose
	]


(* Function to numerically find and evaluate singular strategies given some parameters *)
EvaluateSingularStrategies[gradient_, curvature_, demography_, ninit_, parameters_, xrange_] :=
	
	(* Returns a list of singular strategies with their convergence and evolutionary stability and equilibrium population densities *) 

	(* gradient = an expression for the selection gradient *)
	(* curvature = an expression for the curvature of the fitness function at the singular point *)
	(* demography = the equilibrium condition(s) for the ecological dynamics *)
	(* ninit = initial condition(s) for the ecological equilibrium numerical search *)
	(* xrange = range of trait values explored and step size *)
	
	Block[{xvalues, gradientValues, crossingPoints, singularStrategies, convergences, densities, curvatures},

		(* Trait values to test *)
		xvalues = Range[xrange[[1]], xrange[[2]], xrange[[3]]];

		(* Compute selection gradients along trait axis *)
		gradientValues = ComputeSelectionGradients[gradient, demography, ninit, parameters, xrange];
		
		(* Find where the selection gradient is zero *)
		crossingPoints = ZeroCrossings[gradientValues];
		
		(* What are those singular strategies? *)
		singularStrategies = Map[Mean[xvalues[[#]]]&, crossingPoints];
		
		(* Equilibrium population densities at the singular strategies *)
		densities = Map[FindRoot[demography /. parameters /. x->#, ninit]&, singularStrategies];
		
		(* Are those singular strategies convergence-stable? *)
		convergences = Map[Differences[gradientValues[[#]]]&, crossingPoints];
		
		(* Are they evolutionarily stable? *)
		curvatures = MapThread[curvature /. parameters /. x->#1 /. #2&, {singularStrategies, densities}];
		
		(* Extract numbers from densities *)
		densities = Map[#[[;;, 2]]&, densities];
		
		Map[Flatten[#]&, {singularStrategies, convergences, curvatures, densities} // Transpose]
		
	]
