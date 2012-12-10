QDC (quick direct-method controlled) is an optimized exact implementation of the Gillespie's direct-method.
It is designed for biochemical simulations when there is the need of dynamic parameters whose values can change during the simulation.

CHANGELOG:

2.0.0
-debug mode
-a log file (<model_name>_log.txt) is generated at each run
-stop conditions block: the user can set up conditions to stop the simulation at a specific time point or when the number of molecules of a species reaches a threshold.

1.3.6
-fixed bug in .csv results file generation

1.3.5
-fixed bug libsbml_wrap.cpp (installer)
-the installer works on Ubuntu 11.10
-fixed bug path libsbml in validateSBML.py 

1.3.4
-addition of a guide file pdf for the UI
-introduced several new controls in the parsing stage
-fixed type size check bug
-fixed instant reaction chained execution bug
-modified the selection algorithm for the instant reaction to execute if multiple are selectable: from the selection of the first parsed to the selection of a random one
-fixed a bug for an incorrect data.reactioncounts header line for some instant reaction
