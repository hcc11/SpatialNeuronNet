This folder contains Matlab (R2015a) and C codes for the manuscript: 
Huang C, Ruff DA, Pyle R, Rosenbaum R, Cohen MR and Doiron B (2019) “Circuit models of low dimensional shared variability in cortical networks”, Neuron 101, 337-348, doi: https://doi.org/10.1016/j.neuron.2018.11.034. 

++++++++++++++++++++++++++++++++

To use, first compile all the C codes with the mex compiler in Matlab. Specifically, in the Matlab command line, run the following commands:
mex EIF1DRFfastslowSyn.c
mex spktime2count.c

Scripts of different network simulations are named Simulation_FigX.m, where X is the Figure number. 

Scripts for generating corresponding figures are MakeFigureX.m.  

RF2D3layer.m  is the main simulation function. It contains default parameter values and uses the mex file EIF1DRFfastslowSyn.c for integration. 

demo.m is a demonstration code for two-layer network simulations.   

FA.m is a script for factor analysis, which uses the codes from /fa_Yu/, courtesy of Byron Yu. 

figtools/ folder contains utility functions for figure formatting, courtesy of Bernhard Englitz. 

data/ folder contains part of the data used for plotting figures.  

corr_d.m computes correlation as a function of distance. 

spkcounts.m converts spike time data to spike counts. 

raster2D_ani.m generates movie of spike rasters from 2D spatial networks. 

+++++++++++++++++++++++++++++++++++

One simulation of a two-layer network for 20 sec takes about 1.5 hours CPU time and under 3 gb memory. 
(Simulation_Fig2.m)

One simulation of a three-layer network for 20 sec takes about 3.5 hours CPU time and under 3 gb memory. 
(Simulation_Fig3.m, Simulation_Fig4.m, Simulation_Fig6.m)

