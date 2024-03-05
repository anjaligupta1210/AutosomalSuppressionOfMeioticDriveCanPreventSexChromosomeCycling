# AutosomalSuppressionOfMeioticDriveCanPreventSexChromosomeCycling
Data and Code for Autosomal suppression of meiotic drive can prevent sex chromosome cycling

All Figures are coded in Figures.R


Simulations are in TestXXXXX.R


Simulation filenames based on scenarios presented in paper: 

ScenarioA - TestYsupInvasion_InitialPopDoesNotHaveAsup.R & TestYsupInvasion_InitialPopHasAsup.R
ScenarioB&C - TestAsupInvasion_EqbmPop.R
ScenarioD - TestYsupInvasion_HallCyclingPop.R
EffectOfPopulationSize - TestIfCyclingPeristsUnderFinitePopSize.R


Mathematica file for the recursions for our model is ModelRecursionEquations.nb


Data from simulations is in XXXXXPop.csv

Dataset filenames based on scenarios presented in paper:

ScenarioA - YsupInvasion_AsupAbsentInInitialPop.csv & YsupInvasion_AsupPresentInInitialPop.csv
ScenarioB - AsupInvasion_EqbmPop.csv
ScenarioC - StableCyclingSpaceToTestAsupInvasion.csv
ScenarioD - TestYsupInvasion_HallRegionIV_CyclingPop.csv
EffectofPopulationSize - FinitePopSizeForCyclingPop.csv



R script to verify that we get the same cycling behavior as Hall 2004 in our model: Hall2004CyclingParameterSpace.R
Plots to verify that we get the same cycling behavior as Hall 2004 in our model: Hall2004CyclingParameterSpace.pdf

Mathematica file for the verification of reduced version of our model to Wu1983 is ModelRecursionEquations_ReductionToWuModel.nb
Mathematica file for the verification of reduced version of our model to Hall2004 is ModelRecursionEquations_ReductionToHallModel.nb
