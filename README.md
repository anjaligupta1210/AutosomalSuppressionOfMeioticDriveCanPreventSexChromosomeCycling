### AutosomalSuppressionOfMeioticDriveCanPreventSexChromosomeCycling

**Data and Code for Autosomal suppression of meiotic drive can prevent sex chromosome cycling**

All Figures are coded in ``Figures.R``


Simulations are in ``TestXXXXX.R``






1. **Simulation filenames based on scenarios presented in paper:**

- **ScenarioA** - ``TestYsupInvasion_InitialPopDoesNotHaveAsup.R`` & ``TestYsupInvasion_InitialPopHasAsup.R``

- **ScenarioB&C** - ``TestAsupInvasion_EqbmPop.R``

- **ScenarioD** - ``TestYsupInvasion_HallCyclingPop.R``

- **EffectOfPopulationSize** - ``TestIfCyclingPeristsUnderFinitePopSize.R``







2. Simulations were run separately for some of the figures -

- **Figure1,2,S1,S2,S3,S4,S5** - ``TestYsupInvasion_InitialPopDoesNotHaveAsup_Fig1-2_S1-5.R`` & ``TestYsupInvasion_InitialPopHasAsup_Fig1-2_S1-5.R``

- **Figure3,S6,S7,S8** - ``TestAsupInvasion_EqbmPop_Fig3_S6-8.R``








3. Mathematica file for the recursions for our model is ``ModelRecursionEquations.nb``




4. Data from simulations is in ``XXXXXPop.csv``

**Dataset filenames based on scenarios presented in paper:**

- **ScenarioA** - ``YsupInvasion_AsupAbsentInInitialPop.csv`` & ``YsupInvasion_AsupPresentInInitialPop.csv`` (Unzip .zip file)

- **ScenarioB** - ``AsupInvasion_EqbmPop.csv`` (Unzip .zip file)

- **ScenarioC** - ``StableCyclingSpaceToTestAsupInvasion.csv``

- **ScenarioD** - ``TestYsupInvasion_HallRegionIV_CyclingPop.csv``

- **EffectofPopulationSize** - ``FinitePopSizeForCyclingPop.csv``



5. Data from simulations were run separately for some of the figures -

- **Figure1,2,S1,S2,S3,S4,S5** - ``YsupInvasion_AsupAbsentInInitialPop_Fig1-2_S1-5.csv.zip`` & ``YsupInvasion_AsupPresentInInitialPop_Fig1-2_S1-5.csv.zip`` (Unzip .zip file)

- **Figure3,S6,S7,S8** - ``AsupInvasion_EqbmPop_Fig3_S6-8.csv``



6. R script to verify that we get the same cycling behavior as Hall 2004 in our model: ``Hall2004CyclingParameterSpace.R``



7. Plots to verify that we get the same cycling behavior as Hall 2004 in our model: ``Hall2004CyclingParameterSpace.pdf``



8. Mathematica file for the verification of reduced version of our model to Wu1983 is ``ModelRecursionEquations_ReductionToWuModel.nb``


9. Mathematica file for the verification of reduced version of our model to Hall2004 is ``ModelRecursionEquations_ReductionToHallModel.nb``
