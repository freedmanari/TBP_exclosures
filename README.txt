README for "A unifying theoretical framework for tick-borne disease risk to explain conflicting results of exclosure experiments across scales"

Authors: Ari S. Freedman 1*, Simon A. Levin 1, Stephen A. Felt 2, Giulio A. De Leo 3
1 Department of Ecology and Evolutionary Biology, Princeton University, Princeton, NJ 08544
2 Department of Comparative Medicine, Stanford University, Stanford, CA 94305
3 Hopkins Marine Station, Stanford University, Pacific Grove, CA 93950
* Corresponding author; email: arisf@princeton.edu

Abstract:
Disease ecology has focused greatly on determining how changes to biodiversity may drive infectious disease risk for humans. Fencing off experimental areas (exclosures) has been a common experimental approach to assess how removing large-bodied hosts may affect disease risk, especially with tick-borne pathogens (TBPs). However, exclosure experiments have found conflicting results based on the experiment's scale, with smaller exclosures tending to increase tick densities inside the exclosure and larger exclosures tending to decrease tick densities inside. Previously, we have lacked a unifying theoretical framework able to reconcile the results of exclosure experiments across spatial scales. We present a spatially explicit model of TBP risk incorporating tick dispersal by small competent mammal hosts which can enter the exclosure and by large incompetent mammal hosts excluded from the exclosure. Our model reproduces the scale-dependence and spatial patterning observed in past exclosure experiments while elucidating their causal mechanisms. Specifically, the modeled exclosures produce high densities of infected ticks near their boundaries, with the densities decreasing towards the exclosure's center. Empirical results have found lower tick densities at the exclosure's edge than its center, a pattern we demonstrate can also be produced if we additionally allow ticks in their free-living questing stage to disperse.

Files (to be run in this order):
  - lyme_model.R, encodes the spatial partial different equation model for tick-borne pathogen dynamics with host and tick movement, creates Figs. 2, 4, and 5 in the main text and Figs. S1-S3 in the supplement
  - lyme_bin_search.R, to be run after all of the model functions are defined in lyme_model.R, codes in a binary search to find the "inflection exclosure size" for a variety of different parameter sets, makes Fig. 3 in the main text
 
