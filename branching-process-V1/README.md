# Branching evolution model of cancer incidence

This code simulates a branching process of somatic evolution, with stochastic mutation accumulation and cell
populations undergo logistic growth. Competition between populations is captured through a shared carrying capacity.

## V1 [started on 5 March]
- Intrinsic growth rate given by random samples from a normally-distributed 'g'. Carrying capacity calculated as sum(all other existing populations)-carrying capacity of the focal population.
- Since growth/transition rates are samples from the same distribution, this is the context-independent version of the model.
- Working draft finalised on 17 April.

## np_testing [re-started on 1 May]
- Testing of the whole simulation for 10 values each of p and n.
- Produces absolute counts for each combination of n and p value, and age-adjusted incidence rates standardised to the US 2000 standard population.
