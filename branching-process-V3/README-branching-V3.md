# Branching process V3
## Started ~13 June

Changes made in V3 with respect to V1 and V2 are suggested in Issue #12. V3 is essentially a merger of V1 and V2, with shared effective carrying capacity calculated by subtracting the weighted average of all other population sizes from the focal carrying capacity, with the interaction coefficients being the weights of the corresponding population size.

The carrying capacity is back in the `cell_pop` array, which is now the same as in V1. `alpha_matrix` and the code of setting pairwise `alpha` values is taken from V2.
New function: `get_effective_k()`-takes current populations and pairwise `alpha` values, and returns new carrying capacities for all populations.
