# Branching model of cancer incidence

This code simulates a branching process of somatic evolution, with stochastic mutation accumulation and cell
populations undergo logistic growth. Competition between populations is captured through a shared carrying capacity.

V1 [started on 5 March]:
- Intrinsic growth rate given by random samples from a normally-distributed 'g'. Carrying capacity calculated as sum(all other existing populations)-carrying capacity of the focal population.
- Since growth/transition rates are samples from the same distribution, this is the context-independent version of the model.
- Working draft finalised on 17 April.

V2 [started on 15 May]:
- Different form of branching model, with a common carrying capacity that is constant across cell populations, along with density-dependent interaction coefficients, `alpha`, that carry the effect of inter-clone competition.
- Individual values of carrying capacity do not change between populations.
- `alpha` values are assigned pairwise, to each combination of cell populations, through the function, `set_alpha`. Each population gets an array of `alpha` values corresponding to each pairwise interaction.
    - New values are generated with every new population, both for the new population, and for the new interaction that is added to the old populations. Therefore, if there are `n` populations including the new mutant, all `alpha` arrays would be `n` elements long, with each array updated with the corresponding number of values.
    - **Assumption**: `alpha` values of reciprocal interactions are independent, that is, `alpha(i,j)` and `alpha(j,i)` are assigned independently of each other.
- `k_calculate()` from V1 is dropped in V2, and `generate_mutant_pop()` must include updating the `alpha` matrix with corresponding interaction coefficients.
- `grow_logistically()` calculates size as the sum of intrinsic logistic growth of current population and effect of competitive or facilitative interactions of other populations.
- All other functions are retained without change from V1.

## Important considerations

- The assumption stated above needs to be examined carefully.
- Are there any *a priori* reasons for the shape of the `alpha` distribution to affect cancer incidence?
- `alpha_matrix` will be an interesting parameter to study; at the end of the simulation, it provides a snapshot of the ecological interaction landscape of tumours in the entire population. Technically, with the shared carrying capacity model, we have shown that cancer progression in the model does not *require* mutualistic interactions between clones. There might however be questions worth pursuing about the interactions of a wider range of ecological behaviours and other basic parameters like the mutation rate.
- Ultimately, this might just end up being another, equivalent way of presenting a more selection-centric rather than a mutation-centric picture.
