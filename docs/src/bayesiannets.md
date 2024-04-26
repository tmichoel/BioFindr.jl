```@meta
CurrentModule = BioFindr
```

# Bayesian network learning

[BioFindr](https://github.com/tmichoel/BioFindr.jl) implements the method described in the paper [High-dimensional bayesian network inference from systems genetics data using genetic node ordering](https://doi.org/10.3389/fgene.2019.01196)[^Wang2017] to learn a Bayesian network using a dataframe of posterior probabilities as prior edge weights.

## DAG reconstruction

The first step in the algorithm is to convert a dataframe of posterior probabilities to a directed acyclic graph (DAG). [BioFindr][1] implements the original greedy algorithm where edges are added one by one in descending order of posterior probability and edges that would introduce a cycle in the [`dagfindr_greedy_edges!`](@ref) function. Two additional methods [`dagfindr_greedy_insertion!`](@ref) and [`dagfindr_heuristic_sort!`](@ref) developed by [Kenneth Stoop](https://research.ugent.be/web/person/kenneth-stoop-0/en) and [Pieter Audenaert](https://research.ugent.be/web/person/pieter-audenaert-0/en) in [this paper](https://biblio.ugent.be/publication/8772612) are also implemented. The [`dagfindr!`](@ref) is the main user interface function.

```@docs
dagfindr!
dagfindr_greedy_edges!
dagfindr_heuristic_sort!
dagfindr_greedy_insertion!
greedy_insertions!
edge_weights
names_to_index!
```

[^Wang2019]: Wang L, Audenaert P, Michoel T (2019) [High-dimensional bayesian network inference from systems genetics data using genetic node ordering](https://doi.org/10.3389/fgene.2019.01196). Frontiers in Genetics, Special Topic Machine Learning and Network-Driven Integrative Genomics, 10, 1196.
