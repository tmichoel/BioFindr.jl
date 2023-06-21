```@meta
CurrentModule = Findr
```

# Tests to evaluate

Based on the six [Likelihood ratio tests](@ref), we use the following tests and test combinations for the inference of genetic regulations:

## Coexpression analysis

The correlation test is introduced as a benchmark, against which we can compare other methods involving genotype information. Pairwise correlation is a simple measure for the probability of two genes being functionally related either through direct or indirect regulation, or through coregulation by a third factor. Bayesian inference additionally considers different gene roles. Its predicted posterior probability for regulation is ``P_0``.

Correlation analysis can be performed by calling `findr` with one argument, a matrix or dataframe of gene expression values:

```@docs
findr(X::Matrix{T}) where T<:AbstractFloat
findr(dX::T) where T<:AbstractDataFrame
```

## Differential expression analysis

The secondary linkage test is introduced to test association between genetic variants and gene expression levels, and can be used more generally to analyze differential expression of genes across groups defined by any kind of categorical variable. Its predicted posterior probability for differential expression is ``P_2``.

Differential expression analysis can be performed by calling `findr` with two arguments, matrices or dataframes of gene expression and categorical values, respectively:

```@docs
findr(X::Matrix{T}, G::Array{S}) where {T<:AbstractFloat, S<:Integer}
findr(dX::T, dG::T) where T<:AbstractDataFrame
```

## Causal inference

### Mediation

The traditional causal inference test, as explained in [^Chen2007], suggested that the regulatory relation ``E\to A\to B``can be confirmed with the combination of three separate tests: ``E`` regulates ``A``, ``E`` regulates ``B``, and ``E`` only regulates ``B`` through ``A`` (i.e. ``E`` and ``B`` become independent when conditioning on ``A``). They correspond to the primary, secondary, and conditional independence tests respectively. The regulatory relation ``E\to A\to B`` is regarded positive only when all three tests return positive. The three tests filter the initial hypothesis space of all possible relations between ``E``, ``A``, and ``B``, sequentially to ``E\to A`` (primary test), ``E\to A \wedge E\to B`` (secondary test), and ``E\to A\to B \wedge`` (no confounder for ``A`` and ``B``) (conditional independence test). The resulting test is stronger than ``E\to A\to B`` by disallowing confounders for A and B. So its probability can be broken down as

```math
P_{\text{med}} \equiv P_1P_2P_3
```

Findr expects a set of significant eQTLs and their associated genes as input, and therefore ``P_1=1`` is assured and not calculated separately in Findr. Note that ``P_{\text{med}}`` is the estimated local precision, i.e. the probability that tests 2 and 3 are both true. Correspondinly, its local FDR (the probability that one of them is false) is ``1-P_{\text{med}}``.


[^Chen2007]: Chen L, Emmert-Streib F, Storey J. [Harnessing naturally randomized transcription to infer regulatory relationships among genes](https://doi.org/10.1186/gb-2007-8-10-r219). Genome Biol 8, R219 (2007).

### Instrumental variables

The pleiotropy test is introduced to test if an ``E\to B`` association is not independent of the ``E\to A`` association, that is, if an independent pleiotropic effect of ``E`` on both genes can be excluded. If ``E`` regulates ``A`` (is a cis-eQTL for ``A``), and ``E`` regulates ``B``, and  ``B`` and ``A`` are not independent given ``E``, then we can regard ``E`` as a proxy or instrumental variable for ``A`` and infer a regulatory relation ``E\to A\to B`` from the ``E\to B`` association. The three tests verify the hypothesis that ``B \leftarrow E \to A \wedge \lnot(A ⫫ B | E)``, a superset of ``E\to A\to B``.  Its probability can be broken down as

```math
P_{\text{IV}} \equiv P_1P_2P_5
```

As before, ``P_1=1`` is assured and not calculated separately in Findr. ``P_{\text{IV}}`` is again the estimated local precision, i.e. the probability that tests 2 and 5 are both true, and its local FDR (the probability that one of them is false) is ``1-P_{\text{IV}}``.

### Relevance

The relevance test is introduced to address weak interactions that are undetectable by the secondary test from existing data (``P_2`` close to 0). This term still grants higher-than-null significance to weak interactions, and verifies that ``E\to A \wedge (E\to B \vee A - B)``, also a superset of ``E\to A\to B``. Its probability can be broken down as

```math
P_{\text{relev}} \equiv P_1P_4
```

The original Findr paper proposed to combine the instrumental variable and relevance tests in a novel test whose probability can be broken down as

```math
P_{\text{orig}} \equiv \frac{1}{2} P_1 \bigl( P_4 + P_2P_5) = \frac{1}{2}\bigl( P_{\text{relev}} + P_{\text{IV}} \bigr)
```

In the extreme undetectable limit where ``P_2 = 0`` but ``P_4 \neq 0``, the novel test automatically reduces to one half of the relevance test, which assumes equal probability of either direction and assigns half of the relevance test probability to ``A \to B``.

The composite design of the novel test aims not to miss any genuine regulation whilst distinguishing the full spectrum of possible interactions. When the signal level is too weak for tests 2 and 5, we expect ``P_4`` to still provide distinguishing power better than random predictions. When the interaction is strong, ``P_2 P_5`` is then able to pick up true targets regardless of the existence of hidden confounders.

### Implementation

Causal inference can be performed by calling `findr` with three arguments, matrices or dataframes of gene expression and genotype values, and a mapping of matching ``(E,A)`` pairs; the preferred test can be set through the `combination` parameter:

```@docs
findr(X::Matrix{T},G::Matrix{S},pairGX::Matrix{R}; combination="none") where {T<:AbstractFloat, S<:Integer, R<:Integer}
findr(dX::T, dG::T, dE::T; colG=1, colX=2, combination="IV") where T<:AbstractDataFrame
```

## Bipartite causal inference

In the general case, we assume that there is one set of genes, of which the set of ``A``-genes (genes with matching instrument ``E``) is a subset, and that all possible directed regulations are tested. In some situations we are instead searching for a bipartite network from one set of potential causal factors (e.g. micro-RNAs) to another set of potential targets (e.g. protein-coding genes). In this case, causal inference can be performed by calling `findr` with four arguments that include separate matrices or dataframes of expression values for the potential causes and targets:

```@docs
findr(X1::Matrix{T}, X2::Array{T}, G::Array{S}, pairGX::Matrix{R}; combination="none")  where {T<:AbstractFloat, S<:Integer, R<:Integer}
findr(dX1::T, dX2::T, dG::T, dE::T; colG=1, colX=2, combination="IV") where T<:AbstractDataFrame
```