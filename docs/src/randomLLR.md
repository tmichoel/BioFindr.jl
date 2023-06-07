```@meta
CurrentModule = Findr
```

# Null distributions of the log-likelihood ratios

The null distribution of LLR, ``p(\mathrm{LLR}\mid{\mathcal H}_{\mathrm{null}})``, may be obtained either by simulation or analytically. Simulation, such as random permutations from real data or the generation of random data from statistics of real data, can deal with a much broader range of scenarios in which analytical expressions are unattainable. However, the drawbacks are obvious: simulation can take hundreds of times longer than analytical methods to reach a satisfiable precision. Here we obtained analytical expressions of ``p(\mathrm{LLR}\mid{\mathcal H}_{\mathrm{null}})`` for all the [Likelihood ratio tests](@ref).

## Correlation test

``{\mathcal H}_{\mathrm{null}}^{\mathrm{(0)}}=A\quad B`` indicates no correlation between ``A`` and ``B``. Therefore, we can start from

```math
\tilde{B}_i\sim\mathrm{i.i.d\ }N(0,1).
``` 

In order to simulate the supernormalization step, we normalize ``\tilde{B}_i`` into ``B_i`` with zero mean and unit variance as:

```math
B_i\equiv\frac{\tilde{B}_i-\bar{\tilde{B}}_i}{\sigma_{\tilde{B}}},\hspace{3em}\bar{\tilde{B}}\equiv \frac{1}{n}\sum_{i=1}^n\tilde{B}_i,\hspace{3em}\sigma_{\tilde{B}}^2\equiv\frac{1}{n}\sum_{i=1}^n\left(\tilde{B}_i-\bar{\tilde{B}}\right)^2.
```

Transform the random variables ``\{\tilde{B}_i\}`` by defining  
```math
\begin{aligned}
    X_1 &\equiv \frac{1}{\sqrt{n}}\sum_{i=1}^nA_i\tilde{B}_i,\\
    X_2 &\equiv \frac{1}{\sqrt{n}}\sum_{i=1}^n\tilde{B}_i,\\
    X_3 &\equiv \left(\sum_{i=1}^n\tilde{B}_i^2\right)-X_1^2-X_2^2.
\end{aligned}
``` 

Since ``\tilde{B}_i\sim\mathrm{i.i.d\ }N(0,1)``, we can easily verify that ``X_1,X_2,X_3`` are independent, and  

```math
X_1\sim N(0,1),\hspace{3em}X_2\sim N(0,1),\hspace{3em}X_3\sim\chi^2(n-2).
```

Expressing the LLR for the correlation test (see  [Likelihood ratio tests](@ref)) in terms of ``X_1,X_2,X_3`` gives  

```math
\mathrm{LLR}^{\mathrm{(0)}}=-\frac{n}{2}\ln(1-Y),
``` 

in which

```math
Y\equiv\frac{X_1^2}{X_1^2+X_3}\sim\mathrm{Beta}\Bigl(\frac 12,\frac{n-2}2\Bigr)
```

follows the [Beta distribution](https://en.wikipedia.org/wiki/Beta_distribution).

We define the distribution ``{\mathcal D}(\alpha,\beta)`` as the distribution of a random variable ``Z=-\frac 12\ln(1-Y)`` for ``Y\sim\mathrm{Beta}(\alpha/2,\beta/2)``, i.e.  

```math
Z=-\frac{1}{2}\ln(1-Y)\sim{\mathcal D}(\alpha,\beta).
``` 

The probability density function (PDF) for ``Z\sim{\mathcal D}(\alpha,\beta)`` can be derived as: for ``z>0``, 

```math
p(z\mid \alpha,\beta)=\frac{2}{B(\alpha/2,\beta/2)}\left(1-e^{-2z}\right)^{(\alpha/2-1)}e^{-\beta z},
```

and for ``z\le0``, ``p(z\mid \alpha,\beta)=0``. Here ``B(a,b)`` is the Beta function. Therefore the null distribution for the correlation test is simply 

```math
\mathrm{LLR}^{\mathrm{(0)}}/n\sim{\mathcal D}(1,n-2).
```

## Primary linkage test

``{\mathcal H}_{\mathrm{null}}^{\mathrm{(1)}}=E\hspace{1.6em}A`` indicates no regulation from ``E`` to ``A``. Therefore, similarly with the correlation test, we start from  ``\tilde{A}_i\sim\mathrm{i.i.d\ }N(0,1)`` and normalize them to ``A_i`` with zero mean and unit variance.

The expression of ``\mathrm{LLR}^{\mathrm{(1)}}`` then becomes:

```math
\mathrm{LLR}^{\mathrm{(1)}} = -\frac{n}{2}\ln\left(1-\sum_{j=0}^{n_a}\frac{n_j}{n}\frac{\left(\hat{\tilde{\mu}}_j-\bar{\tilde{A}}\right)^2}{\sigma_{\tilde{A}^2}}\right),
```

where

```math
\hat{\tilde{\mu}}_j\equiv\frac{1}{n_j}\sum_{i=1}^n\tilde{A}_i\delta_{E_ij}.
```

For now, assume all possible genotypes are present, i.e. ``n_j>0`` for ``j=0,\dots,n_a``. Transform ``\{\tilde{A}_i\}`` by defining  

```math
\begin{aligned}
    X_j &\equiv \sqrt{n_j}\,\hat{\tilde{\mu}}_j,\hspace{9em}\mathrm{for\ }j=0,\dots,n_a,\\
    X_{n_a+1} &\equiv \Biggl(\sum_{i=1}^n\tilde{A}_i^2\Biggr)-\Biggl(\sum_{j=0}^{n_a}X_j^2\Biggr).
\end{aligned}
``` 

Then we can similarly verify that ``\{X_i\}`` are pairwise independent, and 

```math 
\begin{aligned}
    X_i &\sim N(0,1),\ \mathrm{for\ }i=0,\dots,n_a,\\
    X_{n_a + 1} &\sim \chi^2(n-n_a-1).
\end{aligned}
```

Again transform ``\{X_i\}`` by defining independent random variables

```math
\begin{aligned}
    Y_1 &\equiv \sum_{j=0}^{n_a}\sqrt{\frac{n_j}{n}}\,X_j\sim N(0,1),\\
    Y_2 &\equiv \left(\sum_{j=0}^{n_a}X_j^2\right)-Y_1^2\sim\chi^2(n_a),\\
    Y_3 &\equiv X_{n_a+1}\sim\chi^2(n-n_a-1).
\end{aligned}
``` 

Some calculation would reveal

```math
\mathrm{LLR}^{\mathrm{(1)}}=-\frac{n}{2}\ln\left(1-\frac{Y_2}{Y_2+Y_3}\right),
```

i.e.

```math
\mathrm{LLR}^{\mathrm{(1)}}/n\sim{\mathcal D}(n_a,n-n_a-1).
```

To account for genotypes that do not show up in the samples, define ``n_v\equiv\sum_{j\in\{j\mid n_j>0\}}1`` as the number of different genotype values across all samples. Then  

```math
\mathrm{LLR}^{\mathrm{(1)}}/n\sim{\mathcal D}(n_v-1,n-n_v).
```

## Secondary linkage test

Since the null hypotheses and LLRs of primary and secondary tests are identical, ``\mathrm{LLR}^{\mathrm{(2)}}`` follows the same null distribution as ``\mathrm{LLR}^{\mathrm{(1)}}``.

## Conditional independence test

The conditional independence test verifies if ``E`` and ``B`` are uncorrelated when conditioning on ``A``, with ``{\mathcal H}_{\mathrm{null}}^{\mathrm{(3)}}=E\rightarrow A\rightarrow B``. For this purpose, we keep ``E`` and ``A`` intact while randomizing ``\tilde{B}_i`` according to ``B``'s correlation with ``A``:

```math
\tilde{B}_i\equiv\hat\rho A_i+\sqrt{1-\hat\rho^2}\,X_i,\hspace{2em}X_i\sim\mathrm{i.i.d\ }N(0,1).
```

Then ``\tilde{B}_i`` is normalized to ``B_i`` as before. The null distribution of ``\mathrm{LLR}^{\mathrm{(3)}}`` can be obtained with similar but more complex computations as before, as

```math
\mathrm{LLR}^{\mathrm{(3)}}/n\sim{\mathcal D}(n_v-1,n-n_v-1).
```

## Relevance test

The null distribution of ``\mathrm{LLR}^{\mathrm{(4)}}`` can be obtained similarly by randomizing ``B_i``, as

```math
\mathrm{LLR}^{\mathrm{(4)}}/n\sim{\mathcal D}(n_v,n-n_v-1).
```

## Pleiotropy test

To compute the null distribution for the controlled test, we start from

```math
\tilde{B}_i=\hat{\nu}_{E_i}+\hat{\sigma}_B X_i,\hspace{2em}X_i\sim N(0,1),
```

and then normalize ``\tilde{B}_i`` into ``B_i``. Some calculation reveals the null distribution as

```math
\mathrm{LLR}^{\mathrm{(5)}}/n\sim{\mathcal D}(1,n-n_v-1).
```

## Implementation

### Null distribution family

The family of distributions ``\mathcal{D}``, termed **LBeta** in the code, is implemented as a continuous, univariate distribution type using the [Distributions](https://juliastats.org/Distributions.jl/stable/) package:

```@docs
LBeta
```

Its PDF is implemented explicitly from the definition above:

```@docs
pdf
logpdf
```

The CDF is implemented by calling the corresponding functions for a [Beta distribution](https://juliastats.org/Distributions.jl/stable/univariate/#Distributions.Beta):

```@docs
cdf
ccdf
logccdf
```

### Test-specific null distributions

To construct a null distribution object for a specific test, use:

```@docs
nulldist
```

Convenience functions are provided that wrap around the [`pdf`](@ref), [`ccdf`](@ref), or [`logccdf`](@ref) functions:

```@docs
nullpdf
nullpval
nulllog10pval
```