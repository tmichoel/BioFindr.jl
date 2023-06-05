```@meta
CurrentModule = Findr
```

# Likelihood ratio tests


Consider correlated genes ``A``, ``B``, and a third variable ``E`` upstream of ``A`` and ``B``, such as a significant eQTL of ``A``. The eQTLs can be obtained either *de novo* using eQTL identification tools such as matrix-eQTL [^Shabalin2012] or kruX [^Qi2014], or from published analyses. Throughout this article, we assume that ``E`` is a significant eQTL of ``A``, whereas extension to other data types is straightforward. We use ``A_i`` and ``B_i`` for the expression levels of gene ``A`` and ``B`` respectively, which are assumed to have gone through the supernormalization, and optionally the genotypes of the best eQTL of ``A`` as ``E_i``, where ``i=1,\dots,n`` across samples. Genotypes are assumed to have a total of ``n_a`` alleles, so ``E_i\in\{0,\dots,n_a\}``. We define the null and alternative hypotheses for a total of six tests, as shown in the table below. 

[^Shabalin2012]: Shabalin AA. [Matrix eQTL: ultra fast eQTL analysis via large matrix operations](https://doi.org/10.1093/bioinformatics/bts163), Bioinformatics, Volume 28, Issue 10, May 2012, Pages 1353–1358, 

[^Qi2014]: Qi J, Foroughi Asl H, Bjorkegren J, Michoel T. [kruX: matrix-based non-parametric eQTL discovery](https://doi.org/10.1186/1471-2105-15-11). BMC Bioinformatics 15, 11 (2014).

!!! note "Six likelihood ratio tests are performed to test the regulation A → B, numbered, named, and defined as shown."
    ![Six likelihood ratio tests](pcbi.1005703.g001.png)

    &copy; 2017 Wang, Michoel. [DOI:10.1371/journal.pcbi.1005703.g001](https://doi.org/10.1371/journal.pcbi.1005703.g001). Reused under the terms of the [Creative Commons Attribution License](http://creativecommons.org/licenses/by/4.0/)

LLRs of every test are calculated separately as follows:

## Correlation test

Define the null hypothesis as ``A`` and ``B`` are independent, and the alternative hypothesis as they are correlated:

```math
{\mathcal H}_{\mathrm{null}}^{\mathrm{(0)}}=A\qquad B,\hspace{4em}{\mathcal H}_{\mathrm{alt}}^{\mathrm{(0)}}=A \to B.
```

The superscript ``(0)`` is the numbering of the test. Both hypotheses are modeled with gene expression levels following bivariate normal distributions, as

```math
\begin{pmatrix}
    A_i\\\\
    B_i
\end{pmatrix} \sim 
N\left(
    \begin{pmatrix}
        0\\\\
        0
    \end{pmatrix},
    \begin{pmatrix}
        \sigma_{A0}^2 & \rho\,\sigma_{A0}\sigma_{B0}\\\\
        \rho\,\sigma_{A0}\sigma_{B0} & \sigma_{B0}^2
    \end{pmatrix}
\right)
```

for ``i=1,\dots,n``. The null hypothesis corresponds to ``\rho=0``.

Maximum likelihood estimators (MLE) for the model parameters ``\rho``, ``\sigma_{A0}``, and ``\sigma_{B0}`` are

```math
\hat{\rho}=\frac{1}{n}\sum_{i=1}^n A_i B_i,\qquad\hat{\sigma}_{A0}=\hat{\sigma}_{B0}=1,
```

and the LLR is simply

```math
\mathrm{LLR}^{\mathrm{(0)}}=-\frac{n}{2}\log(1-\hat{\rho}^2).
```

```@docs
realLLR_col(Y::Matrix{T},Ycol::Vector{T}) where T<:AbstractFloat
```

In the absence of genotype information, we use nonzero correlation between ``A`` and ``B`` as the indicator for ``A\rightarrow B`` regulation, giving the posterior probability

The LLR for the correlation test of a specific gene ``A`` against all other genes ``B`` is implemented as a method of the `realLLR_col` function:

```math
P(A - B)=P(\mathcal{H}_{\mathrm{alt}}^{\mathrm{(0)}} \mid \mathrm{LLR}^{\mathrm{(0)}}).
```


## Primary (linkage) test 

Verify that ``E`` regulates ``A`` from ``{\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}\equiv E\rightarrow A`` and ``{\mathcal H}_{\mathrm{null}}^{\mathrm{(1)}}\equiv E\qquad A``. For ``{\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}``, we model ``E\rightarrow A`` as ``A`` follows a normal distribution whose mean is determined by ``E`` categorically, i.e.

```math
A_i\mid E_i\sim N(\mu_{E_i},\sigma_A^2).
```

From the total likelihood ``p(A\mid E)=\prod_{i=1}^np(A_i\mid E_i)``, we find MLEs for model parameters ``\mu_j,j=0,1,\dots,n_a`` and ``\sigma_A``, as


```math
\hat{\mu}_j=\frac{1}{n_j}\sum_{i=1}^n A_i \delta_{E_i j},\quad \hat{\sigma}_A^2=1-\sum_{j=0}^{n_a}\frac{n_j}{n} \hat{\mu}_j^2,
```

where ``n_j`` is the sample count by genotype category,

```math
n_j \equiv \sum_{i=1}^n \delta_{E_i j}.
```

The Kronecker delta function is defined as ``\delta_{xy}=1`` for ``x=y``, and 0 otherwise. When summing over all genotype values (``j=0,\dots,n_a``), we only pick those that exist (``n_j>0``) throughout this article. Since the null hypothesis is simply that ``A_i`` is sampled from a genotype-independent normal distribution, with MLEs of mean zero and standard deviation one due to the supernormalization (see [General inference algorithm](@ref)), the LLR for test 1 becomes

```math
\mathrm{LLR}^{\mathrm{(1)}}=-\frac{n}{2}\ln\hat{\sigma}_A^2.
```
By favoring a large ``\mathrm{LLR}^{\mathrm{(1)}}``, we select  ``{\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}`` and verify that ``E`` regulates ``A``, with
    
```math 
P(E\rightarrow A)=P({\mathcal H}_{\mathrm{alt}}^{\mathrm{(1)}}\mid\mathrm{LLR}^{\mathrm{(1)}}).
```

The LLR for the primary linkage test can be obtained by selecting the appropriate element from the results of the secondary linkage test.

## Secondary (linkage) test 

The secondary test is identical with the primary test, except it verifies that ``E`` regulates ``B``. Hence repeat the primary test on ``E`` and ``B`` and obtain the MLEs:

```math
\hat{\nu}_j=\frac{1}{n_j}\sum_{i=1}^n B_i \delta_{E_i j},\quad \hat{\sigma}_B^2=1-\sum_{j=0}^{n_a}\frac{n_j}{n}\hat{\nu}_j^2,
```

and the LLR as

```math
\mathrm{LLR}^{\mathrm{(2)}}=-\frac{n}{2}\ln\hat{\sigma}_B^2.
```


``{\mathcal H}_{\mathrm{alt}}^{\mathrm{(2)}}`` is chosen to verify that ``E`` regulates ``B``.


In differential expression analysis, the linkage test is used standalone, and then its LLR for testing a specific group ``E`` against all genes ``B`` is implemented as a method of the `realLLR_col` function:


```@docs
realLLR_col(Y::Matrix{T},E::Vector{S}) where {T<:AbstractFloat, S<:Integer}
```

In causal inference analysis, the LLR for linkage test is computed together with the other tests for efficiency (see below).

## Conditional independence test

Verify that ``E`` and ``B`` are independent when conditioning on ``A``. This can be achieved by comparing 

```math
{\mathcal H}_{\mathrm{alt}}^{\mathrm{(3)}}\equiv B\leftarrow E\rightarrow A\wedge(A\mathrm{\ correlates\ with\ }B)
``` 
against 

```math
{\mathcal H}_{\mathrm{null}}^{\mathrm{(3)}}\equiv E\rightarrow A\rightarrow B.
``` 

LLRs close to zero then prefer ``{\mathcal H}_{\mathrm{null}}^{\mathrm{(3)}}``, and ensure that ``E`` regulates ``B`` only through ``A``:

```math
P(E\perp B\mid A)=P({\mathcal H}_{\mathrm{null}}^{\mathrm{(3)}}\mid\mathrm{LLR}^{\mathrm{(3)}}).
```

For ``{\mathcal H}_{\mathrm{alt}}^{\mathrm{(3)}}``, the bivariate normal distribution dependent on ``E`` can be represented as

```math
\begin{pmatrix}
    A_i\\\\
    B_i
\end{pmatrix}
\mid E_i 
\sim N\left(
    \begin{pmatrix}
        \mu_{E_i}\\\\
        \nu_{E_i}
    \end{pmatrix},
    \begin{pmatrix}
        \sigma_A^2 & \rho\sigma_A\sigma_B \\\\
        \rho\sigma_A\sigma_B & \sigma_B^2
    \end{pmatrix}
    \right).
```
    
For ``{\mathcal H}_{\mathrm{null}}^{\mathrm{(3)}}``, the distributions follow Eq [\[eq-ba-hm-d1\]](#eq-ba-hm-d1){reference-type="ref" reference="eq-ba-hm-d1"}, as well as

```math
B_i\mid A_i\sim N(\rho A_i,\sigma_B^2).
``` 

Substituting parameters ``\mu_j,\nu_j,\sigma_A,\sigma_B,\rho`` of ``{\mathcal H}_{\mathrm{alt}}^{\mathrm{(3)}}`` and ``\mu_j,\rho,\sigma_A,\sigma_B`` of ``{\mathcal H}_{\mathrm{null}}^{\mathrm{(3)}}`` with their MLEs, we obtain the LLR: 

```math
\mathrm{LLR}^{\mathrm{(3)}} =-\frac{n}{2}\ln\left(\hat{\sigma}_A^2\hat{\sigma}_B^2-(\hat{\rho}+\sigma_{AB}-1)^2\right) + \frac{n}{2}\ln\hat{\sigma}_A^2 + \frac{n}{2}\ln(1-\hat{\rho}^2),
``` 

where

```math
\sigma_{AB} \equiv 1-\sum_{j=0}^{n_a}\frac{n_j}{n}\hat{\mu}_j\hat{\nu}_j,
``` 

and ``\hat\rho`` is defined in the [Correlation test](@ref) section.

## Relevance test

Since the indirect regulation ``E\rightarrow B`` tends to be weaker than any of its direct regulation components (``E\rightarrow A`` or ``A\rightarrow B``), we propose to test ``E\rightarrow A\rightarrow B`` with indirect regulation ``E\rightarrow B`` as well as the direct regulation ``A\rightarrow B`` for stronger distinguishing power on weak regulations. We define

```math
{\mathcal H}_{\mathrm{alt}}^{\mathrm{(4)}}\equiv E\rightarrow A\wedge E\rightarrow B\leftarrow A
```

and

```math
{\mathcal H}_{\mathrm{null}}^{\mathrm{(4)}}\equiv E\rightarrow A\qquad B
```
   
This simply verifies that ``B`` is not independent from both ``A`` and ``E`` simultaneously. In the alternative hypothesis, ``B`` is regulated by ``E`` and ``A``, which is modeled as a normal distribution whose mean is additively determined by ``E`` categorically and ``A`` linearly, i.e.

```math
B_i\mid E_i,A_i\sim N(\nu_{E_i}+\rho A_i,\sigma_B^2).
```

We can hence solve its LLR as

```math
\mathrm{LLR}^{\mathrm{(4)}}=-\frac{n}{2}\ln\left(\hat{\sigma}_A^2\hat{\sigma}_B^2-(\hat{\rho}+\sigma_{AB}-1)^2\right)+\frac{n}{2}\ln\hat{\sigma}_A^2,
```

with all MLEs as defined before.

## Pleiotropy test

Based on the positives of the secondary test, we can further distinguish the alternative hypothesis 

```math
{\mathcal H}_{\mathrm{alt}}^{\mathrm{(5)}}\equiv B\leftarrow E\rightarrow A\wedge A\rightarrow B
```

from the null

```math
{\mathcal H}_{\mathrm{null}}^{\mathrm{(5)}}\equiv B\leftarrow E\rightarrow A
```
    
to verify that ``E`` does not regulate ``A`` and ``B`` independently. Its LLR can be solved as

```math
\mathrm{LLR}^{\mathrm{(5)}}=-\frac{n}{2}\ln\left(\hat{\sigma}_A^2\hat{\sigma}_B^2-(\hat{\rho}+\sigma_{AB}-1)^2\right)+\frac{n}{2}\ln\hat{\sigma}_A^2\hat\sigma_B^2.
```

with all MLEs as defined before. Note that this test was called the **controlled** test in the original paper (see figure above), but we now prefer to refer to it as the **pleiotropy** test.

## Implementation

The LLRs for the secondary linkage, conditional independence, relevance, and pleiotrropy tests of a specific gene ``A`` against all other genes ``B`` are implemented as a method of the `realLLR_col` function:


```@docs
realLLR_col(Y::Matrix{T},Ycol::Vector{T},E::Vector{S}) where {T<:AbstractFloat, S<:Integer}
```

The MLEs of the various model parameters are computed in the `llrstats_col` and `groupmeans` functions:

```@docs
llrstats_col
groupmeans
```