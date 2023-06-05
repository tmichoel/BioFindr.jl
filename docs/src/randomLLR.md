```@meta
CurrentModule = Findr
```

# Null distributions of the log-likelihood ratios

The null distribution of LLR, ``p(\mathrm{LLR}\mid{\mathcal H}_{\mathrm{null}})``, may be obtained either by simulation or analytically. Simulation, such as random permutations from real data or the generation of random data from statistics of real data, can deal with a much broader range of scenarios in which analytical expressions are unattainable. However, the drawbacks are obvious: simulation can take hundreds of times longer than analytical methods to reach a satisfiable precision. Here we obtained analytical expressions of ``p(\mathrm{LLR}\mid{\mathcal H}_{\mathrm{null}})`` for all the [Likelihood ratio tests](@ref).

## Correlation test

``{\mathcal H}_{\mathrm{null}}^{\mathrm{(0)}}=A\hspace{1.6em}B`` indicates no correlation between ``A`` and ``B``. Therefore, we can start from

```math
\tilde B_i\sim\mathrm{i.i.d\ }N(0,1).
``` 

In order to simulate the supernormalization step, we normalize ``\tilde B_i`` into ``B_i`` with zero mean and unit variance as:

```math
B_i\equiv\frac{\tilde B_i-\bar{\tilde{B}}_i}{\sigma_{\tilde{B}}},\hspace{3em}\bar{\tilde{B}}\equiv \frac{1}{n}\sum_{i=1}^n\tilde{B}_i,\hspace{3em}\sigma_{\tilde{B}}^2\equiv\frac{1}{n}\sum_{i=1}^n\left(\tilde{B}_i-\bar{\tilde{B}}\right)^2.
```