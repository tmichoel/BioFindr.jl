"""
    fit_mixdist_hist(llr,ns,ng=1,test=:corr)

Return posterior probabilities for a vector of log-likelihood ratio values `llr` for a given BioFindr `test` with sample size `ns` and number of genotype groups `ng` using a histogram-based conversion translated from the original Findr implementation.

The input variable `test` can take the values:

- `:corr` - **correlation test** (test 0)
- `:link` - **linkage test** (test 1/2)
- `:med` - **mediation test** (test 3)
- `:relev` - **relevance test** (test 4)
- `:pleio` - **pleiotropy test** (test 5)

With two input arguments, the correlation test with `ns` samples is used. With three input arguments, or with four arguments and `test == :corr`, the correlation test with `ns` samples is used and the third argument is ignored.

See also [`pi0est`](@ref), [`nulldist`](@ref), [`nullpval`](@ref), [`fit_mixdist_mom`](@ref), [`fit_mixdist_KDE`](@ref).
"""
function fit_mixdist_hist(llr,ns,ng=1,test=:corr)
    dnull = nulldist(ns, ng, test)
    π0 = pi0est(nullpval(llr, ns, ng, test))
    if π0 == 1
        @debug "Estimated prior probability π0=1, returning null distribution"
        return zeros(length(llr))
    end

    dmax = maximum(llr)
    if !(isfinite(dmax) && dmax > 0)
        return zeros(length(llr))
    end

    nd = length(llr)
    nbin = clamp(round(Int, sqrt(nd)), 10, 200)
    edges = collect(range(0.0, dmax * (1 + 1e-6), length=nbin+1))
    widths = diff(edges)

    counts = zeros(Float64, nbin)
    @inbounds for x in llr
        if x < edges[1]
            continue
        elseif x >= edges[end]
            counts[end] += 1
        else
            idx = searchsortedlast(edges, x)
            counts[clamp(idx, 1, nbin)] += 1
        end
    end
    total = sum(counts)
    total == 0 && return zeros(length(llr))
    realdens = counts ./ (total .* widths)
    nb = findlast(>(0), counts)
    isnothing(nb) && return zeros(length(llr))
    nb == 1 && return zeros(length(llr))

    realwork = copy(realdens[1:nb])
    for i in 2:nb
        if realwork[i] == 0
            realwork[i] = realwork[i-1]
        end
    end

    nulldens = similar(realwork)
    @inbounds for i in eachindex(realwork)
        nulldens[i] = (cdf(dnull, edges[i+1]) - cdf(dnull, edges[i])) / widths[i]
    end

    null_over_real = nulldens ./ realwork
    minpos = argmax(null_over_real)
    ndens = inv(null_over_real[minpos])
    if ndens == 0
        return ones(length(llr))
    end

    ratio = 1 .- ndens .* null_over_real
    ratio = clamp.(ratio, 0.0, 1.0)
    ratio[1:minpos] .= 0.0

    pbin = zeros(nbin)
    pbin[1:nb] .= ratio
    if nb < nbin
        pbin[(nb+1):end] .= pbin[nb]
    end

    upper = copy(pbin)
    for i in 2:nbin
        upper[i] = max(upper[i], upper[i-1])
    end
    lower = copy(pbin)
    for i in (nbin-1):-1:1
        lower[i] = min(lower[i], lower[i+1])
    end
    pbin .= 0.5 .* (upper .+ lower)

    ncut = min(nbin ÷ 2, 50)
    if ncut > 0
        sigma = min(ncut / 3, 30.0)
        if sigma > 0
            k = collect(-ncut:ncut)
            kernel = exp.(-(k ./ sigma) .^ 2)
            kernel ./= sum(kernel)
            padded = vcat(fill(pbin[1], ncut), pbin, fill(pbin[end], ncut))
            smooth = similar(pbin)
            for i in eachindex(smooth)
                smooth[i] = dot(@view(padded[i:(i+2ncut)]), kernel)
            end
            pbin .= smooth
        end
    end

    pbin = accumulate(max, clamp.(pbin, 0.0, 1.0))

    mids = 0.5 .* (edges[1:(end-1)] .+ edges[2:end])
    knots = vcat(edges[1], mids, edges[end])
    vals = vcat(pbin[1], pbin, pbin[end])
    itp = LinearInterpolation(knots, vals, extrapolation_bc=Flat())
    pp = itp.(llr)
    pp[llr .<= 0] .= 0
    return clamp.(pp, 0.0, 1.0)
end
