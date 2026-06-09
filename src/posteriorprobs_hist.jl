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

See also [`nulldist`](@ref), [`nullpval`](@ref), [`fit_mixdist_mom`](@ref), [`fit_mixdist_KDE`](@ref).
"""
function _histogram_unequalbins_param_count(n::Integer)
    # Port of tmichoel/findr base/histogram.c: histogram_unequalbins_param_count().
    t = floor(Int, exp(log(float(n)) / 2.5))
    return min(t, 100)
end

function _histogram_equalbins_fromnullpdfs!(binrange::Vector{Float64}, dnull;
                                            nbextend::Int=1000, eps::Float64=1e-5,
                                            maxiter::Int=200)
    # Port of tmichoel/findr base/histogram.c: histogram_equalbins_fromnullpdfs().
    n = length(binrange) - 1
    step = (binrange[end] - binrange[1]) / n
    for i in 2:n
        binrange[i] = binrange[i - 1] + step
    end

    nb = n * nbextend
    vc = zeros(Float64, nb)
    vval = zeros(Float64, nb)
    vstep = ((0:(nbextend - 1)) .+ 0.5) ./ nbextend

    diffmax = 2eps
    niter = 0
    while diffmax > eps && niter < maxiter
        niter += 1

        widths = diff(binrange)
        @inbounds for i in 1:n
            r = ((i - 1) * nbextend + 1):(i * nbextend)
            @views vc[r] .= binrange[i] .+ widths[i] .* vstep
        end

        vval .= pdf.(dnull, vc)
        @inbounds for i in 1:n
            r = ((i - 1) * nbextend + 1):(i * nbextend)
            @views vc[r] .= widths[i] / nbextend
        end
        vval .*= vc
        cumsum!(vval, vval)

        if !(isfinite(vval[end]) && vval[end] > 0)
            break
        end
        vval .*= n / vval[end]

        @inbounds for i in 1:n
            vc[(i - 1) * nbextend + 1] = binrange[i]
        end

        inext = 2
        dnext = 1.0
        dlast = 0.0
        diffmax = 0.0
        i = 1
        while i <= nb && inext <= n
            dnow = vval[i]
            if dnow > dnext
                segstart = ((i - 1) ÷ nbextend) * nbextend + 1
                denom = dnow - dlast
                if denom > 0
                    step = vc[segstart + 1] / denom
                    iaim = min(floor(Int, dnow), n - 1)
                    while (inext - 1) <= iaim && inext <= n
                        t1 = vc[segstart] +
                             vc[segstart + 1] * ((i - 1) % nbextend) +
                             ((inext - 1) - dlast) * step
                        denomw = binrange[inext + 1] - binrange[inext]
                        if denomw > 0
                            diffmax = max(diffmax, abs(t1 - binrange[inext]) / denomw)
                        end
                        binrange[inext] = t1
                        inext += 1
                    end
                    dnext = inext - 1
                end
            end
            dlast = dnow
            i += 1
        end
    end

    return binrange
end

function _histogram_unequalbins_fromequalbins!(binrange::Vector{Float64})
    # Port of tmichoel/findr base/histogram.c: histogram_unequalbins_fromequalbins().
    n = length(binrange) - 1
    twidth = binrange[end] - binrange[1]
    width = log(twidth / n)

    for i in (n + 1):-1:2
        binrange[i] -= binrange[i - 1]
    end

    t2 = n - 1
    for i in 3:n
        w = log(binrange[i]) * (t2 - (i - 2)) + width * (i - 2)
        binrange[i] = exp(w / t2)
    end
    binrange[end] = twidth / n

    s = sum(abs, @view binrange[2:end])
    binrange[2:end] .*= twidth / s
    for i in 2:n
        binrange[i] += binrange[i - 1]
    end
    binrange[end] = binrange[1] + twidth

    return binrange
end

function _fit_mixdist_hist_edges(dnull, dmax::Float64, nd::Int)
    # Port of tmichoel/findr pij/nullhist.c: pij_nullhist_single() + base/histogram.c helpers.
    nbin = _histogram_unequalbins_param_count(nd)
    nbin < 5 && return nothing

    edges = collect(range(0.0, dmax * (1 + 1e-6), length=nbin + 1))
    _histogram_equalbins_fromnullpdfs!(edges, dnull)
    _histogram_unequalbins_fromequalbins!(edges)
    return edges
end

function fit_mixdist_hist(llr,ns,ng=1,test=:corr)
    # Port of tmichoel/findr pij/llrtopij.c + pij/nullhist.c: setup and null histogram construction.
    dnull = nulldist(ns, ng, test)

    dmax = maximum(llr)
    if !(isfinite(dmax) && dmax > 0)
        @debug "Encountered non-positive or non-finite dmax in histogram fit, returning null distribution" dmax
        return zeros(length(llr))
    end

    nd = length(llr)
    edges = _fit_mixdist_hist_edges(dnull, dmax, nd)
    if isnothing(edges)
        @debug "Histogram bin count is too small for conversion, returning null distribution"
        return zeros(length(llr))
    end

    nbin = length(edges) - 1
    widths = diff(edges)

    # Port of tmichoel/findr pij/llrtopij.c: construct real histogram over null-derived unequal bins.
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
    if total == 0
        @debug "Histogram binning produced zero total counts, returning null distribution"
        return zeros(length(llr))
    end
    realdens = counts ./ (total .* widths)
    nb = findlast(>(0), counts)
    if isnothing(nb)
        @debug "No non-zero histogram bins found, returning null distribution"
        return zeros(length(llr))
    end
    if nb == 1
        @debug "Only one non-zero histogram bin found, returning null distribution"
        return zeros(length(llr))
    end

    # Port of tmichoel/findr pij/llrtopij.c: fill zero-density gaps in real histogram.
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

    # Port of tmichoel/findr pij/llrtopij.c: estimate posterior scale from null/real density ratio.
    null_over_real = nulldens ./ realwork
    minpos = argmax(null_over_real)
    ndens = inv(null_over_real[minpos])
    if ndens == 0
        @debug "Histogram scaling produced ndens=0, returning ones posterior"
        return ones(length(llr))
    end

    # Port of tmichoel/findr pij/llrtopij.c: compute raw true-ratio histogram.
    ratio = 1 .- ndens .* null_over_real
    ratio = clamp.(ratio, 0.0, 1.0)
    ratio[1:minpos] .= 0.0

    # Port of tmichoel/findr pij/llrtopij.c: extend ratio to trailing bins.
    pbin = zeros(nbin)
    pbin[1:nb] .= ratio
    if nb < nbin
        pbin[(nb+1):end] .= pbin[nb]
    end

    # Port of tmichoel/findr pij/llrtopij.c: force monotonicity via upper/lower envelopes.
    upper = copy(pbin)
    for i in 2:nbin
        upper[i] = max(upper[i], upper[i-1])
    end
    lower = copy(pbin)
    for i in (nbin-1):-1:1
        lower[i] = min(lower[i], lower[i+1])
    end
    pbin .= 0.5 .* (upper .+ lower)

    # Port of tmichoel/findr pij/llrtopij.c: smooth monotone ratio histogram.
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

    # Port of tmichoel/findr pij/llrtopij.c: final clamp and interpolation to LLR inputs.
    pbin = accumulate(max, clamp.(pbin, 0.0, 1.0))

    mids = 0.5 .* (edges[1:(end-1)] .+ edges[2:end])
    knots = vcat(edges[1], mids, edges[end])
    vals = vcat(pbin[1], pbin, pbin[end])
    itp = LinearInterpolation(knots, vals, extrapolation_bc=Flat())
    pp = itp.(llr)
    pp[llr .<= 0] .= 0
    return clamp.(pp, 0.0, 1.0)
end
