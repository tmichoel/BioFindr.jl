"""
    generate_test_data(nA, nB, fB, ns, ng, maf, bGA, bAB, supernormalize)

Generate test data for Findr with `nA` causal variables, `nB` potential target variables of which a random fraction `fB` are true targets for each causal variable, `ns` samples, `ng` genotype (instrumental variable) groups with minor allele frequence `maf`, and effect sizes `bGA` and `bAB`. Variables are sampled from a linear model with independent Gaussian noise with variance `ϵ` and correlated Gaussian noise with variance `δ` and covariance `δρ`. If `supernormalize` is `true`, the data is supernormalized.
"""
function generate_test_data(nA=2,nB=1000,fB=0.05,ns=100,ng=3,maf=0.3,bGA=1.,bAB=1., ϵ=0.5, δ=0.5, ρ=0.1, supernormalize=true)
    # sample one instrumental variable per causal variable
    G = sum(rand(Bernoulli(maf),ns,nA,ng-1), dims=3)[:,:,1];
    
    # sample causal variables from a linear model defined by its instrument
    XA = bGA .* G + ϵ * randn(ns,nA);
    
    # for each causal variable, randomly assign fB*nB true target variables
    istarget = zeros(Bool,nB,nA);
    for k = axes(istarget,2)
        istarget[shuffle(1:nB)[1:round(Int,fB*nB)],k] .= true;
    end

    # return vectors instead of matrices if nA = 1
    if nA == 1
        G = vec(G);
        XA = vec(XA);
        istarget = vec(istarget);
    end


    # for each target variable, sample from a linear model defined by its causal parents, if any
    XB = bAB .* (XA * istarget') .+ ϵ * randn(ns,nB)    
    
    # add correlated noise to all variables
    μ = zeros(nA+nB);
    Σ = (1-ρ) * I + ρ * ones(nA+nB,nA+nB);

    noise_cor = rand(MvNormal(μ,Σ),ns)';

    XA .+= δ * noise_cor[:,1:nA];
    XB .+= δ * noise_cor[:,nA+1:end];

    # supernormalize data
    if supernormalize
        XA = Findr.supernormalize(XA);
        XB = Findr.supernormalize(XB);
    end
    return G, XA, XB, istarget
end