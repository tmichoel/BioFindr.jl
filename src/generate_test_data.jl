"""
    generate_test_data(nA, nB, ns, ng)

Generate test data for Findr with `nA` causal variables, `nB` potential target variables of which a random fraction `fB` are true targets for each causal variable, `ns` samples, `ng` genotype (instrumental variable) groups with minor allele frequence `maf`, and effect sizes `bGA` and `bAB`.
"""
function generate_test_data(nA=2,nB=1000,fB=0.05,ns=100,ng=3,maf=0.3,bGA=1.,bAB=1.)
    # sample one instrumental variable per causal variable
    G = sum(rand(Bernoulli(maf),ns,nA,ng-1), dims=3)[:,:,1];
    
    # sample causal variables from a linear model defined by its instrument
    XA = zeros(ns,nA);
    for k = axes(XA,2)
        XA[:,k] .= bGA*G[:,k] .+ randn(ns);
    end

    # for each causal variable, randomly assign fB*nB true target variables
    istarget = BitArray(undef,nA,nB);
    for k = axes(istarget,1)
        istarget[k,shuffle(1:nB)[1:round(Int,fB*nB)]] .= true;
    end

    # for each target variable, sample from a linear model defined by its causal parents
    XB = zeros(ns,nB);
    for k = axes(XB,2)
        if sum(istarget[:,k]) == 0
            XB[:,k] .= randn(ns);
        else
            XB[:,k] .= bAB*sum(XA[:,istarget[:,k]],dims=2) + randn(ns);
        end
    end

    return G, XA, XB, istarget
end