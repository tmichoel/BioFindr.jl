Random.seed!(5678)
nA = 3
nB = 100
fB = 0.1
ns = 80
ng = 3
maf = 0.3
bGA = 1.0
bAB = 1.0

G_mat, XA_mat, XB_mat, istarget_mat = BioFindr.generate_test_data(nA, nB, fB, ns, ng, maf, bGA, bAB, true)
# single-source case (nA=1)
G_vec, XA_vec, XB_vec, istarget_vec = BioFindr.generate_test_data(1, nB, fB, ns, ng, maf, bGA, bAB, true)

# Build DataFrames for the DataFrame-based API tests
dXA = DataFrame(XA_mat, :auto)
dXB = DataFrame(XB_mat, :auto)
dG = DataFrame(G_mat, :auto)
BioFindr.coerce_scitypes!(dXA, ScientificTypes.Continuous)
BioFindr.coerce_scitypes!(dXB, ScientificTypes.Continuous)
BioFindr.coerce_scitypes!(dG, ScientificTypes.Multiclass)

@testset "findr_matrix coexpression (1 matrix)" begin
    PP = BioFindr.findr_matrix(XB_mat)
    @test size(PP) == (nB, nB)
    @test all(0 .<= PP .<= 1)
    # diagonal should be 1 (self-probability)
    @test all(PP[i,i] ≈ 1.0 for i in 1:nB)
end

@testset "findr_matrix coexpression subset cols" begin
    cols = [1, 2, 3]
    PP = BioFindr.findr_matrix(XB_mat; cols=cols)
    @test size(PP) == (nB, length(cols))
    @test all(0 .<= PP .<= 1)
end

@testset "findr_matrix coexpression symmetrized" begin
    PP = BioFindr.findr_matrix(XB_mat; combination="prod")
    @test size(PP) == (nB, nB)
    @test PP ≈ PP'
end

@testset "findr_matrix bipartite (2 continuous matrices)" begin
    PP = BioFindr.findr_matrix(XB_mat, XA_mat)
    @test size(PP) == (nB, nA)
    @test all(0 .<= PP .<= 1)
end

@testset "findr_matrix differential expression (integer genotype)" begin
    PP = BioFindr.findr_matrix(XB_mat, G_mat)
    @test size(PP) == (nB, nA)
    @test all(0 .<= PP .<= 1)
end

@testset "findr_matrix causal (with instrument pairs)" begin
    pairGX = hcat(1:nA, 1:nA)  # SNP i instruments gene i
    PP = BioFindr.findr_matrix(XB_mat, G_mat, pairGX; combination="IV")
    @test size(PP) == (nB, nA)
    @test all(0 .<= PP .<= 1)
end

@testset "findr_matrix causal bipartite (X1, X2, G, pairs)" begin
    pairGX = hcat(1:nA, 1:nA)
    PP = BioFindr.findr_matrix(XB_mat, XA_mat, G_mat, pairGX; combination="IV")
    @test size(PP) == (nB, nA)
    @test all(0 .<= PP .<= 1)
end

@testset "findr DataFrame coexpression" begin
    dP = findr(dXB)
    @test dP isa DataFrame
    @test names(dP) == ["Source", "Target", "Probability", "qvalue"]
    @test all(0 .<= dP.Probability .<= 1)
    @test all(0 .<= dP.qvalue .<= 1)
    @test issorted(dP.qvalue)
end

@testset "findr DataFrame coexpression with colnames" begin
    dP = findr(dXB; colnames=names(dXB)[1:2])
    @test dP isa DataFrame
    @test all(dP.Source .∈ Ref(names(dXB)[1:2]))
    @test all(0 .<= dP.Probability .<= 1)
end

@testset "findr DataFrame coexpression FDR filtering" begin
    dP_all  = findr(dXB)
    dP_fdr  = findr(dXB; FDR=0.5)
    @test nrow(dP_fdr) <= nrow(dP_all)
    @test all(dP_fdr.qvalue .<= 0.5)
end

@testset "findr DataFrame bipartite continuous" begin
    dXA_cont = copy(dXA)
    BioFindr.coerce_scitypes!(dXA_cont, ScientificTypes.Continuous)
    dP = findr(dXB, dXA_cont)
    @test dP isa DataFrame
    @test all(0 .<= dP.Probability .<= 1)
end

@testset "findr DataFrame bipartite categorical (Multiclass dG)" begin
    dP = findr(dXB, dG)
    @test dP isa DataFrame
    @test all(0 .<= dP.Probability .<= 1)
end

@testset "findr errors on bad scitype" begin
    dBad = DataFrame(A = ["x","y","x","y"], B = ["a","b","a","b"])
    BioFindr.coerce_scitypes!(dBad, ScientificTypes.Multiclass)
    @test_throws ErrorException findr(dBad)
end
