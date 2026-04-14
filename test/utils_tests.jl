# Define continuous, integer counts and discrete multiclass dataframes for testing
df_continuous = DataFrame(A = 1.0:5.0, B = 6.0:10.0)
df_counts = DataFrame(A = 1:5, B = 6:10)
df_class = DataFrame(A = ["a", "b", "a", "b", "a"], B = ["a", "b", "c", "a", "b"])

@testset "Scitype coercion tests" begin
    BioFindr.coerce_scitypes!(df_continuous, ScientificTypes.Continuous)
    @test all([all(scitype.(col) .== ScientificTypes.Continuous) for col in eachcol(df_continuous)])
    @test BioFindr.test_scitype(df_continuous, ScientificTypes.Continuous)
    @test !BioFindr.test_scitype(df_continuous, ScientificTypes.Count)
    @test !BioFindr.test_scitype(df_continuous, ScientificTypes.Multiclass)


    BioFindr.coerce_scitypes!(df_counts, ScientificTypes.Count)
    @test all([all(scitype.(col) .== ScientificTypes.Count) for col in eachcol(df_counts)])
    @test BioFindr.test_scitype(df_counts, ScientificTypes.Count)
    @test !BioFindr.test_scitype(df_counts, ScientificTypes.Continuous)
    @test !BioFindr.test_scitype(df_counts, ScientificTypes.Multiclass)

    BioFindr.coerce_scitypes!(df_class, ScientificTypes.Multiclass)
    @test all([all(scitype.(col) .<: ScientificTypes.Finite) for col in eachcol(df_class)])
    @test BioFindr.test_scitype(df_class, ScientificTypes.Multiclass)
    @test !BioFindr.test_scitype(df_class, ScientificTypes.Continuous)
    @test !BioFindr.test_scitype(df_class, ScientificTypes.Count)

    BioFindr.coerce_scitypes!(df_class, ScientificTypes.OrderedFactor)
    @test all([all(scitype.(col) .<: ScientificTypes.Finite) for col in eachcol(df_class)])
    @test BioFindr.test_scitype(df_class, ScientificTypes.OrderedFactor)
    @test !BioFindr.test_scitype(df_class, ScientificTypes.Continuous)
    @test !BioFindr.test_scitype(df_class, ScientificTypes.Count)
    @test !BioFindr.test_scitype(df_class, ScientificTypes.Multiclass)
end

@testset "qvalue tests" begin
    # All equal probabilities → all q-values should be 1 - p
    p = fill(0.8, 10)
    q = BioFindr.qvalue(p)
    @test all(q .≈ 0.2)
    @test all(0 .<= q .<= 1)

    # Monotonicity: q-values must be non-decreasing when probabilities are decreasing
    p2 = Float64[0.9, 0.8, 0.7, 0.6, 0.5]
    q2 = BioFindr.qvalue(p2)
    @test issorted(q2[sortperm(p2, rev=true)])
    @test all(0 .<= q2 .<= 1)

    # Single element
    q3 = BioFindr.qvalue([0.6])
    @test q3 ≈ [0.4]
end

@testset "globalfdr tests" begin
    p = Float64[0.9, 0.8, 0.7, 0.3, 0.1]
    idxs, Q = BioFindr.globalfdr(p, 1.0)
    @test length(idxs) == length(p)
    @test all(0 .<= Q .<= 1)

    idxs_strict, _ = BioFindr.globalfdr(p, 0.2)
    @test length(idxs_strict) <= length(idxs)
end

@testset "globalfdr! DataFrame tests" begin
    dP = DataFrame(Source=["A","A","B"], Target=["B","C","C"], Probability=Float64[0.9,0.5,0.1])
    BioFindr.globalfdr!(dP)
    @test "qvalue" ∈ names(dP)
    @test all(0 .<= dP.qvalue .<= 1)
    @test issorted(dP.qvalue)

    # FDR filtering
    dP2 = DataFrame(Source=["A","A","B"], Target=["B","C","C"], Probability=Float64[0.9,0.5,0.1])
    BioFindr.globalfdr!(dP2; FDR=0.2)
    @test all(dP2.qvalue .<= 0.2)

    # Calling again on already-processed df (with qvalue column) should not add a second qvalue column
    ncols_before = ncol(dP)
    BioFindr.globalfdr!(dP)
    @test ncol(dP) == ncols_before
end

@testset "stackprobs tests" begin
    P = [0.8 0.6; 0.4 0.9; 0.5 0.7]
    colnames = ["A", "B"]
    rownames = ["X", "Y", "Z"]
    dP = BioFindr.stackprobs(P, colnames, rownames; nodiag=false)
    @test dP isa DataFrame
    @test names(dP) == ["Source", "Target", "Probability"]
    @test nrow(dP) == 6   # 3 rows × 2 cols

    # nodiag removes rows where Source == Target
    P2 = [0.9 0.6; 0.4 0.8]
    colnames2 = ["A", "B"]
    rownames2 = ["A", "B"]
    dP2 = BioFindr.stackprobs(P2, colnames2, rownames2; nodiag=true)
    @test nrow(dP2) == 2  # only off-diagonal
    @test all(dP2.Source .!= dP2.Target)
end

@testset "symprobs tests" begin
    P = [1.0 0.8; 0.2 1.0]
    # none: no change
    @test BioFindr.symprobs(P; combination="none") == P
    # prod
    Pprod = BioFindr.symprobs(P; combination="prod")
    @test Pprod ≈ Pprod'
    @test Pprod[1,2] ≈ 0.8 * 0.2
    # mean
    Pmean = BioFindr.symprobs(P; combination="mean")
    @test Pmean ≈ Pmean'
    @test Pmean[1,2] ≈ 0.5
    # anti
    Panti = BioFindr.symprobs(P; combination="anti")
    @test Panti[1,2] + Panti[2,1] ≈ 1.0
    # error on non-square with non-none combination
    @test_throws ErrorException BioFindr.symprobs([0.5 0.6 0.7]; combination="prod")
end

@testset "combineprobs tests" begin
    # P has shape (nX, 4, nA): tests 2-5 for each causal variable
    P = ones(5, 4, 2)
    P[:,1,:] .= 0.9   # test 2
    P[:,2,:] .= 0.8   # test 3
    P[:,3,:] .= 0.7   # test 4
    P[:,4,:] .= 0.6   # test 5

    @test BioFindr.combineprobs(P; combination="none") === P
    PIV = BioFindr.combineprobs(P; combination="IV")
    @test size(PIV) == (5, 2)
    @test all(PIV .≈ 0.9 * 0.6)
    Pmed = BioFindr.combineprobs(P; combination="mediation")
    @test all(Pmed .≈ 0.9 * 0.8)
    Porig = BioFindr.combineprobs(P; combination="orig")
    @test all(Porig .≈ 0.5 * (0.9 * 0.6 + 0.7))
end

@testset "getpairs tests" begin
    dX = DataFrame(A=1.0:3.0, B=4.0:6.0, C=7.0:9.0)
    dG = DataFrame(SNP1=1:3, SNP2=4:6)
    dE = DataFrame(SNP=["SNP1","SNP2"], Gene=["A","C"])
    pairs = BioFindr.getpairs(dX, dG, dE; colG=1, colX=2)
    @test size(pairs) == (2, 2)
    # Check actual index values (pairs sorted by X index: A=1 < C=3)
    @test pairs[1,1] == findfirst(==("SNP1"), names(dG))
    @test pairs[1,2] == findfirst(==("A"), names(dX))

    # Genes not in namesX are silently ignored (returns empty pairs)
    dE_filtered = DataFrame(SNP=["SNP1"], Gene=["MISSING"])
    pairs_empty = BioFindr.getpairs(dX, dG, dE_filtered; colG=1, colX=2)
    @test size(pairs_empty, 1) == 0

    # Error on unmatched SNP name when gene IS found
    dE_bad_snp = DataFrame(SNP=["MISSING_SNP"], Gene=["A"])
    @test_throws ErrorException BioFindr.getpairs(dX, dG, dE_bad_snp; colG=1, colX=2)
end
