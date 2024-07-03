# Define continuous, integer counts and discrete multiclass dataframes for testing
df_continuous = DataFrame(A = 1.0:5.0, B = 6.0:10.0)
df_counts = DataFrame(A = 1:5, B = 6:10)
df_class = DataFrame(A = ["a", "b", "a", "b", "a"], B = ["a", "b", "c", "a", "b"])

@testset "Scitype coercion tests" begin
    BioFindr.coerce_scitypes!(df_continuous, ScientificTypes.Continuous)
    @test all([all(scitype.(col) .== ScientificTypes.Continuous) for col in eachcol(df_continuous)])

    BioFindr.coerce_scitypes!(df_counts, ScientificTypes.Count)
    @test all([all(scitype.(col) .== ScientificTypes.Count) for col in eachcol(df_counts)])

    BioFindr.coerce_scitypes!(df_class, ScientificTypes.Multiclass)
    @test all([all(scitype.(col) .<: ScientificTypes.Finite) for col in eachcol(df_class)])

    BioFindr.coerce_scitypes!(df_class, ScientificTypes.OrderedFactor)
    @test all([all(scitype.(col) .<: ScientificTypes.Finite) for col in eachcol(df_class)])
end
