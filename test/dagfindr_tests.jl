Random.seed!(9999)
nA = 2
nB = 50
fB = 0.1
ns = 80
ng = 3
maf = 0.3
bGA = 1.0
bAB = 1.0

G_dag, XA_dag, XB_dag, _ = BioFindr.generate_test_data(nA, nB, fB, ns, ng, maf, bGA, bAB, true)

# Build a DataFrame of edges from findr as input to dagfindr!
pairGX = hcat(1:nA, 1:nA)
PP = BioFindr.findr_matrix(XB_dag, XA_dag, G_dag, pairGX; combination="IV")
dP_dag = BioFindr.stackprobs(PP, names(DataFrame(XA_dag, :auto)), names(DataFrame(XB_dag, :auto)))
BioFindr.globalfdr!(dP_dag; FDR=1.0)

@testset "dagfindr! greedy edges" begin
    dP_copy = copy(dP_dag)
    G, name2idx = dagfindr!(dP_copy; method="greedy edges")
    @test G isa Graphs.SimpleDiGraph
    @test length(name2idx) > 0
    @test !Graphs.is_cyclic(G)
    # result column should be present
    @test "inDAG_greedy_edges" ∈ names(dP_copy)
end

@testset "dagfindr! heuristic sort" begin
    dP_copy = copy(dP_dag)
    G, name2idx = dagfindr!(dP_copy; method="heuristic sort")
    @test G isa Graphs.SimpleDiGraph
    @test !Graphs.is_cyclic(G)
    @test "inDAG_heuristic_sort" ∈ names(dP_copy)
end

@testset "dagfindr! greedy insertion" begin
    dP_copy = copy(dP_dag)
    G, name2idx = dagfindr!(dP_copy; method="greedy insertion")
    @test G isa Graphs.SimpleDiGraph
    @test !Graphs.is_cyclic(G)
    @test "inDAG_greedy_insertion" ∈ names(dP_copy)
end

@testset "dagfindr! unknown method error" begin
    dP_copy = copy(dP_dag)
    @test_throws ErrorException dagfindr!(dP_copy; method="unknown")
end
