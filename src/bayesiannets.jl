
"""
    dagfindr(dP::T) where T<:AbstractDataFrame

Convert a DataFrame of [`findr`](@ref) results (list of edges) to a directed acyclic graph (DAG) represented as a [`SimpleDiGraph`](https://juliagraphs.org/Graphs.jl/stable/core_functions/simplegraphs/#Graphs.SimpleGraphs.SimpleDiGraph) from the [`Graphs`](https://github.com/JuliaGraphs/Graphs.jl) package. Edges are added one-by-one in decreasing order of probability, and only if they do not create a cycle in the graph, using the incremental cycle detection algorithm from the [`Graphs`](https://github.com/JuliaGraphs/Graphs.jl) package.
"""

function dagfindr(dP::T) where T<:AbstractDataFrame
    # Add columns with vertex numbers
    name2num = names2nums!(dP)

    # Split and sort dP into two DataFrames, one with all Source to Source edges, and the other with all Source to non-Source edges
    dP1, dP2 = splitAB(dP)
    
    # Create a directed graph with the right number of vertices
    G = SimpleDiGraph(length(name2num))

    # Add all Source to Source edges, check for cycles
    ict = IncrementalCycleTracker(G)
    for row in eachrow(dP1)
        add_edge_checked!(ict, row.Source_idx, row.Target_idx)
    end

    # Add all Source to non-Source edges, no need to check for cycles
    for row in eachrow(dP2)
        add_edge!(G, row.Source_idx, row.Target_idx)
    end

    return G
end

"""
    splitAB(dP)

Split and sort a DataFrame of [`findr`](@ref) results (list of edges) into two DataFrames, one with all `Source` to `Source` edges, and the other with all `Source` to `non-Source` edges. 
"""
function splitAB(dP::T) where T<:AbstractDataFrame
    source_names = unique(dP.Source);
    dP1 = filter(row -> row.Target ∈ source_names, dP)
    sort!(dP1, :"Probability", rev=true) # sort in decreasing order of probability
    dP2 = filter(row -> row.Target ∉ source_names, dP)
    sort!(dP2, :"Probability", rev=true) # sort in decreasing order of probability
    return dP1, dP2
end

function names2nums!(dP::T) where T<:AbstractDataFrame
    # Create a dictionary to map names to numbers
    vnames = unique(vcat(dP.Source, dP.Target))
    name2num = Dict(name => i for (i, name) in enumerate(vnames))
    # Add columns with vertex numbers
    dP.Source_idx = map(x -> name2num[x], dP.Source)
    dP.Target_idx = map(x -> name2num[x], dP.Target)
    return name2num
end