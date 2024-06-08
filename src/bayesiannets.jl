"""
    dagfindr!(dP::T; method="greedy edges") where T<:AbstractDataFrame

Convert a DataFrame `dP` of [`findr`](@ref) results (list of edges) to a directed acyclic graph (DAG) using the specified `method`. The output is a directed graph
represented as a [`SimpleDiGraph`](https://juliagraphs.org/Graphs.jl/stable/core_functions/simplegraphs/#Graphs.SimpleGraphs.SimpleDiGraph) from the [`Graphs`](https://github.com/JuliaGraphs/Graphs.jl) package. The `method` can be any of

- `"greedy edges"` (default), see [`dagfindr_greedy_edges!`](@ref).
- `"heuristic sort"`, see [`dagfindr_heuristic_sort!`](@ref).
- `"greedy insertion"`, see [`dagfindr_greedy_insertion!`](@ref).
"""
function dagfindr!(dP::T; method="greedy edges") where T<:AbstractDataFrame
    if method == "greedy edges"
        return dagfindr_greedy_edges!(dP)
    elseif method == "heuristic sort"
        return dagfindr_heuristic_sort!(dP)
    elseif method == "greedy insertion"
        return dagfindr_greedy_insertion!(dP)
    else
        error("Method not implemented")
    end
end

"""
    dagfindr_greedy_edges!(dP::T) where T<:AbstractDataFrame

Convert a DataFrame of `dP` of [`findr`](@ref) results (list of edges) to a directed acyclic graph (DAG) represented as a [`SimpleDiGraph`](https://juliagraphs.org/Graphs.jl/stable/core_functions/simplegraphs/#Graphs.SimpleGraphs.SimpleDiGraph) from the [`Graphs`](https://github.com/JuliaGraphs/Graphs.jl) package. This function implements the method of [Wang et al. (2019)](https://doi.org/10.3389/fgene.2019.01196) where edges are added one-by-one in decreasing order of probability, and only if they do not create a cycle in the graph, using the incremental cycle detection algorithm from the [`Graphs`](https://github.com/JuliaGraphs/Graphs.jl) package.
"""
function dagfindr_greedy_edges!(dP::T) where T<:AbstractDataFrame
    # Names of the source vertices
    source_names = unique(dP.Source);
    
    # Sort the input edges in decreasing order of probability
    sort!(dP, :"qvalue")

    # Add columns with vertex numbers
    name2idx = names_to_index!(dP)

    # Add a column to store if an edge is added or not
    dP.inDAG_greedy_edges = trues(nrow(dP))

    # Create a directed graph with the right number of vertices
    G = SimpleDiGraph(length(name2idx))

    # Create an incremental cycle checker for G
    ict = IncrementalCycleTracker(G)
    
    # Add edges one by one, in decreasing order of probability
    for row in eachrow(dP)
        if row.Target ∈ source_names
            ae = add_edge_checked!(ict, row.Source_idx, row.Target_idx)
            if !ae
                row.inDAG_greedy_edges = false
            end
        else
            add_edge!(G, row.Source_idx, row.Target_idx)            
        end
    end

    return G, name2idx
end

"""
    dagfindr_heuristic_sort!(dP::T) where T<:AbstractDataFrame

Convert a DataFrame of `dP` of [`findr`](@ref) results (list of edges) to a directed acyclic graph (DAG) represented as a [`SimpleDiGraph`](https://juliagraphs.org/Graphs.jl/stable/core_functions/simplegraphs/#Graphs.SimpleGraphs.SimpleDiGraph) from the [`Graphs`](https://github.com/JuliaGraphs/Graphs.jl) package. This function implements the heuristic sort method of Stoop et al. (2023) where vertices are sorted by their ratio of out-degree to in-degree, and edges are added only if their source vertex precedes their target vertex in the sorted list. The output is a directed graph and a dictionary to map vertex names to numbers.
"""
function dagfindr_heuristic_sort!(dP::T; epsilon=0.01) where T<:AbstractDataFrame
    # Names of the source vertices
    source_names = unique(dP.Source);
    
    # Sort the input edges in decreasing order of probability
    sort!(dP, :"qvalue")

    # Add columns with vertex numbers
    name2idx = names_to_index!(dP)

    # Add a column to store if an edge is added or not
    dP.inDAG_heuristic_sort = trues(nrow(dP))

    # Create a directed graph with the right number of vertices
    G = SimpleDiGraph(length(name2idx))

    # Create a dictionary to store the out-degree of each vertex
    outdeg = Dict{Int,Float64}()

    # Create a dictionary to store the in-degree of each vertex
    indeg = Dict{Int,Float64}()

    # Create a dictionary to store the ratio of out-degree to in-degree of each vertex
    ratio = Dict{Int,Float64}()

    # Compute the out-degree and in-degree of each vertex
    for row in eachrow(dP)
        if row.Target ∈ source_names
            outdeg[row.Source_idx] = get(outdeg, row.Source_idx, 0) + row.Probability
            indeg[row.Target_idx] = get(indeg, row.Target_idx, 0) + row.Probability
        end
    end

    # Compute the ratio of out-degree to in-degree of each vertex
    for v in keys(name2idx)
        ratio[name2idx[v]] = get(outdeg, name2idx[v], 0) / (get(indeg, name2idx[v], 0) + epsilon)
    end

    # Sort the vertices by their ratio of out-degree to in-degree
    sorted_vertices = sort(collect(keys(ratio)), by=x->ratio[x], rev=true)

    # Add edges one by one, in decreasing order of probability
    for row in eachrow(dP)
        if row.Target ∈ source_names
            #if findfirst(x->x==row.Source_idx, sorted_vertices) < findfirst(x->x==row.Target_idx, sorted_vertices)
            if get(ratio, row.Source_idx, 0) > get(ratio, row.Target_idx, 0)
                add_edge!(G, row.Source_idx, row.Target_idx)
            else
                row.inDAG_heuristic_sort = false
            end
        else
            add_edge!(G, row.Source_idx, row.Target_idx)            
        end
    end

    return G, name2idx
end

"""
    dagfindr_greedy_insertion!(dP::T) where T<:AbstractDataFrame

Convert a DataFrame of `dP` of [`findr`](@ref) results (list of edges) to a directed acyclic graph (DAG) represented as a [`SimpleDiGraph`](https://juliagraphs.org/Graphs.jl/stable/core_functions/simplegraphs/#Graphs.SimpleGraphs.SimpleDiGraph) from the [`Graphs`](https://github.com/JuliaGraphs/Graphs.jl) package. This function implements the greedy insertion method of Stoop et al. (2023) where vertices are sorted iteratively by inserting vertices in the position in the current ordering yields the maximum possible gain of edge weights, where the gain is counted as the difference between the sum of new edges weight included and the sum of old edge weights lost, where edges are counted only if their source vertex precedes their target vertex in the ordering. The output is a directed graph and a dictionary to map vertex names to numbers.
"""
function dagfindr_greedy_insertion!(dP::T) where T<:AbstractDataFrame
    # Names of the source vertices
    source_names = unique(dP.Source);
    
    # Sort the input edges in decreasing order of probability
    sort!(dP, :"qvalue")

    # Add columns with vertex numbers
    name2idx = names_to_index!(dP)

    # Add a column to store if an edge is added or not
    dP.inDAG_greedy_insertion = trues(nrow(dP))

    # Create a directed graph with the right number of vertices
    G = SimpleDiGraph(length(name2idx))

    # Initialize the ordering of the vertices
    n = length(name2idx)
    sorted_vertices = Int64[1:n;]

    # Get the edge weights
    weights = edge_weights(dP)

    # Perform greedy insertions
    greedy_insertions!(sorted_vertices, weights)

    # Get rank of each vertex in the sorted list
    vertex_rank = invperm(sorted_vertices)

    # Add edges one by one, in decreasing order of probability
    for row in eachrow(dP)
        if row.Target ∈ source_names
            if vertex_rank[row.Source_idx] < vertex_rank[row.Target_idx]
                add_edge!(G, row.Source_idx, row.Target_idx)
            else
                row.inDAG_greedy_insertion = false
            end
        else
            add_edge!(G, row.Source_idx, row.Target_idx)            
        end
    end

    return G, name2idx
end

"""
    greedy_insertions!(sorted_vertices, weights)

TBW
"""
function greedy_insertions!(sorted_vertices, weights)
    n = length(sorted_vertices)
    updated = true
    while updated
        updated = false
        for x = 1:n
            ymax = x
            ΔEmax = 0
            ΔE = 0
            for y = x-1:-1:1
                ΔE += get(weights, (sorted_vertices[x], sorted_vertices[y]), 0) - get(weights, (sorted_vertices[y], sorted_vertices[x]), 0)
                if ΔE > ΔEmax
                    ymax = y
                    ΔEmax = ΔE
                end
            end
            ΔE = 0
            for y = x+1:n
                ΔE += get(weights, (sorted_vertices[y], sorted_vertices[x]), 0) - get(weights, (sorted_vertices[x], sorted_vertices[y]), 0)
                if ΔE > ΔEmax
                    ymax = y
                    ΔEmax = ΔE
                end
            end
            if ymax != x
                insert!(sorted_vertices, ymax, popat!(sorted_vertices, x))
                updated = true              
            end
        end
    end
end

"""
    edge_weights(dP::T) where T<:AbstractDataFrame

TBW
"""
function edge_weights(dP::T) where T<:AbstractDataFrame
    weights = Dict{Tuple{Int,Int},Float64}()
    for row in eachrow(dP)
        weights[(row.Source_idx, row.Target_idx)] = row.Probability
    end
    return weights    
end

"""
    names_to_index!(dP::T) where T<:AbstractDataFrame

Add columns with vertex numbers to a DataFrame `dP` of edges. The columns `Source_idx` and `Target_idx` are added to `dP` with the vertex numbers corresponding to the names in the `Source` and `Target` columns, respectively. The function returns a dictionary `name2num` to map vertex names to numbers.
"""
function names_to_index!(dP::T) where T<:AbstractDataFrame
    # Create a dictionary to map names to numbers, BayesNets wants this as a Dict{Symbol,Int}
    vnames = unique(vcat(dP.Source, dP.Target))
    name2idx = Dict(name => i for (i, name) in enumerate(vnames))
    # Add columns with vertex numbers
    dP.Source_idx = map(x -> name2idx[x], dP.Source)
    dP.Target_idx = map(x -> name2idx[x], dP.Target)
    # BayesNets wants the names as Symbols
    symb2idx = Dict(Symbol(name) => i for (name, i) in name2idx)
    return symb2idx
end