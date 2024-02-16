"""
    splitAB(dP)

Split a DataFrame of [`findr`](@ref) results (list of edges) into two DataFrames, one with all `Source` to `Source` edges, and the other with all `Source` to `non-Source` edges. 
"""
function splitAB(dP)
    source_names = unique(dP.Source);
    dP1 = filter(row -> row.Target ∈ source_names, dP)
    dP2 = filter(row -> row.Target ∉ source_names, dP) 
    return dP1, dP2
end
