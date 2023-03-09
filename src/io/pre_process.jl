function remove_longest_Xpercent(X,complete_graph)
    DataFrames.sort!(complete_graph,[:length])
    full_length=Int64(length(complete_graph[!,:length]))
    Xpercent=Int64(ceil(X/100*length(complete_graph[!,:length])))
    reduced_graph=complete_graph[1:full_length-Xpercent,:]
    DataFrames.sort!(reduced_graph,[:from,:to])
    return reduced_graph
end