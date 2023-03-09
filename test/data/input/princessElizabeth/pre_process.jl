using DataFrames, FileIO, CSV, Dates, OrderedCollections, PlotlyJS, XLSX, Geodesy
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
nodes_df = DataFrames.DataFrame(XLSX.readtable(s["rt_ex"]*"input.xlsx", "node_generation")...)
filter!(:type=>x->x==0,nodes_df)
complete_graph=DataFrames.DataFrame(:from=>[],:to=>[],:length=>[])
for (i0,_row0) in enumerate(eachrow(nodes_df))
    for (i1,_row1) in enumerate(eachrow(nodes_df))
        if (i0<i1)
            from_xy=_CBD.utm_gps2xy((_row0[:lat],_row0[:long]))
            to_xy=_CBD.utm_gps2xy((_row1[:lat],_row1[:long]))
            push!(complete_graph, [_row0[:node],_row1[:node],Geodesy.euclidean_distance(from_xy, to_xy, 32, true, Geodesy.wgs84)/1000*1.25])
        end
    end
end

reduced_graph=_CBD.remove_longest_Xpercent(85,complete_graph)
println(reduced_graph)

