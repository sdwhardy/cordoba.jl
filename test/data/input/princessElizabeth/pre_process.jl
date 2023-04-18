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
            #if (_row0[:node]==2.0 && _row1[:node] <=25.0)
            #if (_row0[:node]>2.0 && _row0[:node]<=25.0)
            if (_row0[:node]>25.0)
            #if (_row0[:node]>35.0)         
            #if (issubset([_row0[:node]],[38.0,83.0,104.0,59.0,80.0,91.0,108.0,79.0]))
            #if (issubset([_row0[:node]],[119.0,133.0,95.0,123.0,106.0,74.0,130.0,138.0,43.0,59.0,54.0,79.0]))
            push!(complete_graph, [_row0[:node],_row1[:node],Geodesy.euclidean_distance(from_xy, to_xy, 32, true, Geodesy.wgs84)/1000*1.25])
        end;end
    end
end

reduced_graph=_CBD.remove_longest_Xpercent(80,complete_graph)
println(reduced_graph)

