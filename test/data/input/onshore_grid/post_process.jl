################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

######################### nOBZ to HMD
results=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\nodal_market_NORTH_SEA.jld2")
result_mip_n2z, data_n2z, mn_data_n2z, s_n2z, result_mip_hm_prices_n2z=_CBD.nodal2zonal(results["s"],results["result_mip"],[[4,11],[5,10],[6,12],[1,8,13],[3,9]]);
results_n2z=Dict("result_mip"=>result_mip_n2z,"data"=>data_n2z, "mn_data"=>mn_data_n2z, "s"=>s_n2z, "result_mip_hm_prices"=>result_mip_hm_prices_n2z)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\n2HMD_results_NORTH_SEA_0gap.jld2",results_n2z)

######################### nOBZ to zOBZ
results=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\nodal_market_NORTH_SEA.jld2")
result_mip_n2z, data_n2z, mn_data_n2z, s_n2z, result_mip_hm_prices_n2z=_CBD.nodal2zonal(results["s"],results["result_mip"],[[11,10,12,13,9]]);
results_n2z=Dict("result_mip"=>result_mip_n2z,"data"=>data_n2z, "mn_data"=>mn_data_n2z, "s"=>s_n2z, "result_mip_hm_prices"=>result_mip_hm_prices_n2z)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\n2zOBZ_results_NORTH_SEA_0gap.jld2",results_n2z)

######################### HMD to nOBZ
results=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\zonalHM_results_NORTH_SEA.jld2")
result_mip_n2z, data_n2z, mn_data_n2z, s_n2z=_CBD.zonal2nodal(results["s"],results["result_mip"]);
results_n2z=Dict("result_mip"=>result_mip_n2z,"data"=>data_n2z, "mn_data"=>mn_data_n2z, "s"=>s_n2z)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\HMD2nOBZ_results_NORTH_SEA_0gap.jld2",results_n2z)

######################### zOBZ to nOBZ
results=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\zonalOBZ_results_NORTH_SEA.jld2")
result_mip_n2z, data_n2z, mn_data_n2z, s_n2z=_CBD.zonal2nodal(results["s"],results["result_mip"]);
results_n2z=Dict("result_mip"=>result_mip_n2z,"data"=>data_n2z, "mn_data"=>mn_data_n2z, "s"=>s_n2z)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\zOBZ2nOBZ_results_NORTH_SEA_0gap.jld2",results_n2z)

results_nOBZ=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\nodal_market_NORTH_SEA.jld2")
results_HM=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\zonalHM_results_NORTH_SEA.jld2")
results_zOBZ=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\zonalOBZ_results_NORTH_SEA.jld2")
results_HMD2OBZ=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\HMD2nOBZ_results_NORTH_SEA_0gap.jld2")
results_zOBZ2nOBZ=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\zOBZ2nOBZ_results_NORTH_SEA_0gap.jld2")

results=results_nOBZ
results=results_HM
results=results_zOBZ
results=results_HMD2OBZ
results=results_zOBZ2nOBZ
results=results_n2z
push!(results["s"]["map_gen_types"]["markets"][1],"UK")

######################### Nodal summary
s_nodal, result_mip_nodal, data_nodal, mn_data_nodal=_CBD.summarize_in_s(results);
_CBD.print_solution_wcost_data(result_mip_nodal, s_nodal, data_nodal)#-856896.0245340846
_CBD.print_table_summary(s_nodal)

######################### Zonal summary
results["s"]["cost_summary"]=_CBD.print_solution_wcost_data(results["result_mip"], results["s"], results["data"])
s_z, result_mip_z, data_z, mn_data_z=_CBD.summarize_zonal_in_s(results);
_CBD.print_solution_wcost_data(result_mip_z, s_z, data_z)#-856896.0245340846
_CBD.print_table_summary(s_z)


results_n2z["s"]["cost_summary"]=_CBD.print_solution_wcost_data(results_n2z["result_mip"], results_n2z["s"], results_n2z["data"])
push!(results_n2z["s"]["map_gen_types"]["markets"][1],"UK")
s_n2z, result_mip_n2z, data_n2z, mn_data_n2z=_CBD.summarize_zonal_in_s(results_n2z);
_CBD.print_table_summary(s_n2z)

result_mip_z2n, data_z2n, mn_data_z2n, s_z2n, result_mip_hm_prices_z2n=_CBD.nodal2zonal(results_allhm["s"],results_allhm["result_mip"],[]);
results_z2n=Dict("result_mip"=>result_mip_z2n,"data"=>data_z2n, "mn_data"=>mn_data_z2n, "s"=>s_z2n, "result_mip_hm_prices"=>result_mip_hm_prices_z2n)
s_z2n, result_mip_z2n, data_z2n, mn_data_z2n=_CBD.summarize_zonal_in_s(results_z2n);
_CBD.print_table_summary(s_z2n)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\zonal_results_allhm_k6_VOLL5000_wonshoreZ2N.jld2",results_z2z)












results["result_mip"]["solution"]["nw"]["33"]["bus"]["11"]

_CBD.problemOUTPUT_map(results, 0.1, 0.1,1, 1)
######################################################
#get dictionary of ID, FD, RDispatch
dk_gen_load=_CBD.InitialD_FinalD_ReDispatch(results)
#Get Dataframe of the bus numbers of each generator
df_bus=_CBD.gen_load_values(results["mn_data"]["nw"],"gen_bus")
df_bus=df_bus[!,Symbol.(names(dk_gen_load["FD"]))]
#get dataframes of NPV/Orig clearing prices per node
dk_price=_CBD.bus_values(df_bus,results["result_mip"]["solution"]["nw"],results["s"])
#Get Dataframe of generator NPV hourly values 
push!(dk_price,"GENS"=>_CBD.gen_bid_prices(results["s"]["xd"]["gen"],Symbol.(names(dk_gen_load["FD"]))))
#calculate the final dispatch cost
a=dk_gen_load["ID"].*dk_price["NPV"]
push!(a,sum.(eachcol(a)))
a=_CBD.rename_gen_df_columns(results["s"]["map_gen_types"]["type"],a)
sum(a[end,Not([Symbol(string(_name)) for _name = 167:1:196])])/6
#seperate the up and down regulation
pos,neg=_CBD.decompose_re_dispatch(dk_gen_load["RD"])
#calculate the Re dispatch cost
b=pos.*dk_price["GENS"]
c=neg.*(dk_price["NPV"].-dk_price["GENS"])
rbc=(sum(sum.(eachcol(b)))+sum(sum.(eachcol(c))))/6#191272.64882850088
b=_CBD.rename_gen_df_columns(results["s"]["map_gen_types"]["type"],b)

d=(sum.(eachcol(b)).+sum.(eachcol(c)))./6
for (_c,_v) in enumerate(d);println(names(b)[_c]*":"*string(_v));end

names(b)[argmax(d)]
push!(b,sum.(eachcol(b)))
push!(c,sum.(eachcol(c)))
b=_CBD.rename_gen_df_columns(results["s"]["map_gen_types"]["type"],b)
c=_CBD.rename_gen_df_columns(results["s"]["map_gen_types"]["type"],c)

(sum(b[end,:]+c[end,:]))/6

push!(dk_gen_load["FD"],sum.(eachcol(dk_gen_load["FD"])))
dk_gen_load["FD"]=_CBD.rename_gen_df_columns(results["s"]["map_gen_types"]["type"],dk_gen_load["FD"])
sum(dk_gen_load["FD"][end,[:SLACK,:SLACK_1,:SLACK_2]])/6
print(names(a))

######################################################

results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s)
print_solution_wcost_data(result_mip_z, s_z, data_z)
dk_gen_load=InitialD_FinalD_ReDispatch(results)#generators
#Get Dataframe of the bus numbers of each generator
country="DK"
df_bus=gen_load_values(results["mn_data"]["nw"],"gen_bus")#buses
df_bus[!,Symbol.([k for k=s["map_gen_types"]["countries"][country] if (issubset([k],names(dk_gen_load["FD"])))])]#gen bus numbers for country
dk_gen_load["FD"][!,Symbol.([k for k=s["map_gen_types"]["countries"][country] if (issubset([k],names(dk_gen_load["FD"])))])]#gen bus values for country
df_bus[!,Symbol.([string(k) for k=s["map_gen_types"]["loads"][country] if (issubset([string(k)],names(dk_gen_load["FD"])))])]#loads bus numbers for country
dk_gen_load["FD"][!,Symbol.([string(k) for k=s["map_gen_types"]["loads"][country] if (issubset([string(k)],names(dk_gen_load["FD"])))])]#loads bus values for country
a=sum.(eachrow(dk_gen_load["FD"][!,Symbol.([k for k=s["map_gen_types"]["countries"][country] if (issubset([k],names(dk_gen_load["FD"])))])]))
b=sum.(eachrow(dk_gen_load["FD"][!,Symbol.([string(k) for k=s["map_gen_types"]["loads"][country] if (issubset([string(k)],names(dk_gen_load["FD"])))])]))
c=sum.(eachrow(dk_gen_load["FD"][!,Symbol.([string(k) for k=s["map_gen_types"]["offshore"][country] if (issubset([string(k)],names(dk_gen_load["FD"])))])]))
a+b+c
print_mn_data(mn_data,s)