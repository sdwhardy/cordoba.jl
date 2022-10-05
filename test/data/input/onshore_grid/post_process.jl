################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

results_nodal=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\nodal_market_k6_NL_0gap.jld2")
results_allhm=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_allhm_k6_VOLL5000_wonshore.jld2")
results_567=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_567_k6_VOLL5000b_wonshore.jld2")

s_nodal, result_mip_nodal, data_nodal, mn_data_nodal=_CBD.summarize_in_s(results_nodal);
s_hm, result_mip_hm, data_hm, mn_data_hm=_CBD.summarize_zonal_in_s(results_allhm);
s_567, result_mip_567, data_567, mn_data_567=_CBD.summarize_zonal_in_s(results_567);

s_hm["map_gen_types"]["countries"]["BE"]
dk_gen_load["ID"][!,Symbol.("56")]
result_mip_hm["solution"]["nw"]["1"]["gen"]
_CBD.print_solution_wcost_data(result_mip_nodal, s_nodal, data_nodal)#-856896.0245340846
_CBD.print_solution_wcost_data(result_mip_567, s_567, data_567)#-856896.0245340846
_CBD.print_solution_wcost_data(result_mip_hm, s_hm, data_hm)#-856896.0245340846
_CBD.print_table_summary(s_nodal)
_CBD.print_table_summary(s_hm)
_CBD.print_table_summary(s_24)
_CBD.print_table_summary(s_34)

_CBD.topology_map(s_567,1.48)

result_mip_n2z, data_n2z, mn_data_n2z, s_n2z, result_mip_hm_prices_n2z=_CBD.nodal2zonal(results_nodal["s"],results_nodal["result_mip"],[[2,5],[3,6],[4,7]]);
results_n2z=Dict("result_mip"=>result_mip_n2z,"data"=>data_n2z, "mn_data"=>mn_data_n2z, "s"=>s_n2z, "result_mip_hm_prices"=>result_mip_hm_prices_n2z)
s_n2z, result_mip_n2z, data_n2z, mn_data_n2z=_CBD.summarize_zonal_in_s(results_n2z);
_CBD.print_table_summary(s_n2z)

result_mip_z2n, data_z2n, mn_data_z2n, s_z2n, result_mip_hm_prices_z2n=_CBD.nodal2zonal(results_allhm["s"],results_allhm["result_mip"],[]);
results_z2n=Dict("result_mip"=>result_mip_z2n,"data"=>data_z2n, "mn_data"=>mn_data_z2n, "s"=>s_z2n, "result_mip_hm_prices"=>result_mip_hm_prices_z2n)
s_z2n, result_mip_z2n, data_z2n, mn_data_z2n=_CBD.summarize_zonal_in_s(results_z2n);
_CBD.print_table_summary(s_z2n)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_allhm_k6_VOLL5000_wonshoreZ2N.jld2",results_z2z)

############################
s=results_allhm["s"]
result_mip=results_allhm["result_mip"]
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
s["home_market"]=[[2,5],[3,6],[4,7]]
s["rebalancing"]=true
s["relax_problem"]=true
s["output"]["duals"]=true
mn_data, data, s = _CBD.data_update(s,result_mip);#Build data structure for given options
mn_data, s = _CBD.set_rebalancing_grid(result_mip,mn_data,s);
s, mn_data= _CBD.remove_integers(result_mip,mn_data,data,s);
result_mip_hm_prices = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem

results=Dict("result_mip"=>result_mip_hm_prices,"data"=>data, "mn_data"=>mn_data, "s"=>s, "result_mip_hm_prices"=>result_mip_hm_prices)
dk_gen_load=_CBD.InitialD_FinalD_ReDispatch(results)
df_bus=_CBD.gen_load_values(mn_data["nw"],"gen_bus")
df_bus=df_bus[!,Symbol.(names(dk_gen_load["FD"]))]
#get dataframes of NPV/Orig clearing prices per node
dk_price=_CBD.bus_values(df_bus,results["result_mip"]["solution"]["nw"],results["s"])
dk_price["Orig"][1:144,Symbol.(names(dk_price["Orig"])[25:43])]
############################
delete!(mn_data["scenario"]["1"]["433"],"6")
mn_data_test=deepcopy(mn_data)
for k=433:1:2592
delete!(mn_data["nw"],string(k))
end
s
############################

result_mip_z2n, data_z2n, mn_data_z2n, s_z2n, result_mip_hm_prices_z2n=_CBD.nodal2zonal(results_allhm["s"],results_allhm["result_mip"],[[2,5],[3,6],[4,7]]);
results_z2n=Dict("result_mip"=>result_mip_z2n,"data"=>data_z2n, "mn_data"=>mn_data_z2n, "s"=>s_z2n, "result_mip_hm_prices"=>result_mip_hm_prices_z2n)
s_z2n, result_mip_z2n, data_z2n, mn_data_z2n=_CBD.summarize_zonal_in_s(results_z2n);
_CBD.print_table_summary(s_z2n)
results_allhm["s"]["xd"]["gen"]["221"]
s
s=s_hm;mn_data=mn_data_hm;data=data_hm
_CBD.print_table_summary(s)

_CBD.topology_map(s)

_CBD.plot_cumulative_production_all_scenarios_allWF(s, mn_data)

_CBD.plot_cumulative_income_tl_all_scenarios(s,data)

_CBD.print_solution_wcost_data(result_mip, s, data)

for (k,c) in s["income_summary"]["tso"]["totals"]["ac"];println("AC("*k*"):"*string(c));end

for (k,c) in s["income_summary"]["tso"]["totals"]["dc"];println("DC("*k*"):"*string(c));end

nodal_prices=s_hm["income_summary"]["onshore"]["1"]["2"]["price_not_npv"]
[println(p) for p in nodal_prices] 
hm_prices=results_allhm["s"]["income_summary"]["owpp"]["6"]["2"]["price_not_npv"]
print(keys(results_allhm["s"]))
[println(p) for p in hm_prices]

[println(string(n)*": "*string(results_allhm["result_mip"]["solution"]["nw"]["1"]["gen"][string(n)]["pg"])) for n=keys(results_allhm["result_mip"]["solution"]["nw"]["1"]["gen"])]
results["s"]["map_gen_types"]["type"][224]
results["s"]["map_gen_types"]["countries"]["BE"]
results_allhm["result_mip_hm_prices"]["solution"]["nw"]["100"]["gen"]["225"]
println(dk_gen_load["ID"][1:244,Symbol.(names(dk_price["Orig"])[25:43])])
println(dk_gen_load["ID"][1:244,Symbol.(names(dk_price["Orig"])[85:94])])
println(dk_gen_load["ID"][1:244,Symbol.(names(dk_price["Orig"])[82:82])])
names(dk_price["Orig"])[1:10]
dk_price["Orig"][1:288,Symbol.(names(dk_price["Orig"])[25:43])]
[println(results_allhm["result_mip_hm_prices"]["solution"]["nw"][string(i)]["branchdc"]["7"]["pt"]) for i=1:1:244]
results_allhm["data"]["branch"]["1"]["pt"]
#4 1-2 (10MVA) 10
#7 2-3 (10MVA) -4
println(keys(results_allhm["data"]["scenario"]["1"]))
###############################
######################################################
#get dictionary of ID, FD, RDispatch
results=results_allhm
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


time_series=s["income_summary"]["onshore"]["1"]["3"]["price_not_npv"]
_CBD.plot_clearing_price(time_series)=#

function plot_clearing_price(time_series)
    
    clrs=generation_color_map()
    

        #low_rng=minimum(marginal_prices)
        #high_rng=maximum(marginal_prices)
        scatter_vec_gen=[
            
            PlotlyJS.scatter(
                y=time_series,
                mode="lines", name="BE",
                line=PlotlyJS.attr(width=2, color="red")
            ) ]
            scatter_vec_gen=[
            
            PlotlyJS.scatter(
                y=time_series,
                mode="lines", name="BE",
                line=PlotlyJS.attr(width=2, color="red")
            ) ]
            scatter_vec_gen=[
            
            PlotlyJS.scatter(
                y=time_series,
                mode="lines", name="BE",
                line=PlotlyJS.attr(width=2, color="red")
            ) ]
        PlotlyJS.plot(
            scatter_vec_gen, PlotlyJS.Layout(font_size=35,xaxis_range=(0, 432),yaxis_title="€/MWh",xaxis_title="time steps"))
    end