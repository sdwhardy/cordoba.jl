################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections
##################### File parameters #################################
s = Dict(
"rt_ex"=>pwd()*"\\test\\data\\input\\UK_BE\\",#folder path
"scenario_data_file"=>"C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBE58_2021.jld2",
################# temperal parameters #################
"test"=>false,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>2,
"scenario_names"=>["NT"],#["NT","DE","GA"]
"k"=>365,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k Available: 2, 5, 10, 50, 100
"res_years"=>["2020"],#Options: ["2012","2013","2014","2015","2016"]
"scenario_years"=>["2020"],#Options: ["2020","2030","2040"]
"dr"=>0.04,#discount rate
"yearly_investment"=>100000,
################ electrical parameters ################
"AC"=>"1",#0=false, 1=true
"owpp_mva"=>[5800],#mva of wf in MVA
"conv_lim_onshore"=>5800,#Max Converter size in MVA
"conv_lim_offshore"=>5800,#Max Converter size in MVA
"strg_lim_offshore"=>0.2,
"strg_lim_onshore"=>10,
"candidate_ics_ac"=>[1/10],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1],#DC Candidate Cable sizes (fraction of full MVA)
################## optimization/solver setup options ###################
"output" => Dict("branch_flows" => false),
"eps"=>0.0001,#admm residual (100kW)
"beta"=>5.5,
"relax_problem" => false,
"conv_losses_mp" => true,
"process_data_internally" => false,
"corridor_limit" => true)
########################################################################
#0.0066 - branch
##################################### HM market 
################## Run MIP Formulation ###################
#NOTE only very basic intuitive check passed on functions wgen_type
s["home_market"]=[]
mn_data, data, s = _CBD.data_setup_nodal(s);
s["xd"]["gen"]["58"]["pmax"]=s["xd"]["gen"]["58"]["pmax"].*0
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap"=>0.9e-3)#select solver
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
result_mip["solution"]["nw"]["1"]["branchdc_ne"]
result_mip["solution"]["nw"]["1"]["convdc"]
result_mip["solution"]["nw"]["1"]["convdc_ne"]
result_mip["solution"]["nw"]["1"]["storage"]
######################################################
#get dictionary of ID, FD, RDispatch
results=results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s)
dk_gen_load=_CBD.InitialD_FinalD_ReDispatch(results)
dk_gen_load["FD"]=_CBD.rename_gen_df_columns(results["s"]["map_gen_types"]["type"],dk_gen_load["FD"])
data=FileIO.load(s["scenario_data_file"])
dk_gen_load["FD"]=hcat(dk_gen_load["FD"],data["Demand"]["Base"]["2020"])
dk_gen_load["FD"]=hcat(dk_gen_load["FD"],data["Generation"]["RES"]["Offshore Wind"]["BE"]["2020"][!,:BE_MWh])

dk_gen_load["FD"][!,:BE_demand_total]=sum.(eachrow(dk_gen_load["FD"][!,[Symbol(i) for i=112:1:121]]))
dk_gen_load["FD"][!,:UK_demand_total]=sum.(eachrow(dk_gen_load["FD"][!,[Symbol(i) for i=122:1:131]]))
dk_gen_load["FD"][!,[:UK_MWh,:UK_demand_total]]
dk_gen_load["FD"][!,[:BE_MWh,:BE_demand_total]]
dk_gen_load["FD"]=dk_gen_load["FD"][!,Not([Symbol(i) for i=112:1:131])]
dk_gen_load["FD"][!,:Generation_total]=sum.(eachrow(dk_gen_load["FD"][!,Not([:time_stamp,:BE_MWh,:UK_MWh,:BE_demand_total,:UK_demand_total])]))
dk_gen_load["FD"][!,:Demand_total]=sum.(eachrow(dk_gen_load["FD"][!,[:BE_demand_total,:UK_demand_total]]))
dk_gen_load["FD"][!,[:Generation_total,:Demand_total]]
dk_gen_load["FD"][!,:x1]=dk_gen_load["FD"][!,:x1].*s["owpp_mva"]/100
dk_gen_load["FD"][!,[Symbol("Offshore Wind_2"),:x1]]
ordered_solution=sort!(OrderedCollections.OrderedDict(results["result_mip"]["solution"]["nw"]), by=x->parse(Int64,x))

dk_gen_load["FD"][!,:cable_1]=[ordered_solution[string(i)]["branchdc_ne"]["1"]["pt"] for i=1:1:8784]
dk_gen_load["FD"][!,:cable_2]=[ordered_solution[string(i)]["branchdc_ne"]["2"]["pt"] for i=1:1:8784]
dk_gen_load["FD"][!,Not(Symbol("time_stamp"))]=dk_gen_load["FD"][!,Not(Symbol("time_stamp"))].*100
CSV.write("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\YUSO_data\\DA_auction_2020_G2.csv", dk_gen_load["FD"]);

dk_gen_load["FD"][!,:cable_1]=[ordered_solution[string(i)]["branchdc_ne"]["1"]["pt"] for i=1:1:8760]
dk_gen_load["FD"][!,:cable_2]=[ordered_solution[string(i)]["branchdc_ne"]["2"]["pt"] for i=1:1:8760]
dk_gen_load["FD"][!,:time_stamp]=dk_gen_load["FD"][!,:time_stamp].+Year(1)
dk_gen_load["FD"][!,Not(Symbol("time_stamp"))]=dk_gen_load["FD"][!,Not(Symbol("time_stamp"))].*100
CSV.write("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\YUSO_data\\DA_auction_2021_G2.csv", dk_gen_load["FD"])


print(dk_gen_load["FD"][!,Symbol("Offshore Wind_2")])



sort!()[]["branchdc_ne"]["1"]
results["result_mip"]["solution"]["nw"]
print(names(dk_gen_load["FD"]))
print(names(dk_gen_load["FD"]))
sum.(eachrow(dk_gen_load["FD"][!,[Symbol(i) for i=112:1:121]]))
for c in names(dk_gen_load["FD"])
    if issubset([c],first.(results["s"]["map_gen_types"]["type"]))
    else
        push!(results["s"]["map_gen_types"]["type"],(c,c))
    end
end
_syms=[Symbol(last(results["s"]["map_gen_types"]["type"][parse(Int64,c)])) for c in names(dk_gen_load["FD"])]

rename!(dk_gen_load["FD"],_syms, makeunique=true)


intersect(names(dk_gen_load["FD"]),first.())
rename!(data1_date_time_index, Symbol.(colnames))
#Get Dataframe of the bus numbers of each generator
df_bus=_CBD.gen_load_values(results["mn_data"]["nw"],"gen_bus")
df_bus=df_bus[!,Symbol.(names(dk_gen_load["FD"]))]
#get dataframes of NPV/Orig clearing prices per node
dk_price=_CBD.bus_values(df_bus,results["result_mip"]["solution"]["nw"],results["s"])
#Get Dataframe of generator NPV hourly values 
push!(dk_price,"GENS"=>_CBD.gen_bid_prices(results["s"]["xd"]["gen"],Symbol.(names(dk_gen_load["RD"]))))
#seperate the up and down regulation
pos,neg=_CBD.decompose_re_dispatch(dk_gen_load["RD"])
#calculate the Re dispatch cost
a=pos.*dk_price["GENS"]
b=neg.*(dk_price["NPV"].-dk_price["GENS"])
rbc=(sum(sum.(eachcol(a)))+sum(sum.(eachcol(b))))/6#191272.64882850088
######################################################


#=print_solution_wcost_data(result_mip, s, data)
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
s["rebalancing"]=true
s["relax_problem"]=true
s["output"]["duals"]=true
mn_data, data, s = data_update(s,result_mip);#Build data structure for given options
mn_data, s = set_rebalancing_grid(result_mip,mn_data,s);
s, mn_data= remove_integers(result_mip,mn_data,data,s);
result_mip =  cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
=#   

s["xd"]["gen"]["5"]["pmax"]
10*8.717
10*30.006



