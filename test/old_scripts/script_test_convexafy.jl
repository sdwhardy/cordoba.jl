################## loads external packages ##############################
using Ipopt, Juniper, JuMP, Cbc, Gurobi, Cbc, XLSX, DataFrames, Dates, CSV, FileIO
import cordoba; const _CBD = cordoba#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import FlexPlan; const _FP = FlexPlan
import InfrastructureModels; const _IM = InfrastructureModels
using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
include("../test/data/conv_spec.jl")

################### ENTSO-E scenario ####################################
scenario_names=["EU19","EU20","ST19","ST20","DG19","DG20"]
scenario_years=["2020","2030","2040"]
scenario_planning_horizon=30

################### Input data for script_test.jl ########################
owpp_mva=[4000]#mva of wf (in de)
#interconnectors format: (mva,km)
#must be in same order as .m file ie uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3), wf(4)-uk(1), wf(4)-de(2), wf(4)-dk(3)
ics=[(4000,550),(4000,760),(4000,250),(4000,470),(4000,145),(4000,246)];
ic_lim=maximum([_CBD.DC_cbl(ic,1).num*_CBD.DC_cbl(ic,1).elec.mva for ic in first.(ics)])-0.01

#ics=[(2000,145)];
#location of nodes
markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

################## reads .m input file name ######################
casename = "test_convexafy"
file = "./test/data/input/$casename.m"
data = _PM.parse_file(file)#load data in PM format

#################### PowerModelsACDC settings ############################
_PMACDC.process_additional_data!(data)#add extra DC model data

##################### load time series data ##############################
#scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series
#=scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series

##################### Cluster time series data ###########################
k=50
for (sc,yrs_ts) in scenario_data
    for (yr,ts) in yrs_ts
        daily_ts=_CBD.daily_tss(DateTime.(ts["time_stamp"]))
        daily_dk=_CBD.daily_tss(Float64.(ts["Wnd_MWhDK"]))
        daily_de=_CBD.daily_tss(Float64.(ts["EUR_daDE"]))
        daily_uk=_CBD.daily_tss(Float64.(ts["EUR_daUK"]))

        kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_uk.^2).+(daily_de.^2).+(daily_dk.^2))',axis=1), k)
        ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
        filter!(row -> row.time_stamp in ts, scenario_data[sc][yr])
    end
end=#

scenario_data=load("./test/data/input/EUSTDG_TS_k5.jld2")
##################### Find minimum length scenario
ls=[];for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
push!(ls,length(scenario_data[_yr][_sc].time_stamp))
end;end;ls=minimum(ls)

##################### Make all scenarios the same length
for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
end;end


#################### Multi-period input parameters #######################
all_scenario_data,data,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
extradata,data = _CBD.create_profile_sets_mesh(dim, data, all_scenario_data, markets_wfs, infinite_grid, owpp_mva)

# Scale cost data
_CBD.scale_cost_data_cordoba!(data, scenario)
_CBD.scale_cost_data_cordoba!(extradata, scenario)

# Create data dictionary where time series data is included at the right place
mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "per_unit", "source_version"]))

#select solver
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

#optimation settings
scenes_years=[];for (i,(k,_scene)) in enumerate(scenario["sc_names"]);push!(scenes_years,[]);for (j, _yr) in _scene; push!(scenes_years[i],_yr);end;end
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls, "ic_lim"=>ic_lim)

#run optimization
result_conv = _CBD.cordoba_mp_acdctnepopf_convexafy(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s);









result["solution"]["nw"]["1"]


for (k,n) in result["solution"]["nw"]
    println(k*": pg: "*string(n["gen"]["4"]["pg"]))
end


result["solution"]["nw"]["1"]
result["solution"]["nw"]["50160"]

#save("../../package_results/cordoba/result_ukdedk_stochFULL.jld2","result",result)
#save("../../package_results/cordoba/mn_data_ukdedk_stochFULL.jld2","mn_data",mn_data)

#resultACDC = _CBD.cordoba_mp_acdctnepopf(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

profit_1to2=0
profit_1to3=0
profit_2to3=0
profit_1to4=0
profit_2to4=0
profit_3to4=0
wnd_energy=0
wnd_profit=0
capex=0
mn_data["nw"]["1"]["gen"]
for (i,n) in result["solution"]["nw"]
    #if (parse(Int,i)>5*8360 && parse(Int,i)<6*8360)
        wnd_energy=wnd_energy+n["gen"]["4"]["pg"]
        wnd_profit=wnd_profit+n["gen"]["4"]["pg"]*max(mn_data["nw"][i]["gen"]["1"]["cost"][1],mn_data["nw"][i]["gen"]["2"]["cost"][1],mn_data["nw"][i]["gen"]["3"]["cost"][1])
        #profit_1to2=profit_1to2+n["branchdc_ne"]["3"]["pt"]*(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["1"]["cost"][1])
        profit_1to3=profit_1to3+n["branchdc_ne"]["5"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["1"]["cost"][1])
        profit_2to3=profit_2to3+n["branchdc_ne"]["9"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["2"]["cost"][1])
        profit_1to4=profit_1to4+n["branchdc_ne"]["13"]["pt"]*(mn_data["nw"][i]["gen"]["1"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        profit_2to4=profit_2to4+n["branchdc_ne"]["17"]["pt"]*(mn_data["nw"][i]["gen"]["2"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        #profit_3to4=profit_3to4+n["branchdc_ne"]["23"]["pt"]*(mn_data["nw"][i]["gen"]["3"]["cost"][1]-mn_data["nw"][i]["gen"]["4"]["cost"][1])
        capex=capex+(mn_data["nw"][i]["branchdc_ne"]["5"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["9"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["13"]["cost"]+mn_data["nw"][i]["branchdc_ne"]["17"]["cost"])#+mn_data["nw"][i]["branchdc_ne"]["23"]["cost"])
    #end
end
println(string(profit_1to2)*" "*
string(profit_1to3)*" "*
string(profit_2to3)*" "*
string(profit_1to4)*" "*
string(profit_2to4)*" "*
string(profit_3to4)*" "*
string(capex)*" "*
string(result["objective"])*" "*
string(wnd_energy)*" "*
string(wnd_profit)*" "*
string(8360)*" "*
string(scenario["planning_horizon"]))




uk2dkwind=deepcopy(resultACDC)
uk2dkwind_dta=deepcopy(mn_data)


for (key,nw) in resultACDC["solution"]["nw"]
    for (key2,b) in nw["ne_storage"]
    #println(key*" -  "*key2*" dr:  "*string(b["sd_ne"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["discharge_rating"]))
    #println(key*" -  "*key2*" cr: "*string(b["sc_ne"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["charge_rating"]))
    #println(key*" -  "*key2*" cr: "*string(b["sc_ne"]*mn_data["nw"][key]["ne_storage"][key2]["cost_abs"])*" "*" dr:  "*string(b["sd_ne"]*mn_data["nw"][key]["ne_storage"][key2]["cost_inj"]))
    println(key*" -  "*key2*" cst: "*string(mn_data["nw"][key]["ne_storage"][key2]["inst_cost"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["eq_cost"]))
    #println(key*" -  "*key2*" cst: "*string(mn_data["nw"][key]["ne_storage"][key2]["cost_abs"])*" "*string(mn_data["nw"][key]["ne_storage"][key2]["cost_inj"]))
end;end
