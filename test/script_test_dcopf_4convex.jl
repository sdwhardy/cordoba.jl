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
#scenario_names=["EU19"]
#scenario_years=["2020","2030","2040"]
#scenario_years=["2020"]
scenario_planning_horizon=30

################### Input data for script_test.jl ########################
owpp_mva=[4000]#mva of wf (in de)
#interconnectors format: (mva,km)
#must be in same order as .m file ie uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3), wf(4)-uk(1), wf(4)-de(2), wf(4)-dk(3)
ics=[(4000,550),(4000,760),(4000,250),(4000,470),(4000,145),(4000,246)];
#location of nodes
markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

cost=0
_hrs=scenario["hours"]*scenario["years"]
z_base_dc=(data["busdc"]["1"]["basekVdc"])^2/data["baseMVA"]
konverters=result_conv["solution"]["nw"]["1"]["convdc"]
off_shore_nodes=[4];
kable_indices=[(i,(k["p_rateA"],last(ics[parse(Int64,i)]))) for (i,k) in result_conv["solution"]["nw"]["1"]["branchdc"]]
#for i in [("1","1"),("2","1"),("3","1"),("4","13"),("5","17"),("6","21")]; push!(kables,first(i) => data["branchdc_ne"][last(i)]); end


################## reads .m input file name ######################
casename = "test_opf_dc"
file = "./test/data/input/$casename.m"
data = _PM.parse_file(file)#load data in PM format

################## fill in solution topo #########################
#cables
data,cost=rebuild_system_cables_convex(kable_indices,data,cost,z_base_dc)

data,cost=rebuild_system_converters(konverters,data, off_shore_nodes,cost)#4619



#################### PowerModelsACDC settings ############################
_PMACDC.process_additional_data!(data)#add extra DC model data

##################### load time series data ##############################
#scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series
#=scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series

##################### Cluster time series data ###########################
k=100
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
ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
#scenario_data[_sc][_yr]=scenario_data[_sc][_yr][1:2,:]
push!(ls,length(scenario_data[_sc][_yr].time_stamp))
end;end;ls=minimum(ls)

##################### Make all scenarios the same length
for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    scenario_data[_sc][_yr]=scenario_data[_sc][_yr][1:ls,:]
end;end
#delete!(scenario_data["EU19"], "2040")

#################### Multi-period input parameters #######################
all_scenario_data,data,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data)
scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
extradata,data = _CBD.create_profile_sets_mesh(dim, data, all_scenario_data, markets_wfs, infinite_grid, owpp_mva)

# Scale cost data
_CBD.scale_cost_data_cordoba!(extradata, scenario)

# Create data dictionary where time series data is included at the right place
mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "per_unit", "source_version"]))

#select solver
#gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

#optimation settings
scenes_years=[];for (i,(k,_scene)) in enumerate(scenario["sc_names"]);push!(scenes_years,[]);for (j, _yr) in _scene; push!(scenes_years[i],_yr);end;end
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls, "corridor_limit" => false, "compare_mode" => true)

#run optimizationusing Ipopt

result_opf = _CBD.cordoba_mp_acdctnepopf_stoch(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
println(result_opf["objective"]+cost)
println(result_opf["objective"]+cost)




function rebuild_system_converters(konverters,data_struct, off_shore_nodes,cost)
    for (key,cnv) in konverters
        data_struct["convdc"][key]["Pacmin"]=-1*cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Pacmax"]=cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Qacmin"]=-1*cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Qacmax"]=cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Pacrated"]=cnv["p_pacmax"]*100
        data_struct["convdc"][key]["Qacrated"]=cnv["p_pacmax"]*100
        if (issubset(parse(Int64,key),off_shore_nodes))
            data_struct["convdc"][key]["cost"]=0.29*cnv["p_pacmax"]*100
        else
            data_struct["convdc"][key]["cost"]=0.1925*cnv["p_pacmax"]*100
        end
        cost=cost+data_struct["convdc"][key]["cost"]
    end
    return data_struct,cost
end

function rebuild_system_cables_convex(kable_indices,data,cost,z_base)
    for k in kable_indices
        if (first(last(k))>0)
            cb=_CBD.DC_cbl(first(last(k)),last(last(k)))
            data["branchdc"][first(k)]["rateA"]=data["branchdc"][first(k)]["rateB"]=data["branchdc"][first(k)]["rateC"]=cb.num*cb.elec.mva
            data["branchdc"][first(k)]["cost"]=cb.costs.cpx_i+cb.costs.cpx_p
            data["branchdc"][first(k)]["r"]=((cb.elec.ohm*10^3/cb.num)*cb.length)/z_base
            cost=cost+data["branchdc"][first(k)]["cost"]
        else
            data["branchdc"][first(k)]["rateA"]=data["branchdc"][first(k)]["rateB"]=data["branchdc"][first(k)]["rateC"]=0
            data["branchdc"][first(k)]["cost"]=0
        end
    end
    return data,cost
end
