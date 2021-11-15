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
scenario_years=["2020","2030","20"]
scenario_planning_horizon=30

################### Input data for script_test.jl ########################
owpp_mva=[4000]#mva of wf (in de)
#interconnectors format: (mva,km)
#must be in same order as .m file ie uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3), wf(4)-uk(1), wf(4)-de(2), wf(4)-dk(3)
ics=[(2000,550),(2000,760),(2000,250),(2000,470),(2000,145),(2000,246)];
#location of nodes
markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

################## reads .m input file name ######################
casename = "test_opf"
file = "./test/data/input/$casename.m"
data = _PM.parse_file(file)#load data in PM format40

################## fill in solution topo #########################
#cables
z_base_dc=(data["busdc"]["1"]["basekVdc"])^2/data["baseMVA"]
kables=result_conv["solution"]["nw"]["2160"]["branchdc"]
cost=0;key="1";brn=kables["1"]
for (key,brn) in kables
    if (brn["p_rateA"]>0)
        cb=_CBD.DC_cbl(brn["p_rateA"]*10^2,last(ics[parse(Int64,key)]))
        data["branchdc"][key]["rateA"]=data["branchdc"][key]["rateB"]=data["branchdc"][key]["rateC"]=cb.num*cb.elec.mva;
        data["branchdc"][key]["cost"]=cb.costs.cpx_i+cb.costs.cpx_p
        data["branchdc"][key]["r"]=((cb.elec.ohm*10^3/cb.num)*cb.length)/z_base_dc
        cost=cost+data["branchdc"][key]["cost"]
    else
        delete!(data["branchdc"], key)
    end;end

#=kables=Dict{String, Any}()
for i in [("4","14"),("5","18"),("6","24")]; push!(kables,first(i) => mn_data["nw"]["1"]["branchdc_ne"][last(i)]); end
cost=0
for (key,brn) in kables
    if (brn["rateA"]>0)
        km=last(ics[parse(Int64,key)])
        cb=_CBD.DC_cbl(brn["rateA"]*100, km)
        data["branchdc"][key]["rateA"]=data["branchdc"][key]["rateB"]=data["branchdc"][key]["rateC"]=cb.num*cb.elec.mva;
        data["branchdc"][key]["cost"]=cb.costs.cpx_p+cb.costs.cpx_i
        cost=cost+data["branchdc"][key]["cost"]
    else
        data["branchdc"][key]["rateA"]=data["branchdc"][key]["rateB"]=data["branchdc"][key]["rateC"]=0
        data["branchdc"][key]["cost"]=0
    end;end=#

#converters
konverters=result_conv["solution"]["nw"]["2160"]["convdc"]
offshore=[4]
for (key,cnv) in konverters
    cb=_CBD.DC_cbl(cnv["p_pacmax"]*100, 1)
    data["convdc"][key]["Pacmax"]=4000#cb.num*cb.elec.mva
    data["convdc"][key]["Pacmin"]=-4000#-1*cb.num*cb.elec.mva;
    if (issubset(parse(Int64,key),offshore))
        data["convdc"][key]["cost"]=0.1925*(cb.num*cb.elec.mva)+(0.0975*cb.num*cb.elec.mva)
        0.1925*(4000)+(0.0975*4000)
        cost=cost+data["convdc"][key]["cost"]
    else
        data["convdc"][key]["cost"]=0.1925*(cb.num*cb.elec.mva)
        cost=cost+data["convdc"][key]["cost"]
    end
end

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
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls)

#run optimizationusing Ipopt

result_opf = _CBD.cordoba_mp_acdc_opf(mn_data, _PM.ACPPowerModel, Ipopt.Optimizer, multinetwork=true; setting = s)
#MIP: -32998.6591159048 cost: 5816.865527954101
#convex: -32998.659115904746 cost: 3920.5945318603517
