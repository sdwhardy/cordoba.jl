################## loads external packages ##############################
using Ipopt, Gurobi, JuMP, FileIO, JLD2, Dates, OrderedCollections
using PyCall; ks = pyimport_conda("kshape.core", "kshape.core")
import cordoba; const _CBD = cordoba#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
#import FlexPlan; const _FP = FlexPlan
#import InfrastructureModels; const _IM = InfrastructureModels
include("../aux/post_process/functions.jl")

################### ENTSO-E scenario description ####################################
scenario_names=["EU19"]
#scenario_names=["EU20"]
scenario_years=["2020","2030"]
#scenario_years=["2020"]
scenario_planning_horizon=30

################### Input data for script_test.jl ########################
owpp_mva=[4000]#mva of wf (in de)
#interconnectors format: (mva,km)
#must be in same order as .m file ie uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3), wf(4)-uk(1), wf(4)-de(2), wf(4)-dk(3)
ics=[(2000,550),(2000,760),(2000,250),(2000,470),(2000,145),(2000,246)];
ic_lim=maximum([_CBD.DC_cbl(ic,1).num*_CBD.DC_cbl(ic,1).elec.mva for ic in first.(ics)])

#location of nodes
markets_wfs=[["UK","DE","DK"],["DE"]]#must be in same order as .m file gens
infinite_grid=sum(first.(ics))+sum(owpp_mva)#ensures enough generation and consumption in all markets

##################### load time series data ##############################
k=2
scenario_data = _CBD.get_scenario_year_tss(scenario_names,scenario_years)#Retrieve the scenario time series

##################### Cluster time series data ###########################
for (sc,yrs_ts) in scenario_data
    for (yr,ts) in yrs_ts
        daily_ts=_CBD.daily_tss(DateTime.(ts[!,"time_stamp"]))
        daily_dk=_CBD.daily_tss(Float64.(ts[!,"Wnd_MWhDK"]))
        daily_de=_CBD.daily_tss(Float64.(ts[!,"EUR_daDE"]))
        daily_uk=_CBD.daily_tss(Float64.(ts[!,"EUR_daUK"]))

        kshape_clusters_deuk=ks.kshape(ks.zscore(sqrt.((daily_uk.^2).+(daily_de.^2).+(daily_dk.^2))',axis=1), k)
        ts=Vector{Float64}(); for clusters in last.(kshape_clusters_deuk); if (length(clusters)>0); ts=vcat(ts,daily_ts[:,rand(clusters)+1]); end;end
        filter!(row -> row.time_stamp in ts, scenario_data[sc][yr])
    end
end

#scenario_data=load("./test/data/input/EUSTDG_TS_k"*string(k)*".jld2")

##################### Find minimum length scenario
#scenario_data=load("./test/data/input/EUSTDG_TS_k5.jld2")
ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
scenario_data[_sc][_yr]=scenario_data[_sc][_yr][1:2,:]
push!(ls,length(scenario_data[_sc][_yr].time_stamp))
end;end;ls=minimum(ls)

##################### Make all scenarios the same length
for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
end;end
#save("./test/data/input/EUSTDG_TS_k"*string(k)*".jld2",scenario_data)

#optimation settings
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false, "scenarios_length" => length(scenario_names), "years_length" => length(scenario_years), "hours_length" => ls, "corridor_limit" => true, "hours_length" => ls, "ic_lim"=>ic_lim, "compare_mode" => true)
#select solver
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)


################################################# MIP TNEP #######################################
################## reads .m input file name ######################
casename = "test_convex_conv"
file = "./test/data/input/$casename.m"
data_mip = _PM.parse_file(file)#load data in PM format

#################### Calculates cables for DC lines
z_base_dc=(data_mip["busdc_ne"]["1"]["basekVdc"])^2/data_mip["baseMVA"]
#candidate_ics=[1,4/5,1/2,1/4]#Candidate Cable sizes
#candidate_ics=[1,4/5,3/5,1/2]#Candidate Cable sizes
candidate_ics=[1]#Candidate Cable sizes
data_mip=_CBD.additional_candidatesICS(data_mip,candidate_ics,ics)#adds additional candidates
for (i,bdc) in data_mip["branchdc_ne"]
data_mip["branchdc_ne"][i]=_CBD.candidateIC_cost_impedance(bdc,z_base_dc);end
data_mip["branchdc_ne"]=_CBD.unique_candidateIC(data_mip["branchdc_ne"])#keep only unique candidates
_PMACDC.process_additional_data!(data_mip)#add extra DC model data

#################### Multi-period input parameters #######################
all_scenario_data,data_mip,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data_mip)
scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
extradata,data_mip = _CBD.create_profile_sets_mesh(dim, data_mip, all_scenario_data, markets_wfs, infinite_grid, owpp_mva)

#################### Scale cost data
_CBD.scale_cost_data_per_scenario!(data_mip, scenario)
_CBD.scale_cost_data_per_scenario!(extradata, scenario)

# Create data dictionary where time series data is included at the right place
mn_data_mip = _PMACDC.multinetwork_data(data_mip, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))

#run optimization
result_mip = _CBD.cordoba_acdc_tnep_convex_conv(mn_data_mip, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

#=
################################################# CONVEX TNEP #######################################
################## reads .m input file name ######################
casename = "test_convexafy"
file = "./test/data/input/$casename.m"
data_conv = _PM.parse_file(file)#load data in PM format
_PMACDC.process_additional_data!(data_conv)#add extra DC model data

#################### Multi-period input parameters #######################
all_scenario_data,data_conv,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data_conv)
scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
extradata,data_conv = _CBD.create_profile_sets_mesh(dim, data_conv, all_scenario_data, markets_wfs, infinite_grid, owpp_mva)

# Scale cost data
_CBD.scale_cost_data_per_scenario!(data_conv, scenario)
_CBD.scale_cost_data_per_scenario!(extradata, scenario)

# Create data dictionary where time series data is included at the right place
mn_data_conv = _PMACDC.multinetwork_data(data_conv, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "per_unit", "source_version"]))

#run optimization
result_conv = _CBD.cordoba_mp_acdctnepopf_convexafy(mn_data_conv, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s);

################################################### Re build Convex solution with real cables #########################################
cost=0;off_shore_nodes=[4];
_hrs=scenario["hours"]*scenario["years"]
konverters=result_conv["solution"]["nw"]["1"]["convdc"]
kable_indices=[(i,(k["p_rateA"],last(ics[parse(Int64,i)]))) for (i,k) in result_conv["solution"]["nw"]["1"]["branchdc"]]

################## reads .m input file name ######################
casename = "test_opf_dc"
file = "./test/data/input/$casename.m"
data_opf = _PM.parse_file(file)#load data in PM format

################## Actual cables and converters #########################
data_opf,cost=rebuild_system_cables_convex(kable_indices,data_opf,cost,z_base_dc)#cables
data_opf,cost=rebuild_system_converters(konverters,data_opf, off_shore_nodes,cost)#4619
_PMACDC.process_additional_data!(data_opf)#add extra DC model data

#################### Multi-period input parameters #######################
all_scenario_data,data_opf,scenario, dim = _CBD.multi_period_stoch_year_setup(ls,scenario_years,scenario_names,scenario_data,data_opf)
scenario["planning_horizon"] = scenario_planning_horizon # in years, to scale generation cost
extradata,data_opf = _CBD.create_profile_sets_mesh(dim, data_opf, all_scenario_data, markets_wfs, infinite_grid, owpp_mva)


# Scale cost data
_CBD.scale_cost_data_per_scenario!(extradata, scenario)

# Create data dictionary where time series data is included at the right place
mn_data_opf = _PMACDC.multinetwork_data(data_opf, extradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "per_unit", "source_version"]))

#run optimization
result_opf = _CBD.cordoba_mp_acdctnepopf_stoch(mn_data_opf, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)

############################################## Printing results ##################################################
println()
println("MIP: "*string(result_mip["objective"])*" Convex: "*string(result_opf["objective"]+cost))
println()
for (k,b) in result_conv["solution"]["nw"]["1"]["branchdc"]
    if (b["p_rateA"]>0)
        println(string(mn_data_opf["nw"]["1"]["branchdc"][k]["fbusdc"])*" - "*string(mn_data_opf["nw"]["1"]["branchdc"][k]["tbusdc"])*" CONVEX_MVA: "*string(b["p_rateA"])*" CONV2MIP_MVA: "*string(mn_data_opf["nw"]["1"]["branchdc"][k]["rateA"]*100)*" DIF [%]: "*string(((b["p_rateA"]-mn_data_opf["nw"]["1"]["branchdc"][k]["rateA"]*100)/b["p_rateA"])*100))
    end
end=#
println()
for (k,b) in result_mip["solution"]["nw"]["1"]["branchdc_ne"]
    if (b["isbuilt"]==1)
        println(k*" "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["fbusdc"])*" - "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["tbusdc"])*" MVA: "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["rateA"]*100))
    end
end
println()
[println(string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["fbusdc"])*" - "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["tbusdc"])*" MVA: "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["rateA"]*100)*" Cost: "*string(mn_data_mip["nw"]["1"]["branchdc_ne"][k]["cost"]*scenario["hours"])) for (k,b) in result_mip["solution"]["nw"]["1"]["branchdc_ne"]]
