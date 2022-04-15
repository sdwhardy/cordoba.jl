################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections
##################### Input parameters #################################
rt_ex=pwd()*"\\test\\data\\input\\UK_DE_DK_wGenTypes\\"#folder path
argz = Dict(
"test"=>false,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
#NT2025 is common base case for all so should always be included
"scenario_names"=>["NT2025","DE2030","DE2040"],#["NT2025","DE2030","DE2040","GA2040","NT2040","GA2030","NT2030"],#Options: ["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
"k"=>2,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k Available: 2, 5, 10, 50, 100
#weather data years considered
"res_years"=>["2012","2013","2014"],#Options: ["2012","2013","2014","2015","2016"]
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
"owpp_mva"=>[4000],#mva of wf in MVA
"conv_lim"=>4000,#Max Converter size in MVA
"candidate_ics_ac"=>[6/5,1,4/5,3/5,2/5,1/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[6/5,1,4/5,3/5,2/5,1/5],#DC Candidate Cable sizes (fraction of full MVA)
"dr"=>0.04,#discount rate
"yearly_investment"=>1000000)

################## optimization/solver setup options ###################
s = Dict("output" => Dict("branch_flows" => false),
"home_market"=>[],#nodes within Zonal market
"balancing_reserve"=>0.3,#zonal market must be defined to have any effect
"AC"=>"1",#0=false, 1=true
"eps"=>0.0001,#admm residual (100kW)
"beta"=>5.5,
"relax_problem" => true,
"conv_losses_mp" => false,
"process_data_internally" => false,
"corridor_limit" => true,
"strg_lim_offshore"=>0.2,
"strg_lim_onshore"=>10,
"max_invest_per_year"=>_CBD.max_invest_per_year(argz))
########################################################################
################## Run MIP Formulation ###################
mn_data, data, argz, s = _CBD.main_ACDC_wgen_types(rt_ex,argz, s);#Build data structure for given options

#multi period problem setup
function multi_period_setup_wgen_type(ls,scenario_data,data, markets, infinite_grid, argz, s)
    #################### Multi-period input parameters #######################
    data,scenario, dim = _CBD.multi_period_stoch_year_setup_wgen_type(ls,argz["res_years"],argz["scenario_years"],argz["scenario_names"],data);
    scenario["planning_horizon"] = argz["scenario_planning_horizon"]; # in years, to scale generation cost
    extradata,data =create_profile_sets_mesh_wgen_type(dim, data, scenario_data, markets, infinite_grid, [data["baseMVA"] for wf in argz["owpp_mva"]])
    extradata = create_profile_sets_rest(dim, extradata, data)
    #########################################################################
    #################### Scale cost data
    scale_cost_data_hourly!(extradata, scenario)
    extradata=costs_datas_wREZ(extradata, s, argz)

    ######################################################
    xtradata = Dict{String,Any}()
    xtradata["dim"] = Dict{String,Any}()
    xtradata["dim"] = dim
    ######################################################
    # Create data dictionary where time series data is included at the right place
    mn_data = _PMACDC.multinetwork_data(data, xtradata, Set{String}(["source_type", "scenario", "scenario_prob", "name", "source_version", "per_unit"]))
    return mn_data, extradata
end


#seperates wfs from genz and defines markets/wfs zones
function main_ACDC_wgen_types(rt_ex,argz, s)
    ################# Load topology files ###################################
    _CBD.topology_df(rt_ex, s["relax_problem"], s["AC"])#creates .m file
    data, ics_ac, ics_dc, nodes = _CBD.filter_mfile_cables(rt_ex)#loads resulting topology and filters for candidate cables
    ############### defines size and market of genz and wfs ###################
	scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4UKBEDEDK.jld2")
    infinite_grid, genz, wfz, markets,all_gens = _CBD.gen_types(argz["owpp_mva"],nodes,data, scenario_data)
    push!(argz,"genz"=>genz)
    push!(argz,"wfz"=>wfz)
    #################### Calculates cable options for AC lines
    data=_CBD.AC_cable_options(data,argz["candidate_ics_ac"],ics_ac,data["baseMVA"])
    _CBD.print_topology_data_AC(data,markets)#print to verify
    #################### Calculates cable options for DC lines
    data=_CBD.DC_cable_options(data,argz["candidate_ics_dc"],ics_dc,data["baseMVA"])
    _CBD.additional_params_PMACDC(data)
    _CBD.print_topology_data_DC(data,markets)#print to verify
    ##################### load time series data ##############################
    scenario_data = _CBD.load_time_series_gentypes(rt_ex,argz, scenario_data,markets)
    ##################### multi period setup #################################
	s=_CBD.update_settings(s, argz, data)
    mn_data, xd  = multi_period_setup_wgen_type(ls, scenario_data, data, markets, infinite_grid, argz, s)
	s["xd"]=xd
    return  mn_data, data, argz, s
end

gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
#with Home market in Germany
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print solution
#with OBZ
s["home_market"]=[]#nodes within Zonal market
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print solution

################## Run Convex Formulation ################
s["relax_problem"]=true
mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print convex solution
mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)#Convert convex to candicate cables
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)#In large problems a larger MIPgap (1-5%?) may be desirable
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
_CBD.print_solution_data(result_mip, data, argz)#print MIP solution=#

################## Run ADMM Formulation ################
s["relax_problem"]=true
mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s);#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)
result_mip = _CBD.admm_4_AjAwAgAuAo_main(mn_data, gurobi, s);#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print Convex (ADMM) solution

mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)#Convert convex to candicate cables
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)#In large problems a larger MIPgap (1-5%?) may be desirable
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print MIP solution
