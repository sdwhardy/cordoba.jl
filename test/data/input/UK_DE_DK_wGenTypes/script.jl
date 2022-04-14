################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels

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
#seperates wfs from genz and defines markets/wfs zones

#cost[1], pmax, source_id[2]
#load Time series data
function load_time_series_gentypes(rt_ex, argz, scenario_data,markets)
    #keep only specified scenarios
    d_keys=keys(scenario_data["Generation"]["Scenarios"]);for k in d_keys;if !(issubset([string(k)],argz["scenario_names"]));delete!(scenario_data["Generation"]["Scenarios"],k);end;end
	d_keys=keys(scenario_data["Demand"]);for k in d_keys;if !(issubset([string(k)],argz["scenario_names"]));delete!(scenario_data["Demand"],k);end;end
	#Keep only specified markets
	countries=unique(vcat(markets[1],markets[2]))
	for key in keys(scenario_data["Generation"]["RES"]["Offshore Wind"])
		if !(issubset([string(key)],countries));
			delete!(scenario_data["Generation"]["RES"]["Offshore Wind"],key);
			delete!(scenario_data["Generation"]["RES"]["Onshore Wind"],key);
			delete!(scenario_data["Generation"]["RES"]["Solar PV"],key);
		end;end
	#keep only specified weather years
	for country in keys(scenario_data["Generation"]["RES"]["Offshore Wind"])
		for year in keys(scenario_data["Generation"]["RES"]["Offshore Wind"][country])
			if !(issubset([string(year)],argz["res_years"]));
				delete!(scenario_data["Generation"]["RES"]["Offshore Wind"][country],year);
				delete!(scenario_data["Generation"]["RES"]["Onshore Wind"][country],year);
				delete!(scenario_data["Generation"]["RES"]["Solar PV"][country],year);
			end;end;end
	if (haskey(argz, "test") && argz["test"]==true)
        for k0 in d_keys; for k1 in keys(scenario_data[k0]); scenario_data[k0][k1]=scenario_data[k0][k1][1:2,:];end;end
    end
    ##################### Find minimum length scenario and Make all scenarios the same length
    ls=[];for (_sc, data_by_scenario) in scenario_data; for (_yr, data_by_yr) in data_by_scenario;
    push!(ls,length(scenario_data[_sc][_yr].time_stamp))
    end;end;ls=minimum(ls)

    for (_yr, data_by_yr) in scenario_data; for (_sc, data_by_scenario) in data_by_yr;
    scenario_data[_yr][_sc]=scenario_data[_yr][_sc][1:ls,:]
    end;end
    return scenario_data, ls
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
    scenario_data, ls = load_time_series_gentypes(rt_ex,argz, scenario_data,markets)
    push!(argz,"ls"=>ls)
    ##################### multi period setup #################################
	s=update_settings(s, argz, data)
    mn_data, xd  = multi_period_setup(ls, scenario_data, data, markets, infinite_grid, argz, s)
	s["xd"]=xd
    return  mn_data, data, argz, s=#
    return argz
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
