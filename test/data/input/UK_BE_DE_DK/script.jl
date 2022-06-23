################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, OrderedCollections
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
##################### File parameters #################################
s = Dict(
"rt_ex"=>pwd()*"\\test\\data\\input\\UK_BE_DE_DK\\",#folder path
"scenario_data_file"=>"C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2",
################# temperal parameters #################
"test"=>false,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
"scenario_names"=>["NT","DE","GA"],#["NT","DE","GA"]
"k"=>6,#number of representative days modelled (24 hours per day)//#best for maintaining mean/max is k=6 2014, 2015
"res_years"=>["2014","2015"],#Options: ["2012","2013","2014","2015","2016"]//#best for maintaining mean/max is k=6 2014, 2015
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
"dr"=>0.04,#discount rate
"yearly_investment"=>1000000,
################ electrical parameters ################
"AC"=>"1",#0=false, 1=true
"owpp_mva"=>[2000,4000,4000],#mva of wf in MVA
"conv_lim_onshore"=>3000,#Max Converter size in MVA
"conv_lim_offshore"=>4000,#Max Converter size in MVA
"strg_lim_offshore"=>0.2,
"strg_lim_onshore"=>10,
"candidate_ics_ac"=>[1,4/5,3/5,2/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1,4/5,3/5,2/5],#DC Candidate Cable sizes (fraction of full MVA)
################## optimization/solver setup options ###################
"output" => Dict("branch_flows" => false),
"eps"=>0.0001,#admm residual (100kW)
"beta"=>5.5,
"relax_problem" => false,
"conv_losses_mp" => true,
"process_data_internally" => false,
"corridor_limit" => true)

################## Run MIP Formulation ###################
#NOTE only very basic intuitive check passed on functions wgen_type
s["home_market"]=[]
@time result_mip, data, mn_data, s= _CBD.nodal_market_main(s)
_CBD.print_solution_wcost_data(result_mip, s, data)#-856896.0245340846
results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s)#primal: -565150.39, dual: -565819.65
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\nodal_market_k6_VOLL5000b.jld2",results)
social_welfare = _CBD.SocialWelfare(s, result_mip, mn_data, data)

s["home_market"]=[[5,6,7]]
@time result_mip, data, mn_data, s, result_mip_hm_prices = _CBD.zonal_market_main(s);
s["cost_summary"]=_CBD.print_solution_wcost_data(result_mip, s, data)#-856559.087752747 (MIP)
results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s, "result_mip_hm_prices"=>result_mip_hm_prices)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_567_k6_VOLL5000b_rc.jld2",results)

s, result_mip, data, mn_data=_CBD.summarize_zonal_in_s(results);
_CBD.print_table_summary(s)

#Check balanced
all_gens=[];all_strg=[]
for (n,nw) in result_mip["solution"]["nw"]
push!(all_gens,sum(gen["pg"] for (g,gen) in nw["gen"]));
push!(all_strg,sum(stg["ps"] for (s,stg) in nw["storage"]));end
minimum(all_strg)





################################################
scenario_data=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2")

FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2",scenario_data)


sdgs=deepcopy(scenario_data["Generation"]["Scenarios"])
scenario_data["Generation"]["Scenarios"]=Dict()
push!(scenario_data["Generation"]["Scenarios"],"Base"=>Dict())
push!(scenario_data["Demand"],"NT"=>Dict())
push!(scenario_data["Demand"],"GA"=>Dict())
push!(scenario_data["Demand"],"DE"=>Dict())
push!(scenario_data["Generation"]["Scenarios"]["Base"],"2020"=>sdgs["Base"])
push!(scenario_data["Demand"]["NT"],"2030"=>sdgs["NT2030"])
push!(scenario_data["Demand"]["GA"],"2030"=>sdgs["GA2030"])
push!(scenario_data["Demand"]["DE"],"2030"=>sdgs["DE2030"])
push!(scenario_data["Demand"]["NT"],"2040"=>sdgs["NT2040"])
push!(scenario_data["Demand"]["GA"],"2040"=>sdgs["GA2040"])
push!(scenario_data["Demand"]["DE"],"2040"=>sdgs["DE2040"])






#=function global_optimal_dispatch(s)
    #s["relax_problem"]=true
    data, s = _CBD.get_topology_data(s)#topology.m file
    scenario_data = _CBD.get_scenario_data(s)#scenario time series
	###########################################################################
	all_gens,s = _CBD.gen_types(data,scenario_data,s)
    #################### Calculates cable options for AC lines
    data = _CBD.AC_cable_options(data,s["candidate_ics_ac"],s["ics_ac"],data["baseMVA"])
    #################### Calculates cable options for DC lines
    data = _CBD.DC_cable_options(data,s["candidate_ics_dc"],s["ics_dc"],data["baseMVA"])
    if (haskey(s, "home_market") && length(s["home_market"])>0);data = _CBD.keep_only_hm_cables(s,data);end#if home market reduce to only those in
    _CBD.additional_params_PMACDC(data)
    _CBD.print_topology_data_AC(data,s["map_gen_types"]["markets"])#print to verify
    _CBD.print_topology_data_DC(data,s["map_gen_types"]["markets"])#print to verify
    ##################### load time series data ##############################
    scenario_data = _CBD.load_time_series_gentypes(s, scenario_data)
    ##################### multi period setup #################################
	s = _CBD.update_settings_wgenz(s, data)
    mn_data, s  = _CBD.multi_period_setup_wgen_type(scenario_data, data, all_gens, s);
	push!(s,"max_invest_per_year"=>_CBD.max_invest_per_year(s))
    return  mn_data, data, s
end

function zonal_market_main(s)
    hm=deepcopy(s["home_market"]);
    #s["relax_problem"]=true
    #mn_data, data, s = global_optimal_dispatch(s);#Build data structure for given options
    mn_data, data, s = _CBD.data_setup_zonal(s);#Build data structure for given options
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1);#select solver
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s);#Solve problem
    _CBD.print_solution_wcost_data(result_mip, s, data);
    mn_data, data, s = data_setup_nodal(s);#Build data structure for given options
    mn_data, s = set_inter_zonal_grid(result_mip,mn_data,s);
    s["home_market"]=[]
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)#select solver
    result_mip = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    #print_solution_wcost_data(result_mip, s, data)
    s["home_market"]=hm
    s["rebalancing"]=true
    s["relax_problem"]=true
    s["output"]["duals"]=true
    mn_data, data, s = data_update(s,result_mip);#Build data structure for given options
    mn_data, s = set_rebalancing_grid(result_mip,mn_data,s);
    s, mn_data= remove_integers(result_mip,mn_data,data,s);
    result_mip_hm_prices = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    s["home_market"]=[]
    mn_data, data, s = data_update(s,result_mip);#Build data structure for given options
    mn_data, s = set_rebalancing_grid(result_mip,mn_data,s);
    s, mn_data= remove_integers(result_mip,mn_data,data,s);
    result_mip = cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    result_mip= hm_market_prices(result_mip, result_mip_hm_prices)
    return result_mip, data, mn_data, s
end=#
