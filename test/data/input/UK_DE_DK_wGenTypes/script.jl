
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
"scenario_names"=>["NT","DE"],#["NT","DE","GA"]
"k"=>2,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k Available: 2, 5, 10, 50, 100
#weather data years considered
"res_years"=>["2012","2016"],#Options: ["2012","2013","2014","2015","2016"]
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
"owpp_mva"=>[4000],#mva of wf in MVA
"conv_lim"=>4000,#Max Converter size in MVA
"candidate_ics_ac"=>[6/5,1,4/5,3/5,2/5,1/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[6/5,1,4/5,3/5,2/5,1/5],#DC Candidate Cable sizes (fraction of full MVA)
"dr"=>0.04,#discount rate
"yearly_investment"=>10000)

################## optimization/solver setup options ###################
s = Dict("output" => Dict("branch_flows" => false),
"balancing_reserve"=>0.3,#zonal market must be defined to have any effect
"rebalancing"=>false,#zonal market must be defined to have any effect
"AC"=>"1",#0=false, 1=true
"eps"=>0.0001,#admm residual (100kW)
"beta"=>5.5,
"relax_problem" => false,
"conv_losses_mp" => true,
"process_data_internally" => false,
"corridor_limit" => true,
"strg_lim_offshore"=>0.2,
"strg_lim_onshore"=>10)
#21.7175
########################################################################
################## Run MIP Formulation ###################
#NOTE only very basic intuitive check passed on functions wgen_type
mn_data_wf, data_wf, argz_wf, s_wf, map_gen_types_wf= _CBD.main_ACDC_wgen_types(rt_ex,deepcopy(argz), deepcopy(s));#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
result_mip_wf = _CBD.cordoba_acdc_wf_strg(mn_data_wf, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s_wf)#Solve problem
_CBD.print_solution_wcost_data(result_mip_wf, argz_wf, data_wf)
result_mip_wf=deepcopy(test)

s["rebalancing"]=true
s["relax_problem"]=true
s["output"]["duals"]=true
mn_data, data, argz, s, map_gen_types= _CBD.main_ACDC_wgen_types(rt_ex,argz, s);#Build data structure for given optionsS
mn_data, s = _CBD.set_rebalancing_grid(result_mip_wf,mn_data,s);
s, mn_data=_CBD.remove_integers(result_mip_wf,mn_data,data,s);
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_wcost_data(result_mip, argz, data)
###################################### plotting ################################
gen_consume_summary=_CBD.summarize_generator_solution_data(result_mip, data, argz, map_gen_types,s)#print solution


country="DK"#,"DE","DK"]
scenario="1"
con=gen_consume_summary["onshore_demand"][scenario][country]#[121:144,:]
gen=gen_consume_summary["onshore_generation"][scenario][country]#[121:144,:]
_CBD.plot_generation_profile(deepcopy(gen),deepcopy(con),country*" "*scenario)

country="DE"#,"DE","DK"]
scenario="1"
con=gen_consume_summary["onshore_demand"][scenario][country]#[121:144,:]
gen=gen_consume_summary["offshore_generation"][scenario][country]#[121:144,:]
con=select!(con,:ts)
_CBD.plot_generation_profile(deepcopy(gen),deepcopy(con),country*" "*scenario)

hourly_income_wf=_CBD.owpp_profit_obz(s, result_mip, keys(mn_data["scenario"][scenario]), "4", "163");
_CBD.plot_cumulative_income(hourly_income_wf, s["hours_length"])

hourly_income_tl=_CBD.transmission_line_profits(s, result_mip, keys(mn_data["scenario"][scenario]), data);
_CBD.plot_cumulative_income_tl(hourly_income_tl, s["hours_length"])

result_mip=_CBD.undo_marginal_price_scaling(s,argz,result_mip)
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (1,"UK"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (3,"DK"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (2,"DE"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (4,"WF"))



tss=keys(mn_data["scenario"][scenario])




##################################### HM market 
################## Run MIP Formulation ###################
#NOTE only very basic intuitive check passed on functions wgen_type
s["home_market"]=[]
s["wf_circuit"]=[2,4]
mn_data_wf, data_wf, argz_wf, s_wf, map_gen_types_wf= _CBD.main_ACDC_wgen_types(rt_ex,deepcopy(argz), deepcopy(s));#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
result_mip_wf = _CBD.cordoba_acdc_wf_strg(mn_data_wf, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s_wf)#Solve problem
_CBD.print_solution_wcost_data(result_mip_wf, argz_wf, data_wf)

s["home_market"]=[2,4]
s["wf_circuit"]=[]
mn_data, data, argz, s, map_gen_types = _CBD.main_ACDC_wgen_types(rt_ex,argz, s);#Build data structure for given options
mn_data, s = _CBD.set_zonal_grid(result_mip_wf,mn_data,s);
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
result_mip_hm = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_wcost_data(result_mip_hm, argz, data)

s["home_market"]=[2,4]
s["rebalancing"]=true
s["relax_problem"]=true
s["output"]["duals"]=true
mn_data, data, argz, s, map_gen_types= _CBD.main_ACDC_wgen_types(rt_ex,argz, s);#Build data structure for given options
result_mip= _CBD.combine_solutions(result_mip_hm,result_mip_wf)
mn_data, s = _CBD.set_rebalancing_grid(result_mip,mn_data,s);
s, mn_data=_CBD.remove_integers(result_mip,mn_data,data,s);
result_mip_hm_prices = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_wcost_data(result_mip, argz, data)

s["home_market"]=[]
mn_data, data, argz, s, map_gen_types= _CBD.main_ACDC_wgen_types(rt_ex,argz, s);#Build data structure for given optionsS
result_mip= _CBD.combine_solutions(result_mip_hm,result_mip_wf)
mn_data, s = _CBD.set_rebalancing_grid(result_mip,mn_data,s);
s, mn_data=_CBD.remove_integers(result_mip,mn_data,data,s);
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_wcost_data(result_mip, argz, data)
###################################### plotting ################################
result_mip=_CBD.hm_market_prices(result_mip, result_mip_hm_prices)
result_mip=deepcopy(result_mip_001)
gen_consume_summary=_CBD.summarize_generator_solution_data(result_mip, data, argz, map_gen_types,s)#print solution
#####################################








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
