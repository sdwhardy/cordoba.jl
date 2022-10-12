################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections
##################### File parameters #################################
s = Dict(
"rt_ex"=>pwd()*"\\test\\data\\input\\UK_DE_DK_wOnshore\\",#folder path
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
"owpp_mva"=>[4000],#mva of wf in MVA
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
"corridor_limit" => true,
"onshore_grid"=>true)
########################################################################
#0.0066 - branch
_CBD.problemINPUT_map(data, s)
##################################### HM market 
################## Run MIP Formulation ###################
#NOTE only very basic intuitive check passed on functions wgen_type
s["home_market"]=[]
mn_data, data, s = _CBD.data_setup(s);
@time result = _CBD.nodal_market_main(mn_data, data, s)#-3359431 -33899162 0.89%
_CBD.print_solution_wcost_data(result_mip, s, data)#-856896.0245340846 
results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_DE_DK\\nodal_results_VOLL5000b_onShore30convs.jld2",results)

s["home_market"]=[[3,4]]
@time result_mip, data, mn_data, s, result_mip_hm_prices = _CBD.zonal_market_main(s);
_CBD.print_solution_wcost_data(result_mip, s, data)#-856559.087752747 (MIP)
results=Dict("result_mip"=>result_mip,"data"=>data, "mn_data"=>mn_data, "s"=>s, "result_mip_hm_prices"=>result_mip_hm_prices)
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_DE_DK\\zonal_results_hm34_VOLL5000_onShore30convs.jld2",results)

#result_mip=deepcopy(result_mip_001)
@time gen_consume_summary=_CBD.summarize_generator_solution_data(result_mip, data,s)#print solution
#####################################

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

#hourly_income_wf=_CBD.owpp_profit_obz(s, result_mip, keys(mn_data["scenario"][scenario]), "4", "163");
_CBD.plot_cumulative_income(hourly_income_wf, s["hours_length"])
sum(hourly_income_wf["power"])#3262.67155
sum(hourly_income_wf["income"])#3262.67155

hourly_income_tl=_CBD.transmission_line_profits(s, result_mip, keys(mn_data["scenario"][scenario]), data);
_CBD.plot_cumulative_income_tl(hourly_income_tl, s["hours_length"])
sum(hourly_income_tl["dc"])#3262.67155

result_mip=_CBD.undo_marginal_price_scaling(s,argz,result_mip)
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (1,"UK"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (3,"DK"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (2,"DE"))
_CBD.plot_dual_marginal_price(result_mip, keys(mn_data["scenario"][scenario]), (4,"WF"))

social_welfare = SocialWelfare(s, result_mip, mn_data, data)
social_welfare["totals"]





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
