################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, OrderedCollections, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels


##################### Input parameters #################################
s = Dict(
"rt_ex"=>pwd()*"\\test\\data\\input\\test\\",#folder path if directly
"scenario_data_file"=>"C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4EU.jld2",
################# temperal parameters #################
"test"=>true,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
"scenario_names"=>["NT2025","NT2030","NT2040","DE2030","DE2040","GA2030","GA2040"],#["NT","DE","GA"]#,"NT2030","NT2040","DE2030","DE2040","GA2030","GA2040"
"k"=>4,#number of representative days modelled (24 hours per day)//#best for maintaining mean/max is k=6 2014, 2015
"res_years"=>["2014","2015"],#Options: ["2012","2013","2014","2015","2016"]//#best for maintaining mean/max is k=6 2014, 2015
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
################# Financial parameters ################
"dr"=>0.04,#discount rate
"yearly_investment"=>1000000,#max investment per modelling year
################ electrical parameters ################
"conv_lim_onshore"=>3000,#Max Converter size in MVA
"conv_lim_offshore"=>4000,#Max Converter size in MVA
"strg_lim_offshore"=>0.2,#Max offshore storage capacity
"strg_lim_onshore"=>10,#Max onshore storage capacity
"candidate_ics_ac"=>[1/5,2/5,3/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1/5,2/5,3/5],#DC Candidate Cable sizes (fraction of full MVA)[1,4/5,3/5,2/5]
################## optimization/solver setup options ###################
"relax_problem" => false,#binaries->continuous variables
"corridor_limit" => false,#limit cables in parallel?
"TimeLimit" => 259200,#solver max time in seconds
"MIPGap"=>1e-4,#max gap between MIP and convex solution 
"PoolSearchMode" => 0,#0-single solution, 1- poolsolutions of random quality, 2- poolsolutions of highest quality 
"PoolSolutions" => 1)#number of solutions to find
s=_CBD.hidden_settings(s)
s_z=deepcopy(s)

################## Run nodal Formulation ###################
#nodal data setup
s["home_market"]=[]
mn_data, data, s = _CBD.data_setup(s);
#display inputs
_CBD.problemINPUT_map(data, s)
_CBD.problemINPUT_mapNTCs(data, s)
#solve nodal market
result = _CBD.nodal_market_main(mn_data, data, s)
#display results
s["cost_summary"]=_CBD.print_solution_wcost_data(result["1"]["result_mip"], result["1"]["s"], result["1"]["data"])
pdic=_CBD.problemOUTPUT_map_byTimeStep(result["1"])
PlotlyJS.plot(pdic["trace0"], pdic["layout"])

################## Run zonal Formulation ###################
#zonal data setup
s_z["home_market"]=[[9,10,11,12]]
mn_data_z, data_z, s_z = _CBD.data_setup(s_z);#Build data structure for given options
#display inputs
_CBD.problemINPUT_map(data, s)
_CBD.problemINPUT_mapNTCs(data, s)
#solve zonal market
result_z=_CBD.zonal_market_main(mn_data_z, data_z, s_z)
#display results
s_z["cost_summary"]=_CBD.print_solution_wcost_data(result_z["1"]["result_mip"], result_z["1"]["s"], result_z["1"]["data"])
pdic=_CBD.problemOUTPUT_map_byTimeStep(result_z["1"])
PlotlyJS.plot(pdic["trace0"], pdic["layout"])
