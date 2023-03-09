################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, OrderedCollections, PlotlyJS, XLSX
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
##################### File parameters #################################

s = Dict(
    "rt_ex"=>pwd()*"\\test\\data\\input\\north_sea\\",#folder path if directly
    "scenario_data_file"=>"C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4EU.jld2",
    ################# temperal parameters #################
    "test"=>false,#if true smallest (2 hour) problem variation is built for testing
    "scenario_planning_horizon"=>30,
    "scenario_names"=>["NT2025","GA2030","GA2040"],#["NT","DE","GA"]#,"NT2030","NT2040","DE2030","DE2040","GA2030","GA2040"
    "k"=>4,#number of representative days modelled (24 hours per day)//#best for maintaining mean/max is k=6 2014, 2015
    "res_years"=>["2014","2015"],#Options: ["2012","2013","2014","2015","2016"]//#best for maintaining mean/max is k=6 2014, 2015
    "scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
    "dr"=>0.04,#discount rate
    "yearly_investment"=>10000000,
    ################ electrical parameters ################
    "conv_lim_onshore"=>36000,#Max Converter size in MVA
    "conv_lim_offshore"=>36000,#Max Converter size in MVA
    "strg_lim_offshore"=>0.2,
    "strg_lim_onshore"=>10,
    "candidate_ics_ac"=>[1,2,4,8],#AC Candidate Cable sizes (fraction of full MVA)
    "candidate_ics_dc"=>[1,2,4,8],#DC Candidate Cable sizes (fraction of full MVA)[1,4/5,3/5,2/5]
    ################## optimization/solver setup options ###################
    "relax_problem" => false,
    "corridor_limit" => false,
    "TimeLimit" => 320000,
    "MIPGap"=>1e-4, 
    "PoolSearchMode" => 0, 
    "PoolSolutions" => 1)
    s=_CBD.hidden_settings(s)



######################### Nodal market #########################
s["home_market"]=[]
mn_data, data, s = _CBD.data_setup(s);
#_CBD.problemINPUT_mapNTCs(data, s)
_CBD.problemINPUT_map(data, s)
@time result = _CBD.nodal_market_main(mn_data, data, s)#-3359431 -33899162 0.89%
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\NORTH_SEA_nodal_k4_full.jld2",result)#09gap was good one


######################### Zonal market #########################
#s["home_market"]=[[2,5],[3,6],[4,7]]
#s["home_market"]=[[9,10,11,12,13]]
s["home_market"]=[[4,11],[5,10],[6,12],[1,8,13],[3,9]]
mn_data, data, s = _CBD.data_setup(s);
@time result = _CBD.zonal_market_main(mn_data, data, s)#-3359431 -33899162 0.89%
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\HMD_results_NORTH_SEA_4G.jld2",result)
pdic=_CBD.problemOUTPUT_map_byTimeStep(result["3"])
PlotlyJS.plot(pdic["trace012"], pdic["layout"])
