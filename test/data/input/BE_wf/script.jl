################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, OrderedCollections, FileIO
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels

##################### Input parameters #################################
rt_ex=pwd()*"\\test\\data\\input\\BE_wf\\"#folder path
argz = Dict(
"test"=>false,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
"scenario_names"=>["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"],#Options: ["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
"k"=>5,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
"owpp_mva"=>[4000],#mva of wf in MVA
"conv_lim"=>4000,#Max Converter size in MVA
"candidate_ics_ac"=>[1/10],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1],#,4/5,3/5,1/2],#DC Candidate Cable sizes (fraction of full MVA)
"dr"=>0.04,#discount rate
"yearly_investment"=>1000000)

################## optimization/solver setup options ###################
s = Dict("output" => Dict("branch_flows" => false),
"home_market"=>[],#nodes within Zonal market
"balancing_reserve"=>0.3,#zonal market must be defined to have any effect
"AC"=>"0",#0=false, 1=true
"eps"=>0.0001,#admm residual (100kW)
"relax_problem" => false,
"conv_losses_mp" => false,
"process_data_internally" => false,
"corridor_limit" => true,
"strg_lim_offshore"=>0.2,
"strg_lim_onshore"=>10,
"max_invest_per_year"=>_CBD.max_invest_per_year(argz))
########################################################################
#DE offshore wind installed capacity: 4131(2017), 5051(2018), 6393(2019), 7504(2020)
#scenario_data, ls = _CBD.load_time_series(rt_ex,argz)
  
#file = rt_ex*"topology.m"
#data = PowerModels.parse_file(file)
################## Run MIP Formulation ###################
mn_data, data, argz, s = _CBD.main_ACDC_chandra(rt_ex,argz, s);#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#, "MIPGap" => 3e-2)#select solver
#with OBZ
s["home_market"]=[]#nodes within Zonal market
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print solution
FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\chandra\\result_Bewf_k5.jld2",result_mip)
   

###########################################
#update time series
#installed_capacity_DE=Dict("17"=>4131,"18"=>5051,"19"=>6393,"20"=>7504)
#=scenario_data, ls = _CBD.load_time_series(rt_ex,argz)
installed_capacity_DE=Dict("17"=>4668,"18"=>5153,"19"=>6624,"20"=>6980)
installed_capacity_BE=Dict("17"=>1.020432,"18"=>1.01139,"19"=>1.013661,"20"=>1.001851)

for key1 in keys(scenario_data)
    _year=key1[3:4]
    for key2 in keys(scenario_data[key1])
    #    println(string(key1)*" "*string(key2)*" "*string(maximum(scenario_data[key1][key2][!,:Wnd_MWhDE])))
    #scenario_data[key1][key2][!,:Wnd_MWhBE]=scenario_data[key1][key2][!,:Wnd_MWhBE].*installed_capacity_BE[_year]    
    scenario_data[key1][key2][!,:Wnd_MWhDE]=scenario_data[key1][key2][!,:Wnd_MWhDE]./installed_capacity_DE[_year]
    end
end

for key1 in keys(scenario_data)
    _year=key1[3:4]
    for key2 in keys(scenario_data[key1])
        println(string(key1)*" "*string(key2)*" "*string(maximum(scenario_data[key1][key2][!,:Wnd_MWhBE])))
       # println(string(key1)*" "*string(key2)*" "*string(maximum(scenario_data[key1][key2][!,:Wnd_MWhDE])))
   #     scenario_data[key1][key2][!,:Wnd_MWhDE]=scenario_data[key1][key2][!,:Wnd_MWhDE]./installed_capacity_DE[_year]
    end
end

FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\time_series_k5.jld2", scenario_data)
  =#