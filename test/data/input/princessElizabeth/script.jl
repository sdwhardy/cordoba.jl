################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, OrderedCollections, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels

##################### Input parameters #################################
#1: 1 - 2 MVA: 43.339767 Length: 17 Cost: 423.60074 Status: 1
#2: 1 - 2 MVA: 32.599884 Length: 17 Cost: 272.2895 Status: 1
#3: 1 - 2 MVA: 21.733253 Length: 17 Cost: 185.18484 Status: 1
#4: 1 - 2 MVA: 10.866627 Length: 17 Cost: 94.53569 Status: 1

s = Dict(
"rt_ex"=>pwd()*"\\test\\data\\input\\princessElizabeth\\",#folder path if directly
"scenario_data_file"=>"C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_4EU.jld2",
################# temperal parameters #################
"test"=>true,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>1,
"scenario_names"=>["NT2025"],#["NT","DE","GA"]#,"NT2030","NT2040","DE2030","DE2040","GA2030","GA2040"
"k"=>1,#number of representative days modelled (24 hours per day)//#best for maintaining mean/max is k=6 2014, 2015
"res_years"=>["2014"],#Options: ["2012","2013","2014","2015","2016"]//#best for maintaining mean/max is k=6 2014, 2015
"scenario_years"=>["2020"],#Options: ["2020","2030","2040"]
################# Financial parameters ################
"dr"=>0.04,#discount rate
"yearly_investment"=>1000000,#max investment per modelling year
################ electrical parameters ################
"conv_lim_onshore"=>3000,#Max Converter size in MVA
"conv_lim_offshore"=>4000,#Max Converter size in MVA
"strg_lim_offshore"=>0.2,#Max offshore storage capacity
"strg_lim_onshore"=>10,#Max onshore storage capacity
"candidate_ics_ac"=>[1,3/4,1/2,2/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1],#DC Candidate Cable sizes (fraction of full MVA)[1,4/5,3/5,2/5]
################ collection circuit options ##############
"collection_circuit"=>true,
"no_crossings"=>true,
"collection_voltage"=>66,
"oss_nodes"=>[2],
"max_num_strings_per_oss"=>[15],
"max_num_of_branches_per_turbine"=>1,#1 consider only radial connections >1 branches at turbines 
"cable_losses"=>true,
#"cross_section_limit"=>5,
#"max_turbines_per_string"=>9,#not functional yet
#"no_loops"=>true,#not functional yet
################## optimization/solver setup options ###################
"relax_problem" => false,#binaries->continuous variables
"corridor_limit" => true,#limit cables in parallel?
"TimeLimit" => 75000,#244800,#solver max time in seconds
"MIPGap"=>1e-4,#max gap between MIP and convex solution 
"PoolSearchMode" => 0,#0-single solution, 1- poolsolutions of random quality, 2- poolsolutions of highest quality 
"PoolSolutions" => 1)#number of solutions to find
s=_CBD.hidden_settings(s)

################## Run nodal Formulation ###################
#nodal data setup
s["home_market"]=[]
mn_data, data, s = _CBD.data_setup(s);#9.653269610441e+04

#display inputs
#_CBD.problemINPUT_map(data, s)
#_CBD.problemINPUT_mapNTCs(data, s)
#solve nodal market
result = _CBD.collection_circuit_main(mn_data, data, s)
#display results
pdic=_CBD.problemMIP_OUTPUT_map_byTimeStep(result["1"])
PlotlyJS.plot(pdic["trace0"], pdic["layout"])


FileIO.save("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\collection_circuit_elizabeth_phase2.jld2",result)#09gap was good one


result["1"]["result_mip"]["1"]["branchdc_ne"]["1"]
mn_data["nw"]["1"]["ne_branch"]["346"]["br_r"]=mn_data["nw"]["1"]["ne_branch"]["345"]["br_r"]
mn_data["nw"]["1"]["ne_branch"]["346"]["br_x"]=mn_data["nw"]["1"]["ne_branch"]["345"]["br_x"]

for (i,br) in mn_data["nw"]["1"]["storage"]["1"]
    #if (br["built"]>0.9)
   #     println(i," ",mn_data["nw"]["1"]["ne_branch"][i]["f_bus"]," ",mn_data["nw"]["1"]["ne_branch"][i]["t_bus"]," ",mn_data["nw"]["1"]["ne_branch"][i]["rate_a"])
    #end
    println(br["e_absmax"])
end

result["1"]["result_mip"]["1"]["convdc"]