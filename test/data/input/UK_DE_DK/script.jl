################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, OrderedCollections
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels

##################### Input parameters #################################
rt_ex=pwd()*"\\test\\data\\input\\UK_DE_DK\\"#folder path
argz = Dict(
"test"=>false,#if true smallest (2 hour) problem variation is built for testing
"scenario_planning_horizon"=>30,
"scenario_names"=>["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"],#Options: ["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
"k"=>2,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k
"scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
"owpp_mva"=>[4000],#mva of wf in MVA
"conv_lim"=>4000,#Max Converter size in MVA
"candidate_ics_ac"=>[1/5],#AC Candidate Cable sizes (fraction of full MVA)
"candidate_ics_dc"=>[1,4/5,3/5,1/2],#DC Candidate Cable sizes (fraction of full MVA)
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
#scenario_data, ls = _CBD.load_time_series(rt_ex,argz)
#file = rt_ex*"topology.m"
#data = PowerModels.parse_file(file)
################## Run MIP Formulation ###################
mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s);#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 5e-2)#select solver
#with Home market in Germany
#result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
#_CBD.print_solution_data(result_mip, data, argz)#print solution
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
_CBD.print_solution_data(result_mip, data, argz)#print MIP solution

################## Run ADMM Formulation ################
s["relax_problem"]=true
mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)#Build data structure for given options
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)
result_mip = _CBD.admm_4_AjAwAgAuAo_main(mn_data, gurobi, s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print Convex (ADMM) solution
mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)#Convert convex to candicate cables
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)#In large problems a larger MIPgap (1-5%?) may be desirable
result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
_CBD.print_solution_data(result_mip, data, argz)#print MIP solution

bens=0
for (nw_i,nw) in sort(OrderedCollections.OrderedDict(result_mip["solution"]["nw"]), by=x->parse(Int64,x))
   
   if (issubset([nw_i],string.(year2_keys))) 
     for (g_i,g) in nw["gen"]
        bens=bens+g["pg"]#*mn_data["nw"][nw_i]["gen"][g_i]["cost"][1]
        #println(mn_data["nw"][nw_i]["gen"][g_i]["cost"][1]) 
        end
    end 
end
bens/12

    bens/12
    mn_data["nw"]["1"]
    year1_keys=[];year2_keys=[];year3_keys=[]
    hours_per_timeline=length(mn_data["scenario"]["1"])
    2*hours_per_timeline/3
    for _s in sort(OrderedCollections.OrderedDict(mn_data["scenario"]), by=x->parse(Int64,x))
        length(first(_s))
        for _k in sort(OrderedCollections.OrderedDict(last(_s)), by=x->parse(Int64,x))
            if (parse(Int64,first(_k))<=hours_per_timeline/3)
                push!(year1_keys, last(_k))
            elseif (parse(Int64,first(_k))<=2*hours_per_timeline/3)
                push!(year2_keys, last(_k))
            else
                push!(year3_keys, last(_k))
            end
        end
    end
    println(year3_keys)


    average(scenario_data["EU20"]["2020"][!,:Wnd_MWhDE])

    
function main_ACDC_wstrg(rt_ex,argz, s)
    ################# Load topology files ###################################
    _CBD.topology_df(rt_ex, s["relax_problem"], s["AC"])#creates .m file
    data, ics_ac, ics_dc, nodes = _CBD.filter_mfile_cables(rt_ex)#loads resulting topology and filters for candidate cables
    #ics_dc=[(4000,550),(4000,760),(4000,250),(4000,443),(4000,212),(4000,311)]
    ############### defines size and market of genz and wfs ###################
    infinite_grid, genz, wfz, markets = genz_n_wfs(argz["owpp_mva"],nodes,data["baseMVA"])
    push!(argz,"genz"=>genz)
    push!(argz,"wfz"=>wfz)
    #################### Calculates cable options for AC lines
    data=AC_cable_options(data,argz["candidate_ics_ac"],ics_ac,data["baseMVA"])
    print_topology_data_AC(data,markets)#print to verify
    #################### Calculates cable options for DC lines
    data=DC_cable_options(data,argz["candidate_ics_dc"],ics_dc,data["baseMVA"])
    additional_params_PMACDC(data)
    print_topology_data_DC(data,markets)#print to verify
    ##################### load time series data ##############################
    scenario_data, ls = load_time_series(rt_ex,argz)
    push!(argz,"ls"=>ls)
    ##################### multi period setup #################################
    mn_data = multi_period_setup(ls, scenario_data, data, markets, infinite_grid, argz)
    s=update_settings(s, argz, data)
#   if (haskey(s,"home_market") && length(s["home_market"])>0)
#      mn_data=zonal_adjust(mn_data, s);end
    return  mn_data, data, argz, s
end