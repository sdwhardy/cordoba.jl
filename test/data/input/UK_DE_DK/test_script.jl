############################## DO NOT CHANGE TEST SCRIPT ###################
function main_test()
    ########################################################################
    rez=0.0;
    ##################### Input parameters #################################
    rt_ex=pwd()*"\\data\\input\\UK_DE_DK\\"#folder path
    argz = Dict(
    "test"=>true,#if true smallest (2 hour) problem variation is built for testing
    "scenario_planning_horizon"=>30,
    "scenario_names"=>["EU17"],#Options: ["EU17","EU18","EU19","EU20","ST17","ST18","ST19","ST20","DG17","DG18","DG19","DG20"]
    "k"=>2,#number of representative days modelled (24 hours per day)//Must add clustered time series for each k
    "scenario_years"=>["2020","2030","2040"],#Options: ["2020","2030","2040"]
    "owpp_mva"=>[4000],#mva of wf in MVA
    "conv_lim"=>4000,#Max Converter size in MVA
    "candidate_ics_ac"=>[6/5,1,4/5,3/5,2/5,1/5],#AC Candidate Cable sizes (fraction of full MVA)
    "candidate_ics_dc"=>[6/5,1,4/5,3/5,2/5,1/5],#DC Candidate Cable sizes (fraction of full MVA)
    "dr"=>0.04,#discount rate
    "yearly_investment"=>1000000)

    ################## optimization/solver setup options ###################
    s = Dict("output" => Dict("branch_flows" => false),
    "home_market"=>[2,4],#nodes within Zonal market
    "balancing_reserve"=>0.3,#zonal market must be defined to have any effect
    "AC"=>"1",#0=false, 1=true
    "eps"=>0.1,#admm residual (100kW)
    "relax_problem" => false,
    "conv_losses_mp" => false,
    "process_data_internally" => false,
    "corridor_limit" => true,
    "strg_lim_offshore"=>0.2,
    "strg_lim_onshore"=>10,
    "max_invest_per_year"=>_CBD.max_invest_per_year(argz))
    ########################################################################

    ################## Run MIP Formulation ###################
    mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)#Build data structure for given options
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
    #with Home market in Germany
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    rez=rez+result_mip["objective"]
    #_CBD.print_solution_data(result_mip, data, argz)#print solution
    #with OBZ
    s["home_market"]=[]#nodes within Zonal market
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    #_CBD.print_solution_data(result_mip, data, argz)#print solution
    rez=rez+result_mip["objective"]
    ################## Run Convex Formulation ################
    s["relax_problem"]=true
    mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)#Build data structure for given options
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)#select solver
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    #_CBD.print_solution_data(result_mip, data, argz)#print convex solution
    rez=rez+result_mip["objective"]
    mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)#Convert convex to candicate cables
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)#In large problems a larger MIPgap (1-5%?) may be desirable
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
    rez=rez+result_mip["objective"]
    #_CBD.print_solution_data(result_mip, data, argz)#print MIP solution


    ################## Run ADMM Formulation ################
    s["relax_problem"]=true
    mn_data, data, argz, s = _CBD.main_ACDC_wstrg(rt_ex,argz, s)#Build data structure for given options
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 0)
    result_mip = _CBD.admm_4_AjAwAgAuAo_main(mn_data, gurobi, s)#Solve problem
    rez=rez+result_mip["objective"]
    #_CBD.print_solution_data(result_mip, data, argz)#print Convex (ADMM) solution
    mn_data, data, s = _CBD.convex2mip(result_mip, data, mn_data, s)#Convert convex to candicate cables
    gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1, "MIPGap" => 1e-4)#In large problems a larger MIPgap (1-5%?) may be desirable
    result_mip = _CBD.cordoba_acdc_wf_strg(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)#Solve problem
    #_CBD.print_solution_data(result_mip, data, argz)#print MIP solution
    rez=rez+result_mip["objective"]
    return rez
end