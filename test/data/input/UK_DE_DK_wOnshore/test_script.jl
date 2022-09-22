################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

function main_test()
    s = Dict(
    "rt_ex"=>pwd()*"\\data\\input\\UK_DE_DK\\",#folder path
    "scenario_data_file"=>"C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\scenario_data_for_UKBEDEDK.jld2",
    ################# temperal parameters #################
    "test"=>true,#if true smallest (2 hour) problem variation is built for testing
    "scenario_planning_horizon"=>30,
    "scenario_names"=>["NT"],#["NT","DE","GA"]
    "k"=>2,#number of representative days modelled (24 hours per day)//#best for maintaining mean/max is k=6 2014, 2015
    "res_years"=>["2014"],#Options: ["2012","2013","2014","2015","2016"]//#best for maintaining mean/max is k=6 2014, 2015
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
    "candidate_ics_ac"=>[1,1/2],#AC Candidate Cable sizes (fraction of full MVA)
    "candidate_ics_dc"=>[1,1/2],#DC Candidate Cable sizes (fraction of full MVA)
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
    result_mip_n=_CBD.social_welfare(s)
    s, result_mip_n, data, mn_data=_CBD.summarize_in_s(result_mip_n);
    s["home_market"]=[[2,4]]
    result_mip_z=_CBD.social_welfare(s)
    s, result_mip_z, data, mn_data=_CBD.summarize_in_s(result_mip_z);
    return result_mip_n["objective"]+result_mip_z["objective"]
end