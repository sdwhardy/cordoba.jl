################## loads external packages ##############################
using Ipopt, Juniper, JuMP, Cbc, Gurobi, Cbc, XLSX, DataFrames, Dates, CSV
import cordoba; const _CBD = cordoba#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import FlexPlan; const _FP = FlexPlan
import InfrastructureModels; const _IM = InfrastructureModels
include("../test/data/conv_spec.jl")

################### Input data for script_test.jl ########################
owpp_mva=4000#mva of wf (in de)
#interconnectors format: (mva,(km,type)) type: -1:both ends onshore, 0:both ends offshore, 1:onshore/offshore
#must be in same order as .m file ie uk(1)-de(2), uk(1)-dk(3), de(2)-dk(3), wf(4)-uk(1), wf(4)-de(2), wf(4)-dk(3)
ics=[(2000,(550,-1)),(2000,(760,-1)),(2000,(250,-1)),(2000,(470,1)),(2000,(145,1)),(2000,(246,1))];
#location of nodes
zs=["UK","DE","DK","DE"]#must be in same order as .m file gens
candidate_ics=[1,4/5,3/5,1/2]#Candidate Cables
infinite_grid=sum(first.(ics))+owpp_mva#ensures enough generation and consumption in all markets

################## reads .m input file name ######################
casename = "test"
file = "./test/data/input/$casename.m"
data = _PM.parse_file(file)#load data in PM format
data=_CBD.additional_candidatesICS(data,candidate_ics,ics)#adds additional candidates

####################### These functions are not needed if data is directly entered into .m file but this is a recipe for disaster
#essentially it is just duplicating connections in branchdc_ne and adding real cable values for power flow
#Cost of a candidate is the converter/plarforms/transformers and cables all lumped together
#################### Calculates cables for DC lines
for (i,bdc) in data["branchdc_ne"]
data["branchdc_ne"][i]=_CBD.candidateIC_cost(bdc);end
data["branchdc_ne"]=_CBD.unique_candidateIC(data["branchdc_ne"])#keep only unique candidates
################################################################################

#################### PowerModelsACDC settings ############################
_PMACDC.process_additional_data!(data)#add extra DC model data
converter_parameters_rxb(data)#adjust converter parameters

#################### Multi-period input parameters #######################
n=number_of_hours = 7151 # Number of time points in DE
scenario = Dict{String, Any}("hours" => number_of_hours, "sc_years" => Dict{String, Any}())
scenario["sc_years"]["1"] = Dict{String, Any}()
scenario["sc_years"]["1"]["year"] = 2020#year of data
scenario["sc_years"]["1"]["start"] = Dates.datetime2epochms(DateTime(scenario["sc_years"]["1"]["year"]))   #seconds from year 0 (2020 = 63745056000000)
scenario["sc_years"]["1"]["probability"] = 1
scenario["planning_horizon"] = 30 # in years, to scale generation cost

##################### load time series data ##############################
#time series data files are UKdata.csv, DEdata.csv, DKdata.csv
data, zs_data = _CBD.get_profile_data_sets_mesh(unique(zs),data, n, scenario)

#set problem dimension
n=number_of_hours = length(zs_data.time_stamp)
dim = scenario["hours"] * length(data["scenario"])

 # create a dictionary to pass time series data to data dictionary
extradata,data = _CBD.create_profile_sets_mesh(dim, data, zs_data, zs, infinite_grid, owpp_mva)

# Scale cost data
_FP.scale_cost_data_cordoba!(data, scenario)
_FP.scale_cost_data_cordoba!(extradata, scenario)

# Create data dictionary where time series data is included at the right place
mn_data = _PMACDC.multinetwork_data(data, extradata, Set{String}(["source_type","scenario","name", "source_version", "per_unit"]))

#select solver
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

#optimation settings
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)

#run optimization
resultACDC = _CBD.cordoba_mp_acdctnepopf(mn_data, _PM.DCPPowerModel, gurobi, multinetwork=true; setting = s)
