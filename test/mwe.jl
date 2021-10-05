using JuMP, Gurobi
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
import InfrastructureModels; const _IM = InfrastructureModels

casename = "mwe"
file = "./test/data/input/$casename.m"
data = _PM.parse_file(file)#load data in PM format
_PMACDC.process_additional_data!(data)#add extra DC model data
#select solver
gurobi = JuMP.optimizer_with_attributes(Gurobi.Optimizer,"OutputFlag" => 1)

#settings
s = Dict("output" => Dict("branch_flows" => true), "conv_losses_mp" => false, "process_data_internally" => false)

resultACDC = _PMACDC.run_tnepopf(data, _PM.DCPPowerModel, gurobi, setting = s)
