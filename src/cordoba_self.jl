module Cordoba_self
import JuMP
import JSON
import CSV
import XLSX
import Memento
import Dates
import DataFrames
import OrderedCollections
import FileIO
import JLD2
import Geodesy
import PowerModels;const _PM = PowerModels
import PowerModelsACDC;const _PMACDC = PowerModelsACDC
import FlexPlan;const _FP = FlexPlan
import InfrastructureModels;const _IM = InfrastructureModels
import MathOptInterface;const _MOI = MathOptInterface

include("prob/cordoba_acdc_wf_strg.jl")
include("core/objective.jl")
include("core/constraints.jl")
include("core/admm.jl")
include("core/admm_constraints.jl")
include("core/admm_objective.jl")
include("core/storage.jl")
include("core/variables.jl")
include("io/profile_data.jl")
include("io/functions.jl")
include("io/print_m_file.jl")

##################################################### NOTE ####################################
#economics is an external package under development
#AC_cbl(mva,km) DC_cbl(mva, km) must be defined somewhere if not using external package
include(pwd()*"..//..//..//packages//economics//src//economics.jl");const _ECO = economics
###############################################################################################
end # module
