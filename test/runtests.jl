################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, Test
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
include("data/input/UK_DE_DK/test_script.jl")

@test isapprox(main_test(), -54962.66562984746; atol = 0.1) 
