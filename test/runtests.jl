################## loads external packages ##############################
using Test

include("data/input/test/test_script_auto.jl")
##################### File parameters #################################

@test isapprox(main_test(), -1.5776553277055133e7; atol = 0.1) 
