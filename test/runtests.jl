################## loads external packages ##############################
using Test

include("data/input/test/test_script_auto.jl")
##################### File parameters #################################

@test isapprox(main_test(), -1.5778239351652997e7; atol = 0.1) 
