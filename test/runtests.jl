################## loads external packages ##############################
using Test

include("data/input/test/test_script_full.jl")
##################### File parameters #################################

@test isapprox(main_test(), -1.570149705373684e7; atol = 0.1) 
