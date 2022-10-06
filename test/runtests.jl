################## loads external packages ##############################
using Test

include("data/input/test/test_script_full.jl")
##################### File parameters #################################

@test isapprox(main_test(), -5.808978079364737e6; atol = 0.1) 
