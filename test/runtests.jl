################## loads external packages ##############################
using Test

include("data/input/UK_DE_DK/test_script.jl")
##################### File parameters #################################

@test isapprox(main_test(), -3.4539576650664983e6; atol = 0.1) 
