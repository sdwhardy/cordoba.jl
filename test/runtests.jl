################## loads external packages ##############################
using Test
include("data/input/UK_DE_DK/test_script.jl")

@test isapprox(main_test(), -54821.45904211383; atol = 0.1) 
