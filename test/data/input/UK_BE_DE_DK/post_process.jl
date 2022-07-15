################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections

results_nodal=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\nodal_market_k6_VOLL5000b.jld2")
results_allhm=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_allhm_k6_VOLL5000b_rc.jld2")
results_567=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\UK_BE_DE_DK\\zonal_results_567_k6_VOLL5000b_rc.jld2")