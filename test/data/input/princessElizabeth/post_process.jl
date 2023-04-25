################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS, XLSX
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections


result=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\collection_circuit_elizabeth_phase2.jld2")#09gap was good one
result=FileIO.load("C:\\Users\\shardy\\Documents\\julia\\times_series_input_large_files\\onshore_grid\\collection_circuit_elizabeth_phase1.jld2")#09gap was good one

pdic=_CBD.problemMIP_OUTPUT_map_byTimeStep(result["1"])
PlotlyJS.plot(pdic["trace0"], pdic["layout"])