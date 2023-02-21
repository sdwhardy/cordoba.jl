################## loads external packages ##############################
using Gurobi, JuMP, DataFrames, FileIO, CSV, Dates, PlotlyJS, XLSX
import Cordoba_self; const _CBD = Cordoba_self#Cordoba package backend - under development
import PowerModelsACDC; const _PMACDC = PowerModelsACDC
import PowerModels; const _PM = PowerModels
using OrderedCollections


nodes_n = DataFrames.DataFrame(XLSX.readtable("C:\\Users\\shardy\\Documents\\julia\\packages\\cordoba\\test\\data\\input\\ronne_bank\\ronne_bank_gps.xlsx", "north")...)
nodes_s = DataFrames.DataFrame(XLSX.readtable("C:\\Users\\shardy\\Documents\\julia\\packages\\cordoba\\test\\data\\input\\ronne_bank\\ronne_bank_gps.xlsx", "south")...)

lat_n=[];long_n=[]
for _row in eachrow(nodes_n)
    _lla=_CBD.utm_xy2gps((_row[:_x],_row[:_y]), true, 33)
    push!(lat_n,_lla.lat)
    push!(long_n,_lla.lon)
end
_GPS_n=DataFrames.DataFrame(:lattitude=>lat_n,:longitude=>long_n)

lat_s=[];long_s=[]
for _row in eachrow(nodes_s)
    _lla=_CBD.utm_xy2gps((_row[:_x],_row[:_y]), true, 33)
    push!(lat_s,_lla.lat)
    push!(long_s,_lla.lon)
end
_GPS_s=DataFrames.DataFrame(:lattitude=>lat_s,:longitude=>long_s)
CSV.write("C:\\Users\\shardy\\Documents\\julia\\packages\\cordoba\\test\\data\\input\\ronne_bank\\ronne_bank_north.csv", _GPS_n)
CSV.write("C:\\Users\\shardy\\Documents\\julia\\packages\\cordoba\\test\\data\\input\\ronne_bank\\ronne_bank_south.csv", _GPS_s)
