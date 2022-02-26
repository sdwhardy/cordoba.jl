using Plots, Geodesy, XLSX, CSV, JLD2, FileIO
#UTM zones: UK 30N, NL 31N, Sweden 32N
#verification that using 31N for all will give proper distances
#tronstad=(58.66283058,6.716991259)
#great_island=(52.38114073, -7.017679964)
#Tronstad 32: 449833.9533331686, 31: 449833.95333316666, 30: 449833.95333316876
#Great island 30: 761696.9547920243, 31: 761696.9547920241, 32: 761696.9547920243
interconnectors_df = DataFrame(XLSX.readtable("test/data/input/entso_e_projects/infrastructure.xlsx", "interconnectors")...)
owpps_df = DataFrame(XLSX.readtable("test/data/input/entso_e_projects/infrastructure.xlsx", "owpps")...)
interconnectors_vec=DataFrame("f_country"=>String[],"t_country"=>String[],"f_xy"=>Tuple[],"t_xy"=>Tuple[],"km"=>Float64[])#MWh_up,EUR_up,MWh_dwn,EUR_dwn

centre_of_north_sea_lla=(54.98347204726825,3.6677885444483382)#(lat,long)
centre_of_north_sea_xy=utm_gps2xy(centre_of_north_sea_lla)

for node=1:1:length(first(first(froms_tos)))
    from_xy=utm_gps2xy((interconnectors_df["from_lattitude"][node],interconnectors_df["from_longitude"][node]))
    from_x=Float64((from_xy.x-centre_of_north_sea_xy.x)/1000);from_y=Float64((from_xy.y-centre_of_north_sea_xy.y)/1000)
    to_xy=utm_gps2xy((interconnectors_df["to_lattitude"][node],interconnectors_df["to_Longitude"][node]))
    to_x=Float64((to_xy.x-centre_of_north_sea_xy.x)/1000);to_y=Float64((to_xy.y-centre_of_north_sea_xy.y)/1000)
    push!(interconnectors_vec,(interconnectors_df["from_country"][node],interconnectors_df["to_country"][node],(from_x,from_y),(to_x,to_y),Float64(interconnectors_df["total_route_length_(km)"][node])))
end


function utm_gps2xy(lla,north_south::Bool=true,zone_utm::Int64=31)
    utm_desired = UTMfromLLA(zone_utm, north_south, wgs84)#sets UTM zone
    utm = utm_desired(LLA(first(lla),last(lla)))#coverts to cartesian
    return utm
end

plotly()
p=plot()
width=750
height=500
p=plot(size = (width, height))
#plot!(p,[ps],[cbls_costs],xticks = 0:200:800,xlims=(0,800),ylims=(0,3),color = :black,label="",xaxis = ("MW", font(20, "Courier")),yaxis = ("ME", font(20, "Courier")))
for node=1:1:length(interconnectors_vec["f_xy"])-8
plot!(p,[first(interconnectors_vec["f_xy"][node]),first(interconnectors_vec["t_xy"][node])],[last(interconnectors_vec["f_xy"][node]),last(interconnectors_vec["t_xy"][node])],color = :red,markersize=12)
end;gui()
