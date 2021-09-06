module wind_module
wind_profs=Dict()
function save_wind4_module(wnd,name)
        push!(wind_profs,(name=>wnd))
end
end
