#used
function storage_costs(data,year::Int64=2021)
    if (year==2021)
        csts=[1888,944,472,236,118,59,29.5,14.75,7.375]
    elseif (year==2030)

        x=1
        csts=[1120/x,560/x,280/x,140/x,70/x,35/x,17.5/x,8.75/x,4.375/x]

    else
        println("No battery cost data for specified year defaulting to 2021.")
        csts=[1888,944,472,236,118,59,29.5,14.75,7.375]
    end

    [data["ne_storage"][string(i)]["eq_cost"]=cst for (i,cst) in enumerate(csts)]
    for (i,cst) in enumerate(csts);if (data["ne_storage"][string(i)]["on_off"]==0); cst=0;end;data["ne_storage"][string(i)]["inst_cost"]=cst; end
    return data
end
