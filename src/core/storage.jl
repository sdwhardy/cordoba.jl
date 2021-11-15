#used
function storage_costs(data,year::Int64=2021)
    if (year==2021)
        #csts=[944,472,236,118,59,29.5,14.75,7.375]
        #csts=[472,236,118,59,29.5,14.75,7.375]
        #csts=[236,118,59,29.5,14.75,7.375]
        #csts=[118,59,29.5,14.75,7.375,3.6875,1.84375]
        csts=[1888,944,472,236,118,59,29.5,14.75,7.375]
    elseif (year==2030)
        #csts=[560,280,140,70,35,17.5,8.75,4.375]
        #csts=[280,140,70,35,17.5,8.75,4.375]
        x=1
        #csts=[70/x,35/x,17.5/x,8.75/x,4.375/x,2.1875/x,1.09375/x]
        #csts=[280/x,140/x,70/x,35/x,17.5/x,8.75/x,4.375/x]
        csts=[1120/x,560/x,280/x,140/x,70/x,35/x,17.5/x,8.75/x,4.375/x]
        #csts=[1120/x,4.375/x]
        #csts=[35/x,17.5/x,8.75/x]
        #csts=[1,1,1,1,1,1]
    else
        println("No battery cost data for specified year defaulting to 2021.")
        #csts=[944,472,236,118,59,29.5,14.75,7.375]
        #csts=[472,236,118,59,29.5,14.75,7.375]
        #csts=[236,118,59,29.5,14.75,7.375]
        #csts=[118,59,29.5,14.75,7.375,3.6875,1.84375]
        csts=[1888,944,472,236,118,59,29.5,14.75,7.375]
    end

    [data["ne_storage"][string(i)]["eq_cost"]=cst for (i,cst) in enumerate(csts)]
    for (i,cst) in enumerate(csts);if (data["ne_storage"][string(i)]["on_off"]==0); cst=0;end;data["ne_storage"][string(i)]["inst_cost"]=cst; end
    return data
end
