#=function converter_cost(data)
    for(c,conv) in data["convdc_ne"]
        conv["cost"] = conv["Pacmax"]*0.083 *100+ 28
        #display(conv["cost"])
    end
end=#

#
