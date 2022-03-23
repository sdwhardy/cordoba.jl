
###################################### Storage #########################################
function constraint_storage_losses(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    storage = _PM.ref(pm, nw, :storage, i)

    _PM.constraint_storage_losses(pm, nw, i, storage["storage_bus"], storage["r"], storage["x"], storage["p_loss"], storage["q_loss"])
end


function constraint_storage_thermal_limit(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    ps = _PM.var(pm, nw, :ps, i)
    e_absmax = _PM.var(pm, nw, :e_absmax, i)
    JuMP.@constraint(pm.model, ps-e_absmax/2  <= 0)
    JuMP.@constraint(pm.model, ps+e_absmax/4  >= 0)
end

function constraint_storage_excl_slack(pm::_PM.AbstractPowerModel, i::Int; nw::Int=pm.cnw)
    sc = _PM.var(pm, nw, :sc, i)
    sd = _PM.var(pm, nw, :sd, i)
    e_absmax = _PM.var(pm, nw, :e_absmax, i)

    JuMP.@constraint(pm.model, sc + sd <= e_absmax/2)
end
