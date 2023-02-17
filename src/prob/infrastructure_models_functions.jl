""
function optimize_model!(aim::_IM.AbstractInfrastructureModel; relax_integrality=false, optimizer=nothing, solution_processors=[])
    start_time = time()

    if relax_integrality
        JuMP.relax_integrality(aim.model)
    end

    if JuMP.mode(aim.model) != JuMP.DIRECT && optimizer !== nothing
        if JuMP.backend(aim.model).optimizer === nothing
            JuMP.set_optimizer(aim.model, optimizer)
        else
            Memento.warn(_PM._LOGGER, "Model already contains optimizer, cannot use optimizer specified in `optimize_model!`")
        end
    end

    if JuMP.mode(aim.model) != JuMP.DIRECT && JuMP.backend(aim.model).optimizer === nothing
        Memento.error(_PM._LOGGER, "No optimizer specified in `optimize_model!` or the given JuMP model.")
    end

    _, solve_time, solve_bytes_alloc, sec_in_gc = @timed JuMP.optimize!(aim.model)

    try
        solve_time = JuMP.solve_time(aim.model)
    catch
        Memento.warn(_PM._LOGGER, "The given optimizer does not provide the SolveTime() attribute, falling back on @timed.  This is not a rigorous timing value.");
    end
    
    Memento.debug(_PM._LOGGER, "JuMP model optimize time: $(time() - start_time)")

    start_time = time()
    result = build_result(aim, solve_time; solution_processors=solution_processors)
    Memento.debug(_LOGGER, "solution build time: $(time() - start_time)")

    aim.solution = result["solution"]

    return result
end


function build_result(aim::_IM.AbstractInfrastructureModel, solve_time; solution_processors=[])
    # try-catch is needed until solvers reliably support ResultCount()
    result_count = 1
    try
        result_count = JuMP.result_count(aim.model)
    catch
        Memento.warn(_PM._LOGGER, "the given optimizer does not provide the ResultCount() attribute, assuming the solver returned a solution which may be incorrect.");
    end

    solution = Dict{String,Any}()

    if result_count > 0
        solution = build_solution(aim, post_processors=solution_processors)
    else
        Memento.warn(_PM._LOGGER, "model has no results, solution cannot be built")
    end

    result = Dict{String,Any}(
        "optimizer" => JuMP.solver_name(aim.model),
        "termination_status" => JuMP.termination_status(aim.model),
        "primal_status" => JuMP.primal_status(aim.model),
        "dual_status" => JuMP.dual_status(aim.model),
        "objective" => _guard_objective_value(aim.model),
        "objective_lb" => _guard_objective_bound(aim.model),
        "solve_time" => solve_time,
        "solution" => solution,
    )

    return result
end

""
function _guard_objective_value(model)
    obj_val = NaN

    try
        obj_val = JuMP.objective_value(model)
    catch
    end

    return obj_val
end


""
function _guard_objective_bound(model)
    obj_lb = -Inf

    try
        obj_lb = JuMP.objective_bound(model)
    catch
    end

    return obj_lb
end


""
function build_solution(aim::_IM.AbstractInfrastructureModel; post_processors=[])
    sol = Dict{String, Any}("it" => Dict{String, Any}())
    sol["multiinfrastructure"] = true

    for it in it_ids(aim)
        sol["it"][string(it)] = build_solution_values(aim.sol[:it][it])
        sol["it"][string(it)]["multinetwork"] = true
    end

    _IM.solution_preprocessor(aim, sol)

    for post_processor in post_processors
        _IM.post_processor(aim, sol)
    end

    for it in it_ids(aim)
        it_str = string(it)
        data_it = _IM.ismultiinfrastructure(aim) ? aim.data["it"][it_str] : aim.data

        if _IM.ismultinetwork(data_it)
            sol["it"][it_str]["multinetwork"] = true
        else
            for (k, v) in sol["it"][it_str]["nw"]["$(nw_id_default)"]
                sol["it"][it_str][k] = v
            end

            sol["it"][it_str]["multinetwork"] = false
            delete!(sol["it"][it_str], "nw")
        end

        if !_IM.ismultiinfrastructure(aim)
            for (k, v) in sol["it"][it_str]
                sol[k] = v
            end

            delete!(sol["it"], it_str)
        end
    end

    if !_IM.ismultiinfrastructure(aim)
        sol["multiinfrastructure"] = false
        delete!(sol, "it")
    end

    return sol
end

""
function build_solution_values(var::Dict)
    sol = Dict{String, Any}()

    for (key, val) in var
        sol[string(key)] = build_solution_values(val)
    end

    return sol
end

""
function build_solution_values(var::Array{<:Any,1})
    return [build_solution_values(val) for val in var]
end

""
function build_solution_values(var::Array{<:Any,2})
    return [build_solution_values(var[i, j]) for i in 1:size(var, 1), j in 1:size(var, 2)]
end

""
function build_solution_values(var::Number)
    return var
end

""
function build_solution_values(var::JuMP.VariableRef)
    return JuMP.value(var)
end

""
function build_solution_values(var::JuMP.GenericAffExpr)
    return JuMP.value(var)
end

""
function build_solution_values(var::JuMP.GenericQuadExpr)
    return JuMP.value(var)
end

""
function build_solution_values(var::JuMP.NonlinearExpression)
    return JuMP.value(var)
end

""
function build_solution_values(var::JuMP.ConstraintRef)
    return JuMP.dual(var)
end

""
function build_solution_values(var::Any)
    Memento.warn(_PM._LOGGER, "build_solution_values found unknown type $(typeof(var))")
    return var
end
### Helper functions for working with AbstractInfrastructureModels
function it_ids(aim::_IM.AbstractInfrastructureModel) 
    return keys(aim.ref[:it])
end