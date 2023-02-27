function optimise_β(
    ::Local,
    loss,
    u0)
    optfunc = Optimization.OptimizationFunction(loss, Optimization.AutoForwardDiff())
    optprob = Optimization.OptimizationProblem(optfunc,u0)
    opt_param = Optimization.solve(optprob, BFGS())
    return opt_param.u[1]
end

function optimise_β(
        ::Global,
        loss,
        u0)
        optfunc = Optimization.OptimizationFunction(loss, Optimization.AutoForwardDiff())
        optprob = Optimization.OptimizationProblem(optfunc,u0,lb = [0.0], ub = [10.0])
        opt_param = Optimization.solve(
            optprob, 
            BBO_adaptive_de_rand_1_bin_radiuslimited(), 
            maxiters = 100000,
            maxtime = 1000.0)
        return opt_param.u[1]
end


function optimise_K(
    ::Global,
    loss,
    u0,
    bounds::Tuple)
    optfunc = Optimization.OptimizationFunction(loss, Optimization.AutoForwardDiff())
    optprob = Optimization.OptimizationProblem(optfunc,u0,lb = bounds[1], ub = bounds[2])
    opt_param = Optimization.solve(
        optprob, 
        BBO_adaptive_de_rand_1_bin_radiuslimited(), 
        maxiters = 100000,
        maxtime = 1000.0)
    return opt_param.u[1]
end





function first_analytical_solution(time_vec,N₀,r)
    fun(u) = test_compute_pop_ode_sol(time_vec,N₀,r,u) 
    return fun
end

function myloss(analytical_fun,data)
    loss(u,p=nothing) = sum(abs2.(analytical_fun(u) .- data))
    return loss
end




"""
    Take a dataframe generated from `convert_mat_to_df()`.
    Generate a line of best fit using the lm model from the GLM 
    package.
"""
function compute_lm_predict(df)
    d̂f = []
    R² = []
    for x ∈ names(df[:,Not([:radius,:comp])])
        model_formula = term(:comp) ~ term((x))
        model = lm(model_formula,df)         
        ŷ = predict(model, df)
        push!(R²,r2(model))
        push!(d̂f,ŷ)
    end
    predict_df = DataFrame(d̂f,:auto)

    return (predict_df = predict_df,R² = R²)
end

"""
    Take two vectors, `x` and `y` and compute the line of best fit
    manually. Return the line of `ŷ`.
"""
function fit_line(x,y)
    ybar = mean(y)
    xbar = mean(x)
    xbar_diff = x .- xbar
    ybar_diff = y .- ybar
    slope = sum(xbar_diff .* ybar_diff)/sum(xbar_diff .^ 2)
    intercept = ybar - slope*xbar
    y_line = intercept .+ slope .* x
    return y_line
end 
