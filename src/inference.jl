
mutable struct Paraminf
    settings::Array{Any,2} # DOE runs
    n_runs::Int64 # number of runs
    n_methods::Int64 # number of methods to compare
    max_servers::Int64 # max number of servers
    obs_max::Int64 # max observations available
    time_limit::Int64 # max number of sim reps
    ϵ::Float64 # convergence quality
    window::Int64 # observation window for convergence estimate
    window_detail::Int64 # observation window for error estimate
    step::Int64 # how many observations to skip while calculating convergence
    seed::Int64 # random number seed
end

# Build Distribution for Experiments
# Build appropriate distributions for given traffic intensities
# Input: number of servers, arrival dist, service dist, traffic intensity
# Output: arrival and service distribution
function builddist(c_true::Int64, aname::String, sname::String, rho::Float64)
    lam = 20
    norm_var = 3
    a = 1
    a_lam = 1
    a_mu = 1
    mu = lam / (c_true * rho)
    if aname == "Uniform"
        adist = Uniform(0, 2 * 1 / lam)
    elseif aname == "Exponential"
        adist = Exponential(1 / lam)
    elseif aname == "LogNormal"
        adist = LogNormal(log(1 / lam) - norm_var / 2, sqrt(norm_var))
    elseif aname == "Beta"
        adist = Beta(a_lam, (a_lam / (1 / lam)) - a_lam)
    else
        error("Arrival distribution must be Uniform, Exponential, LogNormal, or Beta")
    end
    if sname == "Uniform"
        sdist = Uniform(0, 2 * 1 / mu)
    elseif sname == "Exponential"
        sdist = Exponential(1 / mu)
    elseif sname == "LogNormal"
        sdist = LogNormal(log(1 / mu) - norm_var / 2, sqrt(norm_var))
    elseif sname == "Beta"
        sdist = Beta(a_mu, (a_mu / (1 / mu)) - a_mu)
    else
        error("Service distribution must be Uniform, Exponential, LogNormal, or Beta")
    end
    adist, sdist
end

# Convergence Estimator
# Calculate aprpoximate convergence
# Input: parametrs, queueing output data, number of servers
# Output: number of observations to approximately converge
function convergence(param_inf::Paraminf, queue_output::DataFrame, c_true::Int64)
    ind = 0
    max_servers = param_inf.max_servers
    lim = param_inf.max_servers + 1
    n_methods = param_inf.n_methods
    obs_max = param_inf.obs_max
    window = param_inf.window
    step = param_inf.step
    ϵ = param_inf.ϵ
    output = queue_output
    bool_true = falses(n_methods)
    bool_meas = falses(n_methods)
    conv_true = zeros(n_methods)
    conv_meas = zeros(n_methods)
    array_size = floor(Int64, (obs_max - 20) / step)
    ests_true = zeros(Int64, array_size, n_methods) # c estimates by obs by run
    ests_meas = zeros(Int64, array_size, n_methods) # c estimates by obs by run
    while lim <= (obs_max - step)
      # println(lim)
      lim = lim + step
      ind += 1
      # true departure
      out_order = sort(output[1:lim,:], :DepartOrder)[:ArriveOrder]
      ro = c_order(out_order)
      (rv, ru, VS) = c_var_unf(output[:Arrive][1:lim], output[:Depart][1:lim], max_servers)
      # time error
      out_order_meas = sort(output[1:lim,:], :DepartMeas)[:ArriveOrder]
      ro_meas = c_order(out_order_meas)
      (rv_meas, ru_meas, VS_meas) = c_var_unf(output[:Arrive][1:lim], output[:DepartMeas][1:lim], max_servers)
      # record results
      ests_true[ind, :] = [ro rv ru]
      ests_meas[ind, :] = [ro_meas rv_meas ru_meas]
      if ind >= window
        for j = 1:n_methods
          if !bool_true[j] && mean(abs.(ests_true[(ind - window + 1):(ind), j] - c_true)) < ϵ
            conv_true[j] = (ind - window + 1) * step
            bool_true[j] = true
          end
          if !bool_meas[j] && mean(abs.(ests_meas[(ind - window + 1):(ind), j]- c_true)) < ϵ
            conv_meas[j] = (ind - window + 1) * step
            bool_meas[j] = true
          end
        end
      end
    end
    (conv_true, conv_meas, ests_true, ests_meas)
end

# Estimation Error Estimator
# Calculate estimation error
# Input: parametrs, queue output data, obs limit, number of servers
# Output: average estimation error in estimation window
function esterror(param_inf::Paraminf, output::DataFrame, obs_limit::Int64, c_true::Int64)
    window_detail = param_inf.window_detail
    n_methods = param_inf.n_methods
    max_servers = param_inf.max_servers
    ind_l = Int64(obs_limit)
    ind_u = Int64(obs_limit  + window_detail - 1)
    detail_true = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
    detail_meas = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
    error_true = zeros(n_methods)
    error_meas = zeros(n_methods)
    ind = 0
    for lim = ind_l:ind_u
      ind += 1
      # true departure
      out_order = sort(output[1:lim,:], :DepartOrder)[:ArriveOrder]
      ro = c_order(out_order)
      (rv, ru, VS) = c_var_unf(output[:Arrive][1:lim], output[:Depart][1:lim], max_servers)
      # time error
      out_order_meas = sort(output[1:lim,:], :DepartMeas)[:ArriveOrder]
      ro_meas = c_order(out_order_meas)
      (rv_meas, ru_meas, VS_meas) = c_var_unf(output[:Arrive][1:lim], output[:DepartMeas][1:lim], max_servers)
      # record results
      detail_true[ind, :] = [ro rv ru]
      detail_meas[ind, :] = [ro_meas rv_meas ru_meas]
    end
    for j = 1:n_methods
        error_true[j] = mean(abs.(detail_true[:, j] - c_true))
        error_meas[j] = mean(abs.(detail_meas[:, j] - c_true))
    end
    (error_true, error_meas, detail_true, detail_meas)
end

# Calculate All Estimates
# Estimate the obs to converge and the estimation error using the Deltamax,
# Order, Variance, Uninformed methods for true, noisy departure times, and
# noisy departure order
# Input: inference parameters
# Output: obs to converge and estimation error for each method and noise level
function infer(param_inf::Paraminf)
    settings = param_inf.settings
    n_runs = param_inf.n_runs
    n_methods = param_inf.n_methods
    max_servers = param_inf.max_servers
    obs_max = param_inf.obs_max
    time_limit = param_inf.time_limit
    ϵ = param_inf.ϵ
    window = param_inf.window
    step = param_inf.step
    seed = param_inf.seed
    srand(seed)
    raw_true = Array{Array{Int64}}(n_runs, n_methods) # c estimates by obs
    raw_meas = Array{Array{Int64}}(n_runs, n_methods) # c estimates by obs
    conv_true = zeros(n_runs, n_methods) + obs_max - window * step # num obs to conv
    conv_meas = zeros(n_runs, n_methods) + obs_max - window * step # num obs to conv
    error_true = zeros(n_runs, n_methods) # c estimate with true departure times
    error_meas = zeros(n_runs, n_methods) # c estimate with noisy departure times
    # caculate each metric for each DOE run
    for i = 1:n_runs
        # initialize constants and data structures
        print("Run $i: ") # status update
        aname = settings[i, 2] # arrival distribution
        sname = settings[i, 3] # service distribution
        obs_limit = settings[i, 4] # number of observations
        c_true = settings[i, 5] # actual number of servers
        rho = settings[i, 6] # traffic intensity
        (adist, sdist) = builddist(c_true, aname, sname, rho)
        q = Queue(c_true, adist, sdist)
        s = Sim(obs_max, time_limit, rand(1:1000))
        q_out = runsim(s,q)
        output = disorder(q_out, c_true) # queue sim results with noise
        array_size = floor(Int64, (obs_max - 20) / step)
        ests_true = zeros(Int64, array_size, n_methods) # c estimates by obs by run
        ests_meas = zeros(Int64, array_size, n_methods) # c estimates by obs by run
        ests_err = zeros(Int64, array_size, n_methods) # c estimates by obs by run
        print("conv, ")
        (conv_true_single, conv_meas_single, ests_single, ests_meas_single) = convergence(param_inf, output, c_true)
        println("err")
        (error_true_single, error_meas_single, detail_single, detail_meas_single) = esterror(param_inf, output, obs_limit, c_true)
        conv_true[i,:] = conv_true_single
        conv_meas[i,:] = conv_meas_single
        error_true[i,:] = error_true_single
        error_meas[i,:] = error_meas_single
        for j = 1:n_methods
            raw_true[i,j] = ests_single[:,j]
            raw_meas[i,j] = ests_meas_single[:,j]
        end
  end
  resp_error = (error_true, error_meas)
  resp_conv = (conv_true, conv_meas)
  resp_raw = (raw_true, raw_meas)
  (resp_error, resp_conv, resp_raw)
end
