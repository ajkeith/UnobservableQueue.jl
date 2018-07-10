using SimJulia, ResumableFunctions
using Distributions, DataFrames, StatsBase
# using BenchmarkTools

mutable struct Queue
    c::Int64 # number of servers
    adist::Distribution # arrival distribution
    sdist::Distribution # service distribution (single server)
end
Queue() = Queue(2, Exponential(1/λ), Exponential(1/μ))

mutable struct Sim
    ncust::Int64 # total number of customers generated
    timelimit::Int64 # time limit for simulation
    seed::Int64 # rng seed
end
Sim() = Sim(1_000, 10_000, 8710)

# convert Paraminf and Paramrun to actual structs
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
end

# define customer behavior
@resumable function customer(env::Environment, server::Resource, id::Integer, t_a::Float64, d_s::Distribution, starttime::Vector{Float64}, departtime::Vector{Float64})
    @yield timeout(env, t_a) # customer arrives
    # println("Customer $id arrived: ", now(env))
    @yield request(server) # customer starts service
    # println("Customer $id entered service: ", now(env))
    starttime[id] = now(env)
    @yield timeout(env, rand(d_s)) # server is busy
    @yield release(server) # customer exits service
    # println("Customer $id exited service: ", now(env))
    departtime[id] = now(env)
end

function runsim(s::Sim, q::Queue)
    ncust = s.ncust
    seed = s.seed
    c = q.c
    adist = q.adist
    sdist = q.sdist

    # output of queue simulation
    # arrive order, arrive time, start service time, depart order, depart time
    arriveorder = collect(1:ncust)
    arrivetime = cumsum(rand(adist, ncust))
    starttime = zeros(ncust)
    departtime = zeros(ncust)

    # setup and run simulation
    srand(seed) # set rng seed
    queuesim = Simulation() # initialize simulation environment
    server = Resource(queuesim, c) # initialize servers
    for i = 1:ncust # initialize customers
        @process customer(queuesim, server, i, arrivetime[i], sdist, starttime, departtime)
    end
    run(queuesim) # run simulation
    depart = hcat(sortrows(hcat(departtime, arriveorder)), collect(1:ncust))
    departorder = convert(Vector{Int64}, sortrows(depart, by = x -> (x[2],x[1],x[3]))[:,3])
    DataFrame(aorder = arriveorder, atime = arrivetime, stime = starttime, dorder = departorder, dtime = departtime)
    # hcat(arrivetime, starttime, departorder, departtime, departtime - starttime, starttime - arrivetime, arriveorder)
end

# Disorder Function
# Create noisy data for robust methods
# Input: queue sim output
# Output: true and noisy arrival and departure data
function disorder(simout::DataFrame, num_servers::Int64)
  # actual data from simulation output
  A = simout[:atime] # arrive time
  D_num = simout[:dorder] # depart_order
  D_sort = simout[:dtime] # depart_time
  D_raw = sort(D_sort) # sorted depart_time

  # constants and data structures
  max_custs = length(A)
  results = DataFrame()

  # Introduce partial observation error on cust depart order and time
  D_num_err = zeros(max_custs)
  D_num_err[1:max_custs] = D_num[1:max_custs] # order error
  D_err = zeros(max_custs) # misordered depart times
  D_err[1:max_custs] = D_sort[1:max_custs] # misordered depart times
  for i = 1:(max_custs - num_servers - 1)
    if rand() < 0.1
      swap = i + rand(1:num_servers)
      temp1 = D_num_err[i]
      temp2 = D_err[i]
      D_num_err[i] = D_num_err[swap]
      D_num_err[swap] = temp1
      D_err[i] = D_err[swap]
      D_err[swap] = temp2
    end
  end
  len = size(D_raw, 1)
  D_raw_meas = zeros(len) # measure error in depart times (original order)
  D_raw_meas[1] = D_raw[1] / 2 + rand() * mean(D_raw[1:2]) - D_raw[1] / 2
  D_raw_meas[len] = D_raw[len]
  for i = 2:(len - 1)
    a = mean(D_raw[(i-1):i])
    b = mean(D_raw[i:(i+1)])
    D_raw_meas[i] = a + rand() * (b - a)
  end
  D_meas = zeros(len)
  incr = 0
  for ind in simout[:dorder]
      incr += 1
      D_meas[incr] = D_raw_meas[ind]
  end

  DataFrame(ArriveOrder = 1:max_custs, Arrive = A[1:max_custs],
    Depart = D_sort[1:max_custs], DepartMeas = D_meas[1:max_custs],
    DepartErr = D_err[1:max_custs], DepartOrder = D_num[1:max_custs],
    DepartOrderErr = D_num_err[1:max_custs])
end


# Deltamax Estimator
# Estimate the number of servers using the Deltamax algorithm
# Input: arrival order sorted by departure order
# Output: estimated number of servers
function c_deltamax(out_order::Vector{Int64})
  N = size(out_order, 1)
  currmax = out_order[1]
  deltamax = 0
  for i = 2:N
    nextmax = max(out_order[i], currmax)
    deltamax = max(deltamax, nextmax - currmax)
    currmax = nextmax
  end
  deltamax
end

# Order-based Estimator
# Estimate the number of servers using the Order-based algorithm
# Input: arrival order sorted by departure order
# Output: estimated number of servers
function c_order_slow(out_order::Vector{Int64})
    N = size(out_order, 1)
    currmax = 0
    nservers = 0
    for i = 1:N
        currmax = max(out_order[i], currmax)
        custall = Set(1:currmax)
        custdeparted = Set(out_order[1:i])
        custservice = setdiff(custall, custdeparted)
        nservers = max(length(custservice) + 1, nservers)
    end
    nservers
end

function c_order(out_order::Vector{Int64})
    N = size(out_order, 1)
    currmax = 0
    noccupied = 0
    nservers = 0
    for i = 1:N
        nextmax = max(out_order[i], currmax)
        if out_order[i] < currmax
            noccupied = noccupied - 1
        else
            noccupied = noccupied + nextmax - currmax - 1
        end
        currmax = nextmax
        nservers = max(nservers, noccupied + 1)
    end
    nservers
end

# Summary Estimator - fast version
# may be incorrect
function c_var_unf_fast(A::Vector{Float64}, D::Vector{Float64}, cmax::Int64)
  # Initilaize data structures and functions
  N = size(D, 1)
  B̂ = zeros(N)
  temp = zeros(N)
  VS = zeros(cmax, 2) # var, unf
  VS[1,:] = [-1,-1]
  max_undef = 1
  Dsort = sort(D)
  ## Grid Search
  for c = 2:cmax
    B̂[1:c] = A[1:c]
    for i = (c + 1):N
      B̂[i] = max(A[i], Dsort[i - c])
    end
    Ŝ = D - B̂
    if any(Ŝ .< 0)
      VS[c, 1:2] = -1
      max_undef = c
    else
      VS[c, 1] = var(Ŝ) * (N-1) # var
    end
  end
  (r1, r2) = (99, 99) # initialize result to 99 for debugging
  # var method
  try
    r1 = findin(VS[:, 1], minimum(filter(x -> x >= 0, VS[:, 1])))[1]
  catch
    r1 = 1
  end
  # uninformed method
  r2 = max_undef + 1
  (r1, r2, VS)
end

# Summary Estimator
# Estimate the number of servers using the Variance and Uninformed methods
# Input: arrival times, departure times, and the max number of servers
# Output: estimated number of servers for each method
# Note: D(k,m) the kth order statistic of the departure times of the first m
# customers is the same as the kth order statistic of the departure times in
# increasing time order
function c_var_unf(A::Vector{Float64}, D::Vector{Float64}, cmax::Int64)
  # Initilaize data structures and functions
  N = size(D, 1)
  B̂ = zeros(N)
  temp = zeros(N)
  VS = zeros(cmax, 2) # var, unf
  max_undef = 0
  ## Grid Search
  for c = 1:cmax
    B̂[1:c] = A[1:c]
    for i = (c + 1):N
      Dₖ = sort(D[1:(i - 1)])
      B̂[i] = max(A[i], Dₖ[i - c])
    end
    Ŝ = D - B̂
    if any(Ŝ .< 0)
      VS[c, 1:2] = -1
      max_undef = c
    else
      VS[c, 1] = var(Ŝ) * (N-1) # var
    end
  end
  (r1, r2) = (99, 99) # initialize result to 99 for debugging
  # var method
  try
    r1 = findin(VS[:, 1], minimum(filter(x -> x >= 0, VS[:, 1])))[1]
  catch
    r1 = 1
  end
  # uninformed method
  r2 = max_undef + 1
  (r1, r2, VS)
end

# separated distribution function
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

# separated convergence function
function convergence(param_inf::Paraminf, queue_output::DataFrame, c_true::Int64)
    ind = 0
    lim = param_inf.max_servers + 1
    n_methods = param_inf.n_methods
    obs_max = param_inf.obs_max
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

# calculate estimation error
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
# Input: NA
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
