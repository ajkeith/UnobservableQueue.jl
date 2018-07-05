using SimJulia, ResumableFunctions
using Distributions, DataFrames, StatsBase
# using BenchmarkTools

mutable struct queue
    c::Int64 # number of servers
    adist::Distribution # arrival distribution
    sdist::Distribution # service distribution (single server)
end
queue() = queue(2, 1/2, 0.9, Exponential(1/λ), Exponential(1/μ))

mutable struct sim
    ncust::Int64 # total number of customers generated
    timelimit::Int64 # time limit for simulation
    seed::Int64 # rng seed
end
sim() = sim(1_000, 10_000, 8710)

# convert infparams and runparams to actuable structs
mutable struct infparams
    settings = readcsv(fn1, header = true)[1] # DOE runs
    # n_runs = length(settings[:,1])
    n_runs = 10 # limit runs for debugging and timing
    n_methods = 6 # number of methods to compare
    MAX_SERVERS = 19 # max number of servers
    obs_max = 1000 # max observations available
    ϵ = 0.05 # convergence quality
    window = 10 # observation window for convergence and estimates
    step = 20 # how many observations to skip while calculating convergence
end

mutable struct runparams
    obs_limit = settings[i, 4] # number of observations
    c_true = settings[i, 5] # actual number of servers
    fn2 = "QueueOutv18\\out" * "$i" * ".csv"
    q_out = readcsv(fn2, header = true)[1] # raw queue simulation results
    output = disorder(q_out, settings[i,:]) # queue sim results with noise
    lim = 20 # need at least 20 obs (since MAX_SERVERS = 19)
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

function runsim(s::sim, q::queue)
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
end

# Disorder Function
# Create noisy data for robust methods
# Input: queue sim output
# Output: true and noisy arrival and departure data
function disorder(r_out::Array{Float64,2}, num_servers::Int64)
  # actual data from simulation output
  A = r_out[:, 1] # arrive_order
  D_num = r_out[:, 3] # depart_order
  D_sort = r_out[:, 4] # depart_time
  D_raw = sort(D_sort) # sorted depart_time
  # D_raw_num = sortrows(r_out[:,3:4], by = x -> x[1])

  # constants and data structures
  NUM_SERVERS = num_servers
  MAX_CUSTS = length(r_out[:, 1])
  results = DataFrame()

  # Introduce partial observation error on cust depart order and time
  D_num_err = zeros(MAX_CUSTS)
  D_num_err[1:MAX_CUSTS] = D_num[1:MAX_CUSTS] # order error
  D_err = zeros(MAX_CUSTS) # misordered depart times
  D_err[1:MAX_CUSTS] = D_sort[1:MAX_CUSTS] # misordered depart times
  for i = 1:(MAX_CUSTS-NUM_SERVERS - 1)
    if rand() < 0.1
      swap = i + rand(1:NUM_SERVERS)
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
  for ind in convert.(Int64, r_out[:,3])
      incr += 1
      D_meas[incr] = D_raw_meas[ind]
  end

  DataFrame(x = 1:MAX_CUSTS, Arrive = A[1:MAX_CUSTS],
    Depart = D_sort[1:MAX_CUSTS], DepartMeas = D_meas[1:MAX_CUSTS],
    DepartErr = D_err[1:MAX_CUSTS], DepartOrder = D_num[1:MAX_CUSTS],
    DepartOrderErr = D_num_err[1:MAX_CUSTS])
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
function c_var_unf(A::Vector{Float64}, D::Vector{Float64}, cmax::Int64)
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
function c_var_unf_slow(A::Vector{Float64}, D::Vector{Float64}, cmax::Int64)
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

# separate estimation error and convergence functions? (yes)

function convergence()
    ind = 0
    lim = 20 # need at least 20 obs (since MAX_SERVERS = 19)
    bool_true = falses(n_methods)
    bool_meas = falses(n_methods)
    bool_err = falses(n_methods)
    while lim <= (obs_max - step)
      # println(lim)
      lim = lim + step
      ind += 1
      # true departure
      rd = c_est_order(Array(output[:DepartOrder])[1:lim])
      (rv, re, rn, rc) = c_est(Array(output[:Arrive])[1:lim], Array(output[:Depart])[1:lim], MAX_SERVERS)
      # rr = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:Depart])[1:lim], MAX_SERVERS)
      rr = 0
      # time error
      rd_meas = copy(rd)
      (rv_meas, re_meas, rn_meas, rc_meas) = c_est(Array(output[:Arrive])[1:lim], Array(output[:DepartMeas])[1:lim], MAX_SERVERS)
      # rr_meas = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:DepartMeas])[1:lim], MAX_SERVERS)
      rr_meas = 0
      # order error
      rd_err = c_est_order(Array(output[:DepartOrderErr])[1:lim])
      (rv_err, re_err, rn_err, rc_err) = c_est(Array(output[:Arrive])[1:lim], Array(output[:DepartErr])[1:lim], MAX_SERVERS)
      # rr_err = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:DepartErr])[1:lim], MAX_SERVERS)
      rr_err = 0
      # record results
      ests_true[ind, :] = [rd rv re rn rc rr]
      ests_meas[ind, :] = [rd_meas rv_meas re_meas rn_meas rc_meas rr_meas]
      ests_err[ind, :] = [rd_err rv_err re_err rn_err rc_err rr_err]
      if ind >= window
        for j = 1:n_methods
          if !bool_true[j] && mean(abs.(ests_true[(ind - window + 1):(ind), j] - c_true)) < ϵ
            conv_true[i, j] = (ind - window + 1) * step
            bool_true[j] = true
          end
          if !bool_meas[j] && mean(abs.(ests_meas[(ind - window + 1):(ind), j]- c_true)) < ϵ
            conv_meas[i, j] = (ind - window + 1) * step
            bool_meas[j] = true
          end
          if !bool_err[j] && mean(abs.(ests_err[(ind - window + 1):(ind), j]- c_true)) < ϵ
            conv_err[i, j] = (ind - window + 1) * step
            bool_err[j] = true
          end
        end
      end
    end
end

# Calculate All Estimates
# Estimate the obs to converge and the estimation error using the Deltamax,
# Order, Variance, Uninformed methods for true, noisy departure times, and
# noisy departure order
# Input: NA
# Output: obs to converge and estimation error for each method and noise level
function infer()
  # Initialize constants and data structures
  fn1 = "QueueOutv18\\design_v18.csv"
  settings = readcsv(fn1, header = true)[1] # DOE runs
  # n_runs = length(settings[:,1])
  n_runs = 10 # limit runs for debugging and timing
  n_methods = 6 # number of methods to compare
  MAX_SERVERS = 19 # max number of servers
  obs_max = 1000 # max observations available
  ϵ = 0.05 # convergence quality
  window = 10 # observation window for convergence and estimates
  step = 20 # how many observations to skip while calculating convergence
  error_true = zeros(n_runs, n_methods) # c estimate with true departure times
  error_meas = zeros(n_runs, n_methods) # c estimate with noisy departure times
  error_err = zeros(n_runs, n_methods) # c estimate with noisy departure order
  conv_true = zeros(n_runs, n_methods) + obs_max - window * step # num obs to conv
  conv_meas = zeros(n_runs, n_methods) + obs_max - window * step # num obs to conv
  conv_err = zeros(n_runs, n_methods) + obs_max - window * step # num obs to conv
  raw_true = Array{Array{Int64}}(n_runs, n_methods) # c estimates by obs
  raw_meas = Array{Array{Int64}}(n_runs, n_methods) # c estimates by obs
  raw_err = Array{Array{Int64}}(n_runs, n_methods) # c estimates by obs

  # caculate each metric for each DOE run
  for i = 1:n_runs
    # initialize constants and data structures
    print("Run $i: ") # status update
    obs_limit = settings[i, 4] # number of observations
    c_true = settings[i, 5] # actual number of servers
    fn2 = "QueueOutv18\\out" * "$i" * ".csv"
    q_out = readcsv(fn2, header = true)[1] # raw queue simulation results
    output = disorder(q_out, settings[i,:]) # queue sim results with noise
    # obs_max = size(output[:,1], 1) # deprecated
    array_size = floor(Int64, (obs_max - 20) / step)
    ests_true = zeros(Int64, array_size, n_methods) # c estimates by obs by run
    ests_meas = zeros(Int64, array_size, n_methods) # c estimates by obs by run
    ests_err = zeros(Int64, array_size, n_methods) # c estimates by obs by run

    # calculate number of observations until convergence
    print("conv, ")
    ind = 0
    lim = 20 # need at least 20 obs (since MAX_SERVERS = 19)
    bool_true = falses(n_methods)
    bool_meas = falses(n_methods)
    bool_err = falses(n_methods)
    while lim <= (obs_max - step)
      # println(lim)
      lim = lim + step
      ind += 1
      # true departure
      rd = c_est_order(Array(output[:DepartOrder])[1:lim])
      (rv, re, rn, rc) = c_est(Array(output[:Arrive])[1:lim], Array(output[:Depart])[1:lim], MAX_SERVERS)
      # rr = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:Depart])[1:lim], MAX_SERVERS)
      rr = 0
      # time error
      rd_meas = copy(rd)
      (rv_meas, re_meas, rn_meas, rc_meas) = c_est(Array(output[:Arrive])[1:lim], Array(output[:DepartMeas])[1:lim], MAX_SERVERS)
      # rr_meas = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:DepartMeas])[1:lim], MAX_SERVERS)
      rr_meas = 0
      # order error
      rd_err = c_est_order(Array(output[:DepartOrderErr])[1:lim])
      (rv_err, re_err, rn_err, rc_err) = c_est(Array(output[:Arrive])[1:lim], Array(output[:DepartErr])[1:lim], MAX_SERVERS)
      # rr_err = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:DepartErr])[1:lim], MAX_SERVERS)
      rr_err = 0
      # record results
      ests_true[ind, :] = [rd rv re rn rc rr]
      ests_meas[ind, :] = [rd_meas rv_meas re_meas rn_meas rc_meas rr_meas]
      ests_err[ind, :] = [rd_err rv_err re_err rn_err rc_err rr_err]
      if ind >= window
        for j = 1:n_methods
          if !bool_true[j] && mean(abs.(ests_true[(ind - window + 1):(ind), j] - c_true)) < ϵ
            conv_true[i, j] = (ind - window + 1) * step
            bool_true[j] = true
          end
          if !bool_meas[j] && mean(abs.(ests_meas[(ind - window + 1):(ind), j]- c_true)) < ϵ
            conv_meas[i, j] = (ind - window + 1) * step
            bool_meas[j] = true
          end
          if !bool_err[j] && mean(abs.(ests_err[(ind - window + 1):(ind), j]- c_true)) < ϵ
            conv_err[i, j] = (ind - window + 1) * step
            bool_err[j] = true
          end
        end
      end
    end

    # calculate estimation error
    println("est")
    window_detail = 50
    ind_l = Int64(obs_limit)
    ind_u = Int64(obs_limit  + window_detail - 1)
    detail_true = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
    detail_meas = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
    detail_err = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
    ind = 0
    for lim = ind_l:ind_u
      ind += 1
      # true departure
      rd = c_est_order(Array(output[:DepartOrder])[1:lim])
      (rv, re, rn, rc) = c_est(Array(output[:Arrive])[1:lim], Array(output[:Depart])[1:lim], MAX_SERVERS)
      # rr = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:Depart])[1:lim], MAX_SERVERS)
      rr = 0
      # time error
      rd_meas = copy(rd)
      (rv_meas, re_meas, rn_meas, rc_meas) = c_est(Array(output[:Arrive])[1:lim], Array(output[:DepartMeas])[1:lim], MAX_SERVERS)
      # rr_meas = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:DepartMeas])[1:lim], MAX_SERVERS)
      rr_meas = 0
      # order error
      rd_err = c_est_order(Array(output[:DepartOrderErr])[1:lim])
      (rv_err, re_err, rn_err, rc_err) = c_est(Array(output[:Arrive])[1:lim], Array(output[:DepartErr])[1:lim], MAX_SERVERS)
      # rr_err = c_est_robust(Array(output[:Arrive])[1:lim], Array(output[:DepartErr])[1:lim], MAX_SERVERS)
      rr_err = 0
      # record results
      detail_true[ind, :] = [rd rv re rn rc rr]
      detail_meas[ind, :] = [rd_meas rv_meas re_meas rn_meas rc_meas rr_meas]
      detail_err[ind, :] = [rd_err rv_err re_err rn_err rc_err rr_err]
    end
    for j = 1:n_methods
      error_true[i, j] = mean(abs.(detail_true[:, j] - c_true))
      error_meas[i, j] = mean(abs.(detail_meas[:, j] - c_true))
      error_err[i, j] = mean(abs.(detail_err[:, j] - c_true))
      raw_true[i, j] = ests_true[:, j]
      raw_meas[i, j] = ests_meas[:, j]
      raw_err[i, j] = ests_err[:, j]
    end
  end

  resp_error = (error_true, error_meas, error_err)
  resp_conv = (conv_true, conv_meas, conv_err)
  resp_raw = (raw_true, raw_meas, raw_err)
  (resp_error, resp_conv, resp_raw)
end
