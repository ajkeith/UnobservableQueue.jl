@resumable function customerLCFS(env::Environment, server::Resource, id::Integer, t_a::Float64, d_s::Distribution, starttime::Vector{Float64}, departtime::Vector{Float64}, NA::Vector{Int}, ND::Vector{Int}, NC::Vector{Int})
    @yield timeout(env, t_a) # customer arrives
    NA[id] = NC[1]
    NC[1] += 1
    @yield request(server, priority = -id) # customer starts service
    starttime[id] = now(env)
    @yield timeout(env, rand(d_s)) # server is busy
    @yield release(server) # customer exits service
    departtime[id] = now(env)
    NC[1] -= 1
    ND[id] = NC[1]
end

function runsimLCFS(s::Sim, q::Queue)
    ncust = s.ncust
    seed = s.seed
    c = q.c
    adist = q.adist
    sdist = q.sdist
    srand(seed) # set rng seed
    # output of queue simulation
    # arrive order, arrive time, start service time, depart order, depart time
    arriveorder = collect(1:ncust)
    arrivetime = cumsum(rand(adist, ncust))
    starttime = zeros(ncust)
    departtime = zeros(ncust)
    NA = zeros(Int, ncust)
    ND = zeros(Int, ncust)
    NC = zeros(Int, 1)
    # setup and run simulation
    queuesim = Simulation() # initialize simulation environment
    server = Resource(queuesim, c) # initialize servers
    for i = 1:ncust # initialize customers
        @process customerLCFS(queuesim, server, i, arrivetime[i], sdist, starttime, departtime, NA, ND, NC)
    end
    run(queuesim) # run simulation
    depart = hcat(sortrows(hcat(departtime, arriveorder)), collect(1:ncust))
    departorder = convert(Vector{Int64}, sortrows(depart, by = x -> (x[2],x[1],x[3]))[:,3])
    DataFrame(aorder = arriveorder, dorder = departorder, atime = arrivetime, stime = starttime, dtime = departtime, na = NA, nd = ND)
end

function lcfs(df::DataFrame, dtype::Symbol)
    arrivals = DataFrame(aorder = df[:aorder], time = df[:atime])
    departs = DataFrame(daorder = df[:aorder], time = df[dtype])
    dfl = join(arrivals, departs, on = :time, kind = :outer)
    sort!(dfl, :time)
    len = size(dfl,1)
    lcfsorder = zeros(Int, len)
    for i in 1:len
        a = dfl[:aorder][i]
        lcfsorder[i] = ismissing(a) ? dfl[:daorder][i] : a
    end
    lcfsorder, dfl
end

# LCFS Disorder Function
# Create noisy data for robust methods
# Input: queue sim output
# Output: true and noisy arrival and departure data
function disorderLCFS(simout::DataFrame)
  len = size(simout, 1)
  D_meas = zeros(len) # measure error in depart times (original order)
  lcfsorder, dfl = lcfs(simout, :dtime)
  len2 = size(dfl,1)
  ismissing(dfl[:daorder][len2]) ? nothing : D_meas[dfl[:daorder][len2]] = dfl[:time][len2]
  for i = 2:(len2 - 1)
      ind = dfl[:daorder][i]
      if !ismissing(ind)
          a = mean(dfl[:time][(i-1):i])
          b = mean(dfl[:time][i:(i+1)])
          D_meas[ind] = a + rand() * (b - a)
      end
  end
  simout[:dmeas] = D_meas
  simout
end

function firstdepart(df::DataFrame, n::Int, t::Float64)
    dfdepart = sort(df, :dtime)
    indt = findfirst(dfdepart[:dtime] .> t, true)
    indt == 0 ? indt = size(dfdepart[:dtime],1) : nothing
    indn = findfirst(dfdepart[:nd][indt:end], n)
    indn == 0 ? indn = size(dfdepart[:nd][indt:end],1) : nothing
    dfdepart[:dtime][indt + indn - 1]
end

function c_order_LCFS(lcfsorder::Vector{Int64})
    nservers = 0
    N = size(lcfsorder, 1)
    est = zeros(Int, N)
    system = Set{Int}()
    service = Set{Int}()
    for i = 1:N
        if lcfsorder[i] in system
            setdiff!(system, lcfsorder[i])
            lcfsorder[i] in service && setdiff!(service, lcfsorder[i])
            if !isempty(system)
                notinservice = setdiff(system, service)
                !isempty(notinservice) && push!(service, maximum(notinservice))
            end
        else
            push!(system, lcfsorder[i])
        end
        est[i] = length(service)
    end
    maximum(est)
end

# Summary Estimator
# Estimate the number of servers using the Variance and Uninformed methods
# Input: arrival times, departure times, and the max number of servers
# Output: estimated number of servers for each method
function c_var_unf_LCFS(df::DataFrame, A::Vector{Float64}, D::Vector{Float64}, cmax::Int64)
  # Initilaize data structures and functions
  N = size(D, 1)
  B̂ = zeros(N)
  temp = zeros(N)
  VS = zeros(cmax, 2) # var, unf
  max_undef = 0
  NA = df[:na]
  ND = df[:nd]
  ## Grid Search
  for c = 1:cmax
    for i = 1:N
        if NA[i] < c
            B̂[i] = A[i]
        else
            B̂[i] = firstdepart(df, NA[i], A[i])
        end
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

# Convergence Estimator
# Calculate aprpoximate convergence
# Input: parametrs, queueing output data, number of servers
# Output: number of observations to approximately converge
function convergenceLCFS(param_inf::Paraminf, indrun::Int)
    settings = param_inf.settings
    aname = settings[indrun, 2] # arrival distribution
    sname = settings[indrun, 3] # service distribution
    c_true = settings[indrun, 5] # actual number of servers
    rho = settings[indrun, 6] # traffic intensity
    max_servers = param_inf.max_servers
    lim = param_inf.max_servers + 1
    n_methods = param_inf.n_methods
    obs_max = param_inf.obs_max
    window = param_inf.window
    step = param_inf.step
    ϵ = param_inf.ϵ
    bool_true = falses(n_methods)
    bool_meas = falses(n_methods)
    conv_true = zeros(n_methods)
    conv_meas = zeros(n_methods)
    array_size = floor(Int64, (obs_max - step) / step)
    ests_true = zeros(Int64, array_size, n_methods) # c estimates by obs by run
    ests_meas = zeros(Int64, array_size, n_methods) # c estimates by obs by run
    (adist, sdist) = builddist(c_true, aname, sname, rho)
    q = Queue(c_true, adist, sdist)
    seed = param_inf.seed
    for (ind, obs) in enumerate(lim:step:(obs_max - step))
        s = Sim(obs, time_limit, seed)
        q_out = runsimLCFS(s,q)
        output = disorderLCFS(q_out) # queue sim results with noise
        # true departure
        lcfsorder, dfl = lcfs(output, :dtime)
        rol = c_order_LCFS(lcfsorder)
        A1, D1 = output[:atime], output[:dtime]
        (rvl, rul, VSl) = c_var_unf_LCFS(output, A1, D1, max_servers)
        # time error
        lcfsorder, dfl = lcfs(output, :dmeas)
        rol_meas = c_order_LCFS(lcfsorder)
        A1, D1 = output[:atime], output[:dmeas]
        (rvl_meas, rul_meas, VSl_meas) = c_var_unf_LCFS(output, A1, D1, max_servers)
        # record results
        ests_true[ind, :] = [rol rvl rul]
        ests_meas[ind, :] = [rol_meas rvl_meas rul_meas]
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

function esterrorLCFS(param_inf::Paraminf, indrun::Int)
    settings = param_inf.settings
    aname = settings[indrun, 2] # arrival distribution
    sname = settings[indrun, 3] # service distribution
    obs_limit = settings[indrun, 4] # number of runs
    c_true = settings[indrun, 5] # actual number of servers
    rho = settings[indrun, 6] # traffic intensity
    window_detail = param_inf.window_detail
    n_methods = param_inf.n_methods
    max_servers = param_inf.max_servers
    ind_l = Int64(obs_limit)
    ind_u = Int64(obs_limit  + window_detail - 1)
    detail_true = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
    detail_meas = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
    error_true = zeros(n_methods)
    error_meas = zeros(n_methods)
    (adist, sdist) = builddist(c_true, aname, sname, rho)
    q = Queue(c_true, adist, sdist)
    seed = param_inf.seed
    for (ind, obs) in enumerate(ind_l:ind_u)
        s = Sim(obs, time_limit, seed)
        q_out = runsimLCFS(s,q)
        output = disorderLCFS(q_out) # queue sim results with noise
        # true departure
        lcfsorder, dfl = lcfs(output, :dtime)
        rol = c_order_LCFS(lcfsorder)
        A1, D1 = output[:atime], output[:dtime]
        (rvl, rul, VSl) = c_var_unf_LCFS(output, A1, D1, max_servers)
        # time error
        lcfsorder, dfl = lcfs(output, :dmeas)
        rol_meas = c_order_LCFS(lcfsorder)
        A1, D1 = output[:atime], output[:dmeas]
        (rvl_meas, rul_meas, VSl_meas) = c_var_unf_LCFS(output, A1, D1, max_servers)
        # record results
        detail_true[ind, :] = [rol rvl rul]
        detail_meas[ind, :] = [rol_meas rvl_meas rul_meas]
    end
    for j = 1:n_methods
        error_true[j] = mean(abs.(detail_true[:, j] - c_true))
        error_meas[j] = mean(abs.(detail_meas[:, j] - c_true))
    end
    (error_true, error_meas, detail_true, detail_meas)
end

function inferLCFS(param_inf::Paraminf)
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
        print("Run $i: ") # status update
        array_size = floor(Int64, (obs_max - 20) / step)
        ests_true = zeros(Int64, array_size, n_methods) # c estimates by obs by run
        ests_meas = zeros(Int64, array_size, n_methods) # c estimates by obs by run
        ests_err = zeros(Int64, array_size, n_methods) # c estimates by obs by run
        print("conv, ")
        (conv_true[i,:], conv_meas[i,:], ests_single, ests_meas_single) = convergenceLCFS(param_inf, i)
        println("err")
        (error_true[i,:], error_meas[i,:], detail_single, detail_meas_single) = esterrorLCFS(param_inf, i)
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
