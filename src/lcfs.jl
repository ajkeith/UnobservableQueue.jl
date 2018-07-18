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
    indn = findfirst(dfdepart[:nd][indt:end], n)
    dfdepart[:dtime][indn]
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
    maximum(est), est
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
#
# # Convergence Estimator
# # Calculate aprpoximate convergence
# # Input: parametrs, queueing output data, number of servers
# # Output: number of observations to approximately converge
# function convergence(param_inf::Paraminf, queue_output::DataFrame, c_true::Int64)
#     ind = 0
#     max_servers = param_inf.max_servers
#     lim = param_inf.max_servers + 1
#     n_methods = param_inf.n_methods
#     obs_max = param_inf.obs_max
#     window = param_inf.window
#     step = param_inf.step
#     ϵ = param_inf.ϵ
#     output = queue_output
#     bool_true = falses(n_methods)
#     bool_meas = falses(n_methods)
#     conv_true = zeros(n_methods)
#     conv_meas = zeros(n_methods)
#     array_size = floor(Int64, (obs_max - 20) / step)
#     ests_true = zeros(Int64, array_size, n_methods) # c estimates by obs by run
#     ests_meas = zeros(Int64, array_size, n_methods) # c estimates by obs by run
#     while lim <= (obs_max - step)
#       # println(lim)
#       lim = lim + step
#       ind += 1
#       # true departure
#       out_order = sort(output[1:lim,:], :DepartOrder)[:ArriveOrder]
#       ro = c_order(out_order)
#       (rv, ru, VS) = c_var_unf_LCFS(output[:Arrive][1:lim], output[:Depart][1:lim], max_servers)
#       # time error
#       out_order_meas = sort(output[1:lim,:], :DepartMeas)[:ArriveOrder]
#       ro_meas = c_order(out_order_meas)
#       (rv_meas, ru_meas, VS_meas) = c_var_unf(output[:Arrive][1:lim], output[:DepartMeas][1:lim], max_servers)
#       # record results
#       ests_true[ind, :] = [ro rv ru]
#       ests_meas[ind, :] = [ro_meas rv_meas ru_meas]
#       if ind >= window
#         for j = 1:n_methods
#           if !bool_true[j] && mean(abs.(ests_true[(ind - window + 1):(ind), j] - c_true)) < ϵ
#             conv_true[j] = (ind - window + 1) * step
#             bool_true[j] = true
#           end
#           if !bool_meas[j] && mean(abs.(ests_meas[(ind - window + 1):(ind), j]- c_true)) < ϵ
#             conv_meas[j] = (ind - window + 1) * step
#             bool_meas[j] = true
#           end
#         end
#       end
#     end
#     (conv_true, conv_meas, ests_true, ests_meas)
# end
#
# # Estimation Error Estimator
# # Calculate estimation error
# # Input: parametrs, queue output data, obs limit, number of servers
# # Output: average estimation error in estimation window
# function esterror(param_inf::Paraminf, output::DataFrame, obs_limit::Int64, c_true::Int64)
#     window_detail = param_inf.window_detail
#     n_methods = param_inf.n_methods
#     max_servers = param_inf.max_servers
#     ind_l = Int64(obs_limit)
#     ind_u = Int64(obs_limit  + window_detail - 1)
#     detail_true = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
#     detail_meas = zeros(Int64, window_detail, n_methods) - 1 # precise c estimates
#     error_true = zeros(n_methods)
#     error_meas = zeros(n_methods)
#     ind = 0
#     for lim = ind_l:ind_u
#       ind += 1
#       # true departure
#       out_order = sort(output[1:lim,:], :DepartOrder)[:ArriveOrder]
#       ro = c_order(out_order)
#       (rv, ru, VS) = c_var_unf(output[:Arrive][1:lim], output[:Depart][1:lim], max_servers)
#       # time error
#       out_order_meas = sort(output[1:lim,:], :DepartMeas)[:ArriveOrder]
#       ro_meas = c_order(out_order_meas)
#       (rv_meas, ru_meas, VS_meas) = c_var_unf(output[:Arrive][1:lim], output[:DepartMeas][1:lim], max_servers)
#       # record results
#       detail_true[ind, :] = [ro rv ru]
#       detail_meas[ind, :] = [ro_meas rv_meas ru_meas]
#     end
#     for j = 1:n_methods
#         error_true[j] = mean(abs.(detail_true[:, j] - c_true))
#         error_meas[j] = mean(abs.(detail_meas[:, j] - c_true))
#     end
#     (error_true, error_meas, detail_true, detail_meas)
# end
#
# # Calculate All Estimates
# # Estimate the obs to converge and the estimation error using the Deltamax,
# # Order, Variance, Uninformed methods for true, noisy departure times, and
# # noisy departure order
# # Input: inference parameters
# # Output: obs to converge and estimation error for each method and noise level
# function infer(param_inf::Paraminf)
#     settings = param_inf.settings
#     n_runs = param_inf.n_runs
#     n_methods = param_inf.n_methods
#     max_servers = param_inf.max_servers
#     obs_max = param_inf.obs_max
#     time_limit = param_inf.time_limit
#     ϵ = param_inf.ϵ
#     window = param_inf.window
#     step = param_inf.step
#     seed = param_inf.seed
#     srand(seed)
#     raw_true = Array{Array{Int64}}(n_runs, n_methods) # c estimates by obs
#     raw_meas = Array{Array{Int64}}(n_runs, n_methods) # c estimates by obs
#     conv_true = zeros(n_runs, n_methods) + obs_max - window * step # num obs to conv
#     conv_meas = zeros(n_runs, n_methods) + obs_max - window * step # num obs to conv
#     error_true = zeros(n_runs, n_methods) # c estimate with true departure times
#     error_meas = zeros(n_runs, n_methods) # c estimate with noisy departure times
#     # caculate each metric for each DOE run
#     for i = 1:n_runs
#         # initialize constants and data structures
#         print("Run $i: ") # status update
#         aname = settings[i, 2] # arrival distribution
#         sname = settings[i, 3] # service distribution
#         obs_limit = settings[i, 4] # number of observations
#         c_true = settings[i, 5] # actual number of servers
#         rho = settings[i, 6] # traffic intensity
#         (adist, sdist) = builddist(c_true, aname, sname, rho)
#         q = Queue(c_true, adist, sdist)
#         s = Sim(obs_max, time_limit, rand(1:1000))
#         q_out, lcfsorder = runsimLCFS(s,q)
#         output = disorder(q_out, c_true) # queue sim results with noise
#         array_size = floor(Int64, (obs_max - 20) / step)
#         ests_true = zeros(Int64, array_size, n_methods) # c estimates by obs by run
#         ests_meas = zeros(Int64, array_size, n_methods) # c estimates by obs by run
#         ests_err = zeros(Int64, array_size, n_methods) # c estimates by obs by run
#         print("conv, ")
#         (conv_true_single, conv_meas_single, ests_single, ests_meas_single) = convergence(param_inf, output, c_true)
#         println("err")
#         (error_true_single, error_meas_single, detail_single, detail_meas_single) = esterror(param_inf, output, obs_limit, c_true)
#         conv_true[i,:] = conv_true_single
#         conv_meas[i,:] = conv_meas_single
#         error_true[i,:] = error_true_single
#         error_meas[i,:] = error_meas_single
#         for j = 1:n_methods
#             raw_true[i,j] = ests_single[:,j]
#             raw_meas[i,j] = ests_meas_single[:,j]
#         end
#   end
#   resp_error = (error_true, error_meas)
#   resp_conv = (conv_true, conv_meas)
#   resp_raw = (raw_true, raw_meas)
#   (resp_error, resp_conv, resp_raw)
# end
