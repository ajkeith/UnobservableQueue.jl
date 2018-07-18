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
    df = DataFrame(aorder = arriveorder, dorder = departorder, atime = arrivetime, stime = starttime, dtime = departtime, na = NA, nd = ND)
    arrivals = DataFrame(aorder = df[:aorder], time = df[:atime])
    departs = DataFrame(dorder = df[:aorder], time = df[:dtime])
    dfl = join(arrivals, departs, on = :time, kind = :outer)
    sort!(dfl, :time)
    lcfsorder = zeros(Int, ncust)
    for i in 1:ncust
        a = dfl[:aorder][i]
        lcfsorder[i] = ismissing(a) ? dfl[:dorder][i] : a
    end
    df, lcfsorder
end

function firstdepart(df::DataFrame, n::Int, t::Float64)
    dfdepart = sort(df, :dtime)
    indt = findfirst(dfdepart[:dtime] .> t, true)
    indn = findfirst(dfdepart[:nd][indt:end], n)
    dfdepart[:dtime][indn]
end

function c_order_LCFS(lcfsorder::Vector{Int64}, ncust::Int)
    nservers = 0
    est = zeros(Int, ncust)
    N = size(lcfsorder, 1)
    system = Set{Int}()
    service = Set{Int}()
    for i = 1:N
        if lcfsorder[i] in system
            # @show lcfsorder[i]
            # @show system
            setdiff!(system, lcfsorder[i])
            lcfsorder[i] in service && setdiff!(service, lcfsorder[i])
            if !isempty(system)
                notinservice = setdiff(system, service)
                !isempty(notinservice) && push!(service, maximum(notinservice))
            end
        else
            push!(system, lcfsorder[i])
        end
        # @show service
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
  (r1, r2, VS, B̂)
end
