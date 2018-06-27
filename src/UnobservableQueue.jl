using SimJulia, ResumableFunctions
using Distributions, DataFrames
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
# Input: departure order
# Output: estimated number of servers
function c_deltamax(out_order::Vector{Number})
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
