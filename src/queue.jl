
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
    srand(seed) # set rng seed
    # output of queue simulation
    # arrive order, arrive time, start service time, depart order, depart time
    arriveorder = collect(1:ncust)
    arrivetime = cumsum(rand(adist, ncust))
    starttime = zeros(ncust)
    departtime = zeros(ncust)
    # setup and run simulation
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
  D_raw_meas[1] = D_raw[1] / 2 + rand() * (mean(D_raw[1:2]) - D_raw[1] / 2)
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
