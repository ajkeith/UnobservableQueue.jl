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
