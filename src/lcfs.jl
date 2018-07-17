@resumable function customerLCFS(env::Environment, server::Resource, id::Integer, t_a::Float64, d_s::Distribution, starttime::Vector{Float64}, departtime::Vector{Float64})
    @yield timeout(env, t_a) # customer arrives
    # println("Customer $id arrived: ", now(env))
    @yield request(server, priority = -id) # customer starts service
    # println("Customer $id entered service: ", now(env))
    starttime[id] = now(env)
    @yield timeout(env, rand(d_s)) # server is busy
    @yield release(server) # customer exits service
    # println("Customer $id exited service: ", now(env))
    departtime[id] = now(env)
end

function runsimLCFS(s::Sim, q::Queue)
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
        @process customerLCFS(queuesim, server, i, arrivetime[i], sdist, starttime, departtime)
    end
    run(queuesim) # run simulation
    depart = hcat(sortrows(hcat(departtime, arriveorder)), collect(1:ncust))
    departorder = convert(Vector{Int64}, sortrows(depart, by = x -> (x[2],x[1],x[3]))[:,3])
    df = DataFrame(aorder = arriveorder, dorder = departorder, atime = arrivetime, stime = starttime, dtime = departtime)
    arrivals = DataFrame(aorder = df[:aorder], time = df[:atime])
    departs = DataFrame(dorder = df[:aorder], time = df[:dtime])
    dfl = join(arrivals, departs, on = :time, kind = :outer)
    sort!(dfl, :time)
    lcfsorder = zeros(Int, ncust)
    for i in 1:ncust
        a = dfl[:aorder][i]
        lcfsorder[i] = ismissing(a) ? dfl[:dorder][i] : a
    end
    lcfsorder
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
