include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
using Base.Test
using Distributions
using Plots; gr()

@testset "Unobservable Queue" begin
    @testset "Queue Simulation" begin
        # M/M/2 queue
        c = 2 # number of servers
        μ = 1 / c # service rate (mean service time is 1/μ) (single server)
        λ = 0.9 # arrival rate (mean interarrival time is 1/λ)
        adist = Exponential(1/λ) # arrival distribution
        sdist = Exponential(1/μ) # service distribution (single server)
        ncust = 10_000 # total number of customers generated
        timelimit = 1_000 # time limit for simulation
        seed = 8710 # rng seed
        s1 = Sim(ncust, timelimit, seed)
        q1 = Queue(c, adist, sdist)
        df1 = runsim(s1, q1)
        @test mean(df1[:dtime] - df1[:stime]) ≈ (1 / μ) atol = 0.1
        @test df1[:atime] |> diff |> mean ≈ (1 / λ) atol = 0.1
        @test all(diff(df1[:atime]) .> 0)
        @test all(diff(sort(df1, :dorder)[:dtime]) .> 0)

        # G/G/15 queue
        c = 15 # number of servers
        μa = (1 / 2) * 0.99
        μs = 1 / 2
        σ = 0.9
        adist = LogNormal(exp(1/μa), σ) # arrival distribution
        sdist = LogNormal(exp(1/μs), σ) # service distribution (single server)
        ncust = 10_000 # total number of customers generated
        timelimit = 10_000 # time limit for simulation
        seed = 3001 # rng seed
        s2 = Sim(ncust, timelimit, seed)
        q2 = Queue(c, adist, sdist)
        df2 = runsim(s2, q2)
        @test mean(df2[:dtime] - df2[:stime]) ≈ mean(rand(sdist, 10_000)) atol = 100
        @test df2[:atime] |> diff |> mean ≈ mean(rand(adist, 10_000)) atol = 100
        @test all(diff(df2[:atime]) .> 0)
        @test all(diff(sort(df2, :dorder)[:dtime]) .> 0)
    end

    @testset "Server Estimation" begin
        # Introduce observation error into true data
        c = 2 # number of servers
        μ = 1 / c # service rate (mean service time is 1/μ) (single server)
        λ = 0.9 # arrival rate (mean interarrival time is 1/λ)
        adist = Exponential(1/λ) # arrival distribution
        sdist = Exponential(1/μ) # service distribution (single server)
        ncust = 10_000 # total number of customers generated
        timelimit = 1_000 # time limit for simulation
        seed = 8710 # rng seed
        s1 = Sim(ncust, timelimit, seed)
        q1 = Queue(c, adist, sdist)
        df1 = runsim(s1, q1)
        noiseout = disorder(df1, 2)
        @test count(noiseout[:DepartOrderErr] - noiseout[:DepartOrder] .== 0.0) / size(noiseout,1) ≈ 0.9 atol = 0.1
        @test all(sort(noiseout, :DepartMeas)[:DepartOrder] .== collect(1:size(noiseout,1)))
        @test all(sort(noiseout, :DepartMeas)[:DepartOrder] .== sort(noiseout, :Depart)[:DepartOrder])

        # Deltamax algorithm
        o1 = [1,2,3,5,4]
        o2 = [1,3,6,5,4]
        @test c_deltamax(o1) == 2 && c_deltamax(o2) == 3

        # Order-based algorithm
        c = 15 # number of servers
        μa = (1 / c) * 0.99
        μs = 1 / c
        σ = 0.9
        adist = LogNormal(exp(1/μ), σ) # arrival distribution
        sdist = LogNormal(exp(1/μ), σ) # service distribution (single server)
        ncust = 10_000 # total number of customers generated
        timelimit = 10_000 # time limit for simulation
        seed = 3001 # rng seed
        s2 = Sim(ncust, timelimit, seed)
        q2 = Queue(c, adist, sdist)
        df2 = runsim(s2, q2)
        outorder1 = sort(df1, :dorder)[:aorder]
        outorder2 = sort(df2, :dorder)[:aorder]
        outorder3 = randperm(1_000)
        @test c_order(outorder3) == c_order_slow(outorder3)
        @test c_order([1,3,5]) == 3
        @test c_order(outorder1) == 2

        # Variance algorithm
        cmax = 19
        A1, A2 = df1[:atime][1:1000], df2[:atime][1:1000]
        D1, D2 = df1[:dtime][1:1000], df2[:dtime][1:1000]
        (cvar1, cunf1, VS1) = c_var_unf(A1, D1, cmax)
        (cvar2, cunf2, VS2) = c_var_unf(A2, D2, cmax)
        @test cvar1 == 2 && cunf1 == 2

        # Build distributions
        (adist, sdist) = builddist(2, "Uniform", "LogNormal", 0.99)
        (adist2, sdist2) = builddist(15, "Exponential", "Beta", 0.9)
        @test (1/mean(adist)) / (2 * (1/mean(sdist))) ≈ 0.99 atol = 0.0001
        @test (1/mean(adist2)) / (15 * (1/mean(sdist2))) ≈ 0.9 atol = 0.0001

        # Convergence estimation
        settings = ["Order" "Exponential" "Exponential" 400 9 0.9;
                    "Order" "Uniform" "LogNormal" 400 15 0.9]
        n_runs = 10 # limit runs for debugging and timing
        n_methods = 3 # number of methods to compare
        max_servers = 19 # max number of servers
        obs_max = 1_000 # max observations available
        time_limit = 10_000 # max simulation time
        ϵ = 0.05 # convergence quality
        window = 20 # observation window for convergence estimate
        window_detail = 50 # observation window for error estimate
        step = 20 # how many observations to skip while calculating convergence
        seed = 8710
        param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
        c = settings[2,5] # number of servers
        aname = convert(String, settings[2,2])
        sname = convert(String, settings[2,3])
        rho = settings[2,6]
        (adist, sdist) = builddist(c, aname, sname, rho)
        ncust = 1_000 # total number of customers generated
        timelimit = 1_000 # time limit for simulation
        seed = 8710 # rng seed
        s1 = Sim(ncust, timelimit, seed)
        q1 = Queue(c, adist, sdist)
        df1 = runsim(s1, q1)
        output = disorder(df1, c)
        (conv, conv_meas, ests, ests_meas) = convergence(param_inf, output, c)
        @test (conv[1] > 0) && (conv[1] < 400) && (conv[2] == 0)
        @test (conv[1] == conv_meas[1]) && (conv_meas[2] == 0)

        # Error estimation
        settings = ["Order" "Exponential" "Exponential" 400 9 0.9;
                    "Order" "Uniform" "LogNormal" 400 15 0.9;
                    "Order" "LogNormal" "LogNormal" 60 9 0.9]
        n_runs = 10 # limit runs for debugging and timing
        n_methods = 3 # number of methods to compare
        max_servers = 19 # max number of servers
        obs_max = 1_000 # max observations available
        time_limit = 10_000 # max simulation time
        ϵ = 0.05 # convergence quality
        window = 20 # observation window for convergence estimate
        window_detail = 50 # observation window for error estimate
        step = 20 # how many observations to skip while calculating convergence
        seed = 8710
        param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
        ncust = 1_000 # total number of customers generated
        timelimit = 1_000 # time limit for simulation
        seed = 8710 # rng seed
        ind = 3
        c = settings[ind,5] # number of servers
        aname = convert(String, settings[ind,2])
        sname = convert(String, settings[ind,3])
        obs_limit = settings[ind,4]
        rho = settings[ind,6]
        (adist, sdist) = builddist(c, aname, sname, rho)
        s1 = Sim(ncust, timelimit, seed)
        q1 = Queue(c, adist, sdist)
        df1 = runsim(s1, q1)
        output = disorder(df1, c)
        (err, err_meas, detail, detail_meas) = esterror(param_inf, output, obs_limit, c)
        @test err[1] == 0 && err[2] >= 1
        @test err_meas[1] == 0 && err_meas[2] >= 8

        # Comparison
        settings = ["Order" "Exponential" "Exponential" 400 9 0.9;
                    "Order" "Uniform" "LogNormal" 400 15 0.9;
                    "Order" "LogNormal" "LogNormal" 60 9 0.9]
        n_runs = size(settings,1) # experimental runs
        n_methods = 3 # number of methods to compare
        max_servers = 19 # max number of servers
        obs_max = 1_000 # max observations available
        time_limit = 10_000 # max simulation time
        ϵ = 0.05 # convergence quality
        window = 20 # observation window for convergence estimate
        window_detail = 50 # observation window for error estimate
        step = 20 # how many observations to skip while calculating convergence
        seed = 8710
        param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
        (rerr, rconv, rraw) = infer(param_inf)
        @test mean(rerr[1][:,1]) < mean(rerr[1][:,2])
        @test mean(rerr[2][:,1]) < mean(rerr[1][:,2])
        @test (rconv[1][1,1] < rconv[1][1,2]) || (rconv[1][1,2] == 0)
        @test (rconv[2][1,1] < rconv[1][1,2]) || (rconv[2][1,2] == 0)
    end
end
