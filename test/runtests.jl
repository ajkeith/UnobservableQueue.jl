include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
using Base.Test
using Distributions

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
        s1 = sim(ncust, timelimit, seed)
        q1 = queue(c, adist, sdist)
        df1 = runsim(s1, q1)
        @test mean(df1[:dtime] - df1[:stime]) ≈ (1 / μ) atol = 0.1
        @test df1[:atime] |> diff |> mean ≈ (1 / λ) atol = 0.1
        @test all(diff(df1[:atime]) .> 0)
        @test all(diff(sort(df1, :dorder)[:dtime]) .> 0)

        # G/G/15 queue
        c = 15 # number of servers
        μa = (1 / c) * 0.99
        μs = 1 / c
        σ = 0.9
        adist = LogNormal(exp(1/μ), σ) # arrival distribution
        sdist = LogNormal(exp(1/μ), σ) # service distribution (single server)
        ncust = 10_000 # total number of customers generated
        timelimit = 10_000 # time limit for simulation
        seed = 3001 # rng seed
        s2 = sim(ncust, timelimit, seed)
        q2 = queue(c, adist, sdist)
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
        s1 = sim(ncust, timelimit, seed)
        q1 = queue(c, adist, sdist)
        df1 = runsim(s1, q1)
        simout = hcat(df1[:atime], df1[:stime], df1[:dorder], df1[:dtime], df1[:dtime] - df1[:stime], df1[:stime] - df1[:atime], df1[:aorder])
        noiseout = disorder(simout, 2)
        @test count(noiseout[:DepartOrderErr] - noiseout[:DepartOrder] .== 0.0) / size(noiseout,1) ≈ 0.9 atol = 0.1
        @test all(sort(noiseout, :DepartMeas)[:DepartOrder] .== collect(1:size(noiseout,1)))

        # Deltamax algorithm
        o1 = [1,2,3,5,4]
        o2 = [1,3,6,5,4]
        @test c_deltamax(o1) == 2 && c_deltamax(o2) == 3

        # Variance algorithm

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
        s2 = sim(ncust, timelimit, seed)
        q2 = queue(c, adist, sdist)
        df2 = runsim(s2, q2)
        outorder1 = sort(df1, :dorder)[:aorder]
        outorder2 = sort(df2, :dorder)[:aorder]
        outorder3 = randperm(1_000)
        @test c_order(outorder3) == c_order_slow(outorder3)
        @test c_order([1,3,5]) == 3
        @test c_order(outorder1) == 2
        
        # Comparison
    end
end
