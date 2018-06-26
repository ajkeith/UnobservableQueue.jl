include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
using Base.Test
using Distributions

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
    df = runsim(s1, q1)
    @test mean(df[:dtime] - df[:stime]) ≈ (1 / μ) atol = 0.1
    @test df[:atime] |> diff |> mean ≈ (1 / λ) atol = 0.1
    @test all(diff(df[:atime]) .> 0)
    @test all(diff(sort(df, :dorder)[:dtime]) .> 0)

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
    df = runsim(s2, q2)
    @test mean(df[:dtime] - df[:stime]) ≈ mean(rand(sdist,10_000)) atol = 100
    @test df[:atime] |> diff |> mean ≈ mean(rand(adist, 10_000)) atol = 100
    @test all(diff(df[:atime]) .> 0)
    @test all(diff(sort(df, :dorder)[:dtime]) .> 0)
end
