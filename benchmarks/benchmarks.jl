include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
using BenchmarkTools

c = 7 # number of servers
λ = 0.99 # arrival rate (mean interarrival time is 1/λ)
μ = 1 / c # service rate (mean service time is 1/μ) (single server)
adist = Exponential(1/λ) # arrival distribution
sdist = Exponential(1/μ) # service distribution (single server)
ncust = 10000 # total number of customers generated
timelimit = 10_000 # time limit for simulation
seed = 1302 # rng seed
cmax = 19
s2 = Sim(ncust, timelimit, seed)
q2 = Queue(c, adist, sdist)
df1 = runsim(s2, q2)
A1, D1 = df1[:atime], df1[:dtime]
outorder1 = sort(df1, :dorder)[:aorder]

function f(df1::DataFrame)
    outorder1 = sort(df1, :dorder)[:aorder]
    c_order(outorder1)
end

bench_var = @benchmark c_var_unf($A1, $D1, $cmax)
bench_order = @benchmark f($df1)
