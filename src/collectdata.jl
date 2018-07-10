include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
using Plots; gr()

# Small trial run
cd("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl")
fn = "data\\design.csv"
input = readcsv(fn, header = true)[1]
settings = Array{Any,2}(size(input))
settings[:, 1:3] = convert.(String, input[:, 1:3])
settings[:, 4:6] = input[:,4:6]
n_methods = 3 # number of methods to compare
max_servers = 19 # max number of servers
obs_max = 1_000 # max observations available
time_limit = 10_000 # max simulation time
ϵ = 0.05 # convergence quality
window = 20 # observation window for convergence estimate
window_detail = 50 # observation window for error estimate
step = 20 # how many observations to skip while calculating convergence
n_runs = 10 # trial runs
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step)
(rerr, rconv, rraw) = infer(param_inf)
