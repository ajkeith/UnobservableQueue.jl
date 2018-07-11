include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
# using Plots; gr()

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
seed = 8710 # random number seed
n_runs = size(settings,1) # trial runs
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
(rerr, rconv, rraw) = infer(param_inf)

using CSV, DataFrames
err = DataFrame(rerr[1])
CSV.write("data\\err.csv", err)
err_meas = DataFrame(rerr[2])
CSV.write("data\\err_meas.csv", err_meas)
conv = DataFrame(rconv[1])
CSV.write("data\\conv.csv", conv)
conv_meas = DataFrame(rconv[2])
CSV.write("data\\conv_meas.csv", conv_meas)

using JLD
save("data\\output.jld","err", err, "err_meas", err_meas, "conv", conv, "conv_meas", conv_meas, "raw", rraw[1], "raw_meas", rraw[2])
data = load("data\\output.jld")

using Distributions, StatPlots
h = histogram(err[1])
histogram!(err[2])
