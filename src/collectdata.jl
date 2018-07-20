using CSV, DataFrames, JLD
include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
# using Plots; gr()

# Initial settings
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
window = 10 # observation window for convergence estimate
window_detail = 50 # observation window for error estimate
step = 20 # how many observations to skip while calculating convergence
trialruns = 10 # how many runs to check data transfer

# FCFS Data Collection Part 1
seed = 8710
n_runs = trialruns
# n_runs = size(settings,1)
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
@time (rerr, rconv, rraw) = infer(param_inf)
err = DataFrame(rerr[1])
CSV.write("data\\err11.csv", err)
err_meas = DataFrame(rerr[2])
CSV.write("data\\err_meas11.csv", err_meas)
conv = DataFrame(rconv[1])
CSV.write("data\\conv11.csv", conv)
conv_meas = DataFrame(rconv[2])
CSV.write("data\\conv_meas11.csv", conv_meas)
save("data\\output11.jld","err", err, "err_meas", err_meas, "conv", conv, "conv_meas", conv_meas, "raw", rraw[1], "raw_meas", rraw[2])
# data11 = load("data\\output11.jld")

# FCFS Data Collection Part 2
seed2 = 3637
n_runs = trialruns
# n_runs = size(settings,1)
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed2)
@time (rerr, rconv, rraw) = infer(param_inf)
err = DataFrame(rerr[1])
CSV.write("data\\err12.csv", err)
err_meas = DataFrame(rerr[2])
CSV.write("data\\err_meas12.csv", err_meas)
conv = DataFrame(rconv[1])
CSV.write("data\\conv12.csv", conv)
conv_meas = DataFrame(rconv[2])
CSV.write("data\\conv_meas12.csv", conv_meas)
save("data\\output12.jld","err", err, "err_meas", err_meas, "conv", conv, "conv_meas", conv_meas, "raw", rraw[1], "raw_meas", rraw[2])
# data12 = load("data\\output12.jld")

# LCFS Data Collection Part 1
seed = 8710
n_runs = trialruns
# n_runs = size(settings,1)
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
@time (rerr, rconv, rraw) = inferLCFS(param_inf)
err = DataFrame(rerr[1])
CSV.write("data\\err13.csv", err)
err_meas = DataFrame(rerr[2])
CSV.write("data\\err_meas13.csv", err_meas)
conv = DataFrame(rconv[1])
CSV.write("data\\conv13.csv", conv)
conv_meas = DataFrame(rconv[2])
CSV.write("data\\conv_meas13.csv", conv_meas)
save("data\\output13.jld","err", err, "err_meas", err_meas, "conv", conv, "conv_meas", conv_meas, "raw", rraw[1], "raw_meas", rraw[2])
# data13 = load("data\\output13.jld")
