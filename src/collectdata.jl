using CSV, DataFrames, JLD
using UnobservableQueues

# Initial settings
datapath = joinpath(@__DIR__, "..", "data")
fn = joinpath(datapath, "design.csv")
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
trialruns = 3 # how many runs to check data transfer

# FCFS Data Collection Part 1
seed = 8710
n_runs = trialruns
# n_runs = size(settings,1)
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
@time (rerr, rconv, rraw) = infer(param_inf)
err = DataFrame(rerr[1])
CSV.write(joinpath(datapath,"err21.csv"), err)
err_meas = DataFrame(rerr[2])
CSV.write(joinpath(datapath,"err_meas21.csv"), err_meas)
conv = DataFrame(rconv[1])
CSV.write(joinpath(datapath,"conv21.csv"), conv)
conv_meas = DataFrame(rconv[2])
CSV.write(joinpath(datapath,"conv_meas21.csv"), conv_meas)
save(joinpath(datapath,"output21.jld"),"err", err, "err_meas", err_meas, "conv", conv, "conv_meas", conv_meas, "raw", rraw[1], "raw_meas", rraw[2])
# data11 = load(joinpath(datapath,"output21.jld"))

# FCFS Data Collection Part 2
seed2 = 3637
n_runs = trialruns
# n_runs = size(settings,1)
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed2)
@time (rerr, rconv, rraw) = infer(param_inf)
err = DataFrame(rerr[1])
CSV.write(joinpath(datapath,"err22.csv"), err)
err_meas = DataFrame(rerr[2])
CSV.write(joinpath(datapath,"err_meas22.csv"), err_meas)
conv = DataFrame(rconv[1])
CSV.write(joinpath(datapath,"conv22.csv"), conv)
conv_meas = DataFrame(rconv[2])
CSV.write(joinpath(datapath,"\conv_meas22.csv"), conv_meas)
save(joinpath(datapath,"output22.jld"),"err", err, "err_meas", err_meas, "conv", conv, "conv_meas", conv_meas, "raw", rraw[1], "raw_meas", rraw[2])
# data12 = load(joinpath(datapath,"output22.jld"))

# LCFS Data Collection Part 1
seed = 8710
n_runs = trialruns
# n_runs = size(settings,1)
param_inf = Paraminf(settings, n_runs, n_methods, max_servers, obs_max, time_limit, ϵ, window, window_detail, step, seed)
@time (rerr, rconv, rraw) = inferLCFS(param_inf)
err = DataFrame(rerr[1])
CSV.write(joinpath(datapath,"err23.csv"), err)
err_meas = DataFrame(rerr[2])
CSV.write(joinpath(datapath,"err_meas23.csv"), err_meas)
conv = DataFrame(rconv[1])
CSV.write(joinpath(datapath,"conv23.csv"), conv)
conv_meas = DataFrame(rconv[2])
CSV.write(joinpath(datapath,"conv_meas23.csv"), conv_meas)
save(joinpath(datapath,"output23.jld"),"err", err, "err_meas", err_meas, "conv", conv, "conv_meas", conv_meas, "raw", rraw[1], "raw_meas", rraw[2])
# data13 = load(joinpath(datapath,"output23.jld"))
