# tests two different implementations of unf (and var)
# I think they are not the same, but I'm not positive

using BenchmarkTools
include("C:\\Users\\op\\Documents\\Julia Projects\\UnobservableQueue.jl\\src\\UnobservableQueue.jl")
c = 8 # number of servers
μ = 1 / c # service rate (single server)
λ = 0.9 # arrival rate (mean interarrival time is 1/λ)
adist = Exponential(1/λ) # arrival distribution
sdist = Exponential(1/μ) # service distribution (single server)
ncust = 1_000 # total number of customers generated
timelimit = 1_000 # time limit for simulation
seed = 8710 # rng seed
s1 = sim(ncust, timelimit, seed)
q1 = queue(c, adist, sdist)
df1 = runsim(s1, q1)
cmax = 19
A1 = df1[:atime]
D1 = df1[:dtime]
outorder1 = sort(df1, :dorder)[:aorder]

c_order(outorder1)
c_var_unf(A1,D1,cmax)
c_var_unf_slow(A1,D1,cmax)
c_var_unf_slow(A1,D1,cmax)[4]
c_var_unf(A1,D1,cmax)[4]
c_var_unf_slow(A1,D1,cmax)[4][1] - c_var_unf(A1,D1,cmax)[4][1] .|> abs |> sum
c_var_unf_slow(A1,D1,cmax)[4][2] - c_var_unf(A1,D1,cmax)[4][2]  .|> abs |> sum
c_var_unf_slow(A1,D1,cmax)[4][3] - c_var_unf(A1,D1,cmax)[4][3] .|> abs |> sum
# why is Ŝ different between slow and fast, but D and B̂ are same
# isn't shouldn't Ŝ = D - B̂ ??

@btime c_order(outorder1)
@btime c_var_unf(A1,D1,cmax)
@btime c_var_unf_slow(A1,D1,cmax)

# td = sort(df1,:dorder)[:dtime]
# ta = sort(df1,:aorder)[:dtime]
# to = zeros(998)
# for i = 1:998
#     to[i] = sort(ta[1:(i+2)])[i]
# end
# to - td[1:998] |> sum
#


# Check that both methods of sorting match
# different orderings for c hat < c true
# not sure why, or if it necessarily has to matter
D = df1[:dtime]
Dsort = sort(D)
rt = Vector{Bool}(cmax)
for c = 1:2
        results = zeros(1000,2)
        for i = c+1 : 1000
            results[i,1] = Dsort[i - c]
            results[i,2] = sort(D[1:(i - 1)])[i - c]
        end
        rt[c] = all(results[:,1] .== results[:,2])
end
rt

c = 8
results = zeros(1000,2)
for i = c+1 : 1000
    results[i,1] = Dsort[i - c]
    results[i,2] = sort(D[1:(i - 1)])[i - c]
end
all(results[:,1] .== results[:,2])
