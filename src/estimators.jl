# Deltamax Estimator
# Estimate the number of servers using the Deltamax algorithm
# Input: arrival order sorted by departure order
# Output: estimated number of servers
function c_deltamax(out_order::Vector{Int64})
  N = size(out_order, 1)
  currmax = out_order[1]
  deltamax = 0
  for i = 2:N
    nextmax = max(out_order[i], currmax)
    deltamax = max(deltamax, nextmax - currmax)
    currmax = nextmax
  end
  deltamax
end

# Order-based Estimator (Sets)
# Estimate the number of servers using the Order-based algorithm
# Input: arrival order sorted by departure order
# Output: estimated number of servers
function c_order_slow(out_order::Vector{Int64})
    N = size(out_order, 1)
    currmax = 0
    nservers = 0
    for i = 1:N
        currmax = max(out_order[i], currmax)
        custall = Set(1:currmax)
        custdeparted = Set(out_order[1:i])
        custservice = setdiff(custall, custdeparted)
        nservers = max(length(custservice) + 1, nservers)
    end
    nservers
end

# Order-based Estimator (Scalar)
# Estimate the number of servers using the Order-based algorithm
# Input: arrival order sorted by departure order
# Output: estimated number of servers
function c_order(out_order::Vector{Int64})
    N = size(out_order, 1)
    currmax = 0
    noccupied = 0
    nservers = 0
    for i = 1:N
        nextmax = max(out_order[i], currmax)
        if out_order[i] < currmax
            noccupied = noccupied - 1
        else
            noccupied = noccupied + nextmax - currmax - 1
        end
        currmax = nextmax
        nservers = max(nservers, noccupied + 1)
    end
    nservers
end

# WORK IN PROGRESS: Summary Estimator - fast version
# WARNING: PRODUCES INCORRECT OUTPUT
# function c_var_unf_fast(A::Vector{Float64}, D::Vector{Float64}, cmax::Int64)
#   # Initilaize data structures and functions
#   N = size(D, 1)
#   B̂ = zeros(N)
#   temp = zeros(N)
#   VS = zeros(cmax, 2) # var, unf
#   VS[1,:] = [-1,-1]
#   max_undef = 1
#   Dsort = sort(D)
#   ## Grid Search
#   for c = 2:cmax
#     B̂[1:c] = A[1:c]
#     for i = (c + 1):N
#       B̂[i] = max(A[i], Dsort[i - c])
#     end
#     Ŝ = D - B̂
#     if any(Ŝ .< 0)
#       VS[c, 1:2] = -1
#       max_undef = c
#     else
#       VS[c, 1] = var(Ŝ) * (N-1) # var
#     end
#   end
#   (r1, r2) = (99, 99) # initialize result to 99 for debugging
#   # var method
#   try
#     r1 = findin(VS[:, 1], minimum(filter(x -> x >= 0, VS[:, 1])))[1]
#   catch
#     r1 = 1
#   end
#   # uninformed method
#   r2 = max_undef + 1
#   (r1, r2, VS)
# end

# Summary Estimator
# Estimate the number of servers using the Variance and Uninformed methods
# Input: arrival times, departure times, and the max number of servers
# Output: estimated number of servers for each method
function c_var_unf(A::Vector{Float64}, D::Vector{Float64}, cmax::Int64)
  # Initilaize data structures and functions
  N = size(D, 1)
  B̂ = zeros(N)
  temp = zeros(N)
  VS = zeros(cmax, 2) # var, unf
  max_undef = 0
  ## Grid Search
  for c = 1:cmax
    B̂[1:c] = A[1:c]
    for i = (c + 1):N
      Dₖ = sort(D[1:(i - 1)])
      B̂[i] = max(A[i], Dₖ[i - c])
    end
    Ŝ = D - B̂
    if any(Ŝ .< 0)
      VS[c, 1:2] = -1
      max_undef = c
    else
      VS[c, 1] = var(Ŝ) * (N-1) # var
    end
  end
  (r1, r2) = (99, 99) # initialize result to 99 for debugging
  # var method
  try
    r1 = findin(VS[:, 1], minimum(filter(x -> x >= 0, VS[:, 1])))[1]
  catch
    r1 = 1
  end
  # uninformed method
  r2 = max_undef + 1
  (r1, r2, VS)
end
