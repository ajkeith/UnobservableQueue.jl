"""
A basic implementation of the order-based algorithm for queues with unobservable service.
"""
module UnobservableQueues

using SimJulia, ResumableFunctions
using Distributions, DataFrames, StatsBase

export
    Paraminf,
    infer,
    inferLCFS

include("queue.jl")
include("estimators.jl")
include("inference.jl")
include("lcfs.jl")

end
