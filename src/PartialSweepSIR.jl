module PartialSweepSIR

import Base: vec, push!, length, size, getindex

using Chain
using LinearAlgebra
using OrdinaryDiffEq
using Parameters
using StatsBase


const PSS = PartialSweepSIR
export PSS

include("objects.jl")
include("dynamics.jl")
include("tools.jl")
include("analytics.jl")

end
