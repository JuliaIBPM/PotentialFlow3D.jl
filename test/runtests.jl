using PotentialFlow3D
using Test

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All"
    include("types.jl")
    include("operations.jl")
  end