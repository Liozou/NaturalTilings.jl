module NaturalTilings

using Graphs, StaticArrays, PeriodicGraphs, PeriodicGraphEmbeddings
using PeriodicGraphs: EdgeDict, DistanceRecord, IterativeGaussianEliminationDecomposition,
                      VertexPair, strong_erings, normalize_cycle!, retrieve_track!,
                      convert_to_ering, gaussian_elimination!,
                      IterativeGaussianEliminationLength
using LinearAlgebra: norm, dot

include("types.jl")
include("geometry.jl")
include("phantomedges.jl")
include("tiling.jl")
include("io.jl")

end
