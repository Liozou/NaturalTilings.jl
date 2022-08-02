module NaturalTilings

using Graphs, StaticArrays, PeriodicGraphs, PeriodicGraphEmbeddings
using PeriodicGraphs: EdgeDict, DistanceRecord, IterativeGaussianEliminationDecomposition,
                      VertexPair, strong_erings, normalize_cycle!, retrieve_track!,
                      OffsetVertexIterator, convert_to_ering, gaussian_elimination!,
                      IterativeGaussianEliminationLength
using LinearAlgebra: norm, dot
using Statistics: mean

include("types.jl")
include("geometry.jl")
include("phantomedges.jl")
include("tiling.jl")
include("io.jl")

end
