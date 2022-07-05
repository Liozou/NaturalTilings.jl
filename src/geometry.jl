
include("TriangleIntersect.jl")

cycle_center(cycle::Vector{T}, positions) where {T} = sum(positions[x] for x in cycle; init=zero(T)) / length(cycle)

"""
    TriangleIterator{T} <: AbstractVector{NTuple{3,SVector{3,T}}}

Iterator over a triangulation of a cycle.

If the cycle is a triangle itself, simply yield the cycle. Otherwise, yield the triangles
starting from the center of the cycle and going through two consecutive vertices.
"""
struct TriangleIterator{T} <: AbstractVector{NTuple{3,SVector{3,T}}}
    cycle::Vector{PeriodicVertex3D}
    positions::Vector{SVector{3,T}}
    center::SVector{3,T}
end

Base.size(ti::TriangleIterator) = begin n = length(ti.cycle); n == 3 ? (1,) : (n,) end
Base.IndexStyle(::Type{TriangleIterator{T}}) where {T} = Base.IndexLinear()
function Base.getindex(ti::TriangleIterator{T}, i::Int) where T
    pos = ti.positions
    i == 1 && length(ti.cycle) == 3 && return (pos[ti.cycle[1]], pos[ti.cycle[2]], pos[ti.cycle[3]])
    i == length(ti.cycle) && return (pos[ti.center], pos[ti.cycle[end]], pos[ti.cycle[1]])
    return (pos[ti.center], pos[ti.cycle[i]], pos[ti.cycle[i+1]])
end

