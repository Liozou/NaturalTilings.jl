
include("TriangleIntersect.jl")
using .TriangleIntersect

"""
    TriangleIterator{T} <: AbstractVector{NTuple{3,SVector{3,T}}}

Iterator over a triangulation of a cycle.

If the cycle is a triangle itself, simply yield the cycle. Otherwise, yield the triangles
starting from the center of the cycle and going through two consecutive vertices.
"""
struct TriangleIterator{T,U<:AbstractVector{PeriodicVertex3D}} <: AbstractVector{NTuple{3,SVector{3,T}}}
    pge::PeriodicGraphEmbedding3D{T}
    cycle::U
    center::SVector{3,T}
end

Base.size(ti::TriangleIterator) = begin n = length(ti.cycle); n == 3 ? (1,) : (n,) end
Base.IndexStyle(::Type{TriangleIterator{T}}) where {T} = Base.IndexLinear()
function Base.getindex(ti::TriangleIterator{T}, i::Int) where T
    pos = ti.pge
    
    i == 1 && length(ti.cycle) == 3 && return (pos[ti.cycle[1]], pos[ti.cycle[2]], pos[ti.cycle[3]])
    i == length(ti.cycle) && return (ti.center, pos[ti.cycle[end]], pos[ti.cycle[1]])
    return (ti.center, pos[ti.cycle[i]], pos[ti.cycle[i+1]])
end

function intersect_triangleiter(line::Line, cycle::TriangleIterator)
    kind = 0
    above = false
    for (a,b,c) in cycle
        inter = intersect_triangle(line, Triangle(a, b, c))
        if inter.is_intersection
            if 0 < inter.dist < 1
                return 2 # the segment crosses the cycle
            elseif kind == 0
                kind = copysign(1, inter.dist)
            elseif (kind > 0 && inter.dist < 0) || (kind < 0 && inter.dist > 0)
                above = true # both ends of the line cross the cycle
            end
        end
    end
    return ifelse(above, -2, kind)
end

"""
    intersect_polyedra(tiling::Tiling{D}, tile::Vector{PeriodicVertex{D}}, line::Line) where D

Compute the intersection between a `line` ``(a,b)`` and a candidate `tile` given as a list
of its rings. Return:
* -2 if the segment appears to not intersect the tile but is in its convex hull.
* -1 if there is an intersection in the ``(-∞,a)`` part of the line.
* 0 if there is no intersection.
* 1 if there is an intersection in the ``(b,∞)`` part of the line.
* 2 if the segment crosses through a ring.
* 3 if the segment crosses through the tile
"""
function intersect_polyedra(tiling::Tiling{D}, line::Line, tile::Vector{PeriodicVertex{D}}) where D
    kind = 0
    for (x, ofs) in tile
        cycle = PeriodicGraphs.OffsetVertexIterator{D}(ofs, tiling.rings[x])
        iter = TriangleIterator(tiling.pge, cycle, tiling.ringcenters[x] + ofs)
        inter = intersect_triangleiter(line, iter)
        if inter != 0
            if abs(inter) > 1
                return inter
            elseif kind == 0
                kind = inter
            elseif (kind == 1 && inter < 0) || (kind == -1 && inter > 0)
                return 3
            end
        end
    end
    return kind
end
