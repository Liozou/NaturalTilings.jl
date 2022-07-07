
"""
    PseudoTriangleIterator{T} <: AbstractVector{NTuple{3,SVector{3,T}}}

Iterator over a triangulation of a cycle.

If the cycle is a triangle itself, simply yield the cycle. Otherwise, yield the triangles
starting from the center of the cycle and going through two consecutive vertices of the
cycle.
"""
struct PseudoTriangleIterator{T,U<:AbstractVector{PeriodicVertex3D}} <: AbstractVector{NTuple{3,SVector{3,T}}}
    pge::PeriodicGraphEmbedding3D{T}
    cycle::U
    center::SVector{3,T}
end

Base.size(pti::PseudoTriangleIterator) = begin n = length(pti.cycle); n == 3 ? (1,) : (n,) end
Base.IndexStyle(::Type{PseudoTriangleIterator{T}}) where {T} = Base.IndexLinear()
function Base.getindex(pti::PseudoTriangleIterator, i::Int)
    pos = pti.pge
    
    i == 1 && length(pti.cycle) == 3 && return (pos[pti.cycle[1]], pos[pti.cycle[2]], pos[pti.cycle[3]])
    i == length(pti.cycle) && return (pti.center, pos[pti.cycle[end]], pos[pti.cycle[1]])
    return (pti.center, pos[pti.cycle[i]], pos[pti.cycle[i+1]])
end

struct TriangleIterator{T,U<:AbstractVector{PeriodicVertex3D}} <: AbstractVector{Triangle{T}}
    pti::PseudoTriangleIterator{T,U}
end
TriangleIterator(pge, cycle, center) = TriangleIterator(PseudoTriangleIterator(pge, cycle, center))
Base.size(ti::TriangleIterator) = size(ti.pti)
Base.IndexStyle(::Type{TriangleIterator{T}}) where {T} = Base.IndexLinear()
Base.getindex(ti::TriangleIterator, i::Int) = Triangle(ti.pti[i])


struct OffsetTriangleIterator{T,U<:AbstractVector{Triangle{T}},S} <: AbstractVector{Triangle{T}}
    triangles::U
    ofs::S
end
Base.size(oti::OffsetTriangleIterator) = (length(oti.triangles),)
Base.IndexStyle(::Type{OffsetTriangleIterator{T,U,S}}) where {T,U,S} = Base.IndexLinear()
Base.getindex(oti::OffsetTriangleIterator, i::Int) = Triangle(oti.triangles[i], oti.ofs)


function intersect_triangleiter(line::Line, triangleiter::AbstractVector{Triangle{T}}) where T
    kind = 0
    above = false
    for t in triangleiter
        inter = intersect_triangle(line, t)
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
        inter = intersect_triangleiter(line, OffsetTriangleIterator(tiling.triangles[x], ofs))
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
