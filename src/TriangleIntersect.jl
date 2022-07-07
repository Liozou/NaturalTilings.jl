# TriangleIntersect.jl adapted from https://github.com/JuliaGeometry/TriangleIntersect.jl

# MIT "Expat" License:
# Copyright (c) 2014: Ariel Keselman.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


module TriangleIntersect

using StaticArrays
export Triangle, Ray, Line, Intersection, intersect_triangle, shortline

const Point{T} = SVector{3,T} # type alias

dot(p1::Point, p2::Point) = p1.x*p2.x+p1.y*p2.y+p1.z*p2.z
cross(p1::Point{T}, p2::Point{T}) where {T} = Point{T}(p1.y*p2.z-p1.z*p2.y, -p1.x*p2.z+p1.z*p2.x, p1.x*p2.y-p1.y*p2.x)

struct Triangle{T}
    a::Point{T}
    v1::Point{T}
    v2::Point{T}
    normal::Point{T}
    v1v1::T
    v2v2::T
    v1v2::T
    denom::T
end

function Triangle(a::Point{T}, b::Point{T}, c::Point{T}) where {T}
    v1 = b-a
    v2 = c-a
    normal = cross(v1, v2)
    v1v1 = dot(v1, v1)
    v2v2 = dot(v2, v2)
    v1v2 = dot(v1, v2)
    denom = v1v2*v1v2 - v1v1*v2v2
    Triangle{T}(a, v1, v2, normal, v1v1, v2v2, v1v2, denom)
end
Triangle((a, b, c)) = Triangle(Point(a), Point(b), Point(c))

function Triangle(t::Triangle{T}, ofs::Point{S}) where {T,S}
    Triangle{T}(t.a + ofs, t.v1, t.v2, t.normal, t.v1v1, t.v2v2, t.v1v2, t.denom)
end

function Base.show(io::IO, triangle::Triangle)
    print(io, Triangle, "(((")
    join(io, triangle.a, ", ")
    print(io, "), (")
    join(io, triangle.a + triangle.v1, ", ")
    print(io, "), (")
    join(io, triangle.a + triangle.v2, ", ")
    print(io, ")))")
end

struct Intersection{T}
    dist::T # intersecting distance
    is_intersection::Bool
end
Intersection{T}() where {T} = Intersection{T}(zero(T), false)

# struct Ray{T}
#     origin::Point{T}
#     direction::Point{T}
#     Ray{T}(a::Point{T}, direction::Point{T}) where {T} = new{T}(a, direction)
# end
# Ray(a::Point{T}, b::Point{T}) where {T} = Ray{T}(a, unitize(b))

# function intersect(r::Ray{T}, t::Triangle{T}) where {T}
#     denom = t.normal*r.direction
#     iszero(denom) && return Intersection{T}()
#     ri = t.normal*(t.a - r.origin) / denom
#     ri < 0 && return Intersection{T}()
#     plane_intersection =  ri * r.direction + r.origin
#     w = plane_intersection - t.a
#     wv1 = w*t.v1
#     wv2 = w*t.v2
#     s_intersection = (t.v1v2*wv2 - t.v2v2*wv1) / t.denom
#     s_intersection <= 0 && return Intersection{T}()
#     s_intersection >= 1 && return Intersection{T}()
#     t_intersection = (t.v1v2*wv1 - t.v1v1*wv2) / t.denom
#     t_intersection <= 0 && return Intersection{T}()
#     t_intersection >= 1 && return Intersection{T}()
#     s_intersection + t_intersection >= 1 && return Intersection{T}()
#     # Intersection(t.a + s_intersection*t.v1+t_intersection*t.v2, ri, true)
#     Intersection{T}(ri, true)
# end

struct Line{T}
    a::Point{T}
    b::Point{T}
    unit::T
    direction::Point{T}
    function Line(a::Point{T}, b::Point{T}) where {T}
        p = b - a
        unit = inv(dot(p, p))
        new{T}(a, b, unit, unit*p)
    end
end
Line(a, b) = Line(Point(a), Point(b))
Line(line::Line{T}, ofs) where {T} = Line(Point{T}(line.a + ofs), Point{T}(line.b + ofs))

function Base.show(io::IO, line::Line)
    print(io, Line, "((")
    join(io, line.a, ", ")
    print(io, "), (")
    join(io, line.b, ", ")
    print(io, "))")
end

function shortline(a::Point, b::Point)
    d = b - a
    δ = d/(1000*dot(d, d))
    return Line(a+δ, b-δ)
end
shortline(a, b) = shortline(Point(a), Point(b))

function intersect_triangle(l::Line{T}, t::Triangle) where {T}
    denom = dot(t.normal, l.direction)
    iszero(denom) && return Intersection{T}()
    ri = dot(t.normal, (t.a - l.a)) / denom
    origin = l.a
    direction = l.direction
    sign = false
    if ri < 0
        origin = l.b
        direction = -l.direction
        ri = - dot(t.normal, (t.a - l.b)) / denom
        sign = true
        @assert ri ≥ 0
    end
    plane_intersection =  ri * direction + origin
    w = plane_intersection - t.a
    wv1 = dot(w, t.v1)
    wv2 = dot(w, t.v2)
    s_intersection = (t.v1v2*wv2 - t.v2v2*wv1) / t.denom
    s_intersection < 0 && return Intersection{T}()
    s_intersection > 1 && return Intersection{T}()
    t_intersection = (t.v1v2*wv1 - t.v1v1*wv2) / t.denom
    t_intersection < 0 && return Intersection{T}()
    t_intersection > 1 && return Intersection{T}()
    s_intersection + t_intersection > 1 && return Intersection{T}()
    if !(t.a + s_intersection*t.v1 + t_intersection*t.v2 ≈ origin + ri*direction)
        @show t
        @show l
        @show l.a + ri*l.direction
        @show origin + abs(ri)*direction
        @show denom, ri, wv1, wv2, s_intersection, t_intersection
        @assert false
    end
    rescaled_ri = ri*l.unit
    return Intersection{T}(ifelse(sign, one(T) - rescaled_ri, rescaled_ri), true)
end

end # module TriangleIntersect
