# TriangleIntersect.jl adapted from https://github.com/JuliaGeometry/TriangleIntersect.jl

# MIT "Expat" License:
# Copyright (c) 2014: Ariel Keselman.
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


module TriangleIntersect

import Base: -, *, +, /#, intersect, ≈

export Point, Triangle, Ray, Line, Intersection, intersect_triangle

struct Point{T}
    x::T
    y::T
    z::T
end
Point((x, y, z)) = Point(x, y, z)

*(p::Point{T}, n::Number) where {T} = Point{T}(p.x*n, p.y*n, p.z*n)
*(n::Number, p::Point) = p*n
*(p1::Point, p2::Point) = p1.x*p2.x+p1.y*p2.y+p1.z*p2.z
-(p1::Point{T}, p2::Point{T}) where {T} = Point{T}(p1.x-p2.x, p1.y-p2.y, p1.z-p2.z)
-(p::Point{T}) where {T} = Point{T}(-p.x, -p.y, -p.z)
+(p1::Point{T}, p2::Point{T}) where {T} = Point{T}(p1.x+p2.x, p1.y+p2.y, p1.z+p2.z)
cross(p1::Point{T}, p2::Point{T}) where {T} = Point{T}(p1.y*p2.z-p1.z*p2.y, -p1.x*p2.z+p1.z*p2.x, p1.x*p2.y-p1.y*p2.x)
/(p::Point{T}, n::Number) where {T} = Point{T}(p.x/n, p.y/n, p.z/n)
# unitize(p::Point) = p/(p*p)
# ≈(p1::Point, p2::Point) = p1.x ≈ p2.x && p1.y ≈ p2.y && p1.z ≈ p2.z

struct Triangle{T}
    a::Point{T}
    v1::Point{T}
    v2::Point{T}
    normal::Point{T}
    v1v1::T
    v2v2::T
    v1v2::T
    denom::T
    function Triangle(a::Point{T}, b::Point{T}, c::Point{T}) where {T}
        v1 = b-a
        v2 = c-a
        normal = cross(v1, v2)
        v1v1 = v1*v1
        v2v2 = v2*v2
        v1v2 = v1*v2
        denom = v1v2*v1v2 - v1v1*v2v2
        new{T}(a, v1, v2, normal, v1v1, v2v2, v1v2, denom)
    end
end
Triangle(a, b, c) = Triangle(Point(a), Point(b), Point(c))

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
        p = b-a
        unit = inv(p*p)
        new{T}(a, b, unit, p*unit)
    end
end
Line(a, b) = Line(Point(a), Point(b))

function intersect_triangle(l::Line{T}, t::Triangle) where {T}
    denom = t.normal*l.direction
    iszero(denom) && return Intersection{T}()
    ri = t.normal*(t.a - l.a) / denom
    origin, direction = ri ≥ 0 ? (l.a, l.direction) : (l.b, -l.direction)
    plane_intersection =  abs(ri) * direction + origin
    w = plane_intersection - t.a
    wv1 = w*t.v1
    wv2 = w*t.v2
    s_intersection = (t.v1v2*wv2 - t.v2v2*wv1) / t.denom
    s_intersection < 0 && return Intersection{T}()
    s_intersection > 1 && return Intersection{T}()
    t_intersection = (t.v1v2*wv1 - t.v1v1*wv2) / t.denom
    t_intersection < 0 && return Intersection{T}()
    t_intersection > 1 && return Intersection{T}()
    s_intersection + t_intersection > 1 && return Intersection{T}()
    # @assert t.a + s_intersection*t.v1+t_intersection*t.v2 ≈ l.a + ri*l.direction
    return Intersection{T}(ri*l.unit, true)
end

end # module TriangleIntersect
