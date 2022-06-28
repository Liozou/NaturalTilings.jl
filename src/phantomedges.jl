
"""
    ref_vertexpair(kp::EdgeDict{D}, x::Int) where D

Return `(y, ofs)` where `y::Int` refers to the same edge as `x` offset by `-ofs`, such that
all representatives of `x` yield the same `y`.
"""
function ref_vertexpair(kp::EdgeDict{D}, x::Int) where D
    (v1, ofs1), (v2, ofs2) = kp[x]
    iszero(ofs1) && return x, ofs1
    return kp[(PeriodicVertex{D}(v1), PeriodicVertex{D}(v2, ofs2 .- ofs1))], ofs1
end

struct ERingAttributions{D}
    erings::Vector{Vector{Int}}
    attrs::Vector{Vector{Tuple{Int,Int}}}
    kp::EdgeDict{D}

    function ERingAttributions(erings, kp::EdgeDict{D}) where D
        n = length(kp.direct)
        attrs = [Tuple{Int,Int}[] for _ in 1:n]
        for (i, er) in enumerate(erings)
            for (j, x) in enumerate(er)
                ref, = ref_vertexpair(kp, x)
                @assert ref ≤ n
                push!(attrs[ref], (i, j))
            end
        end
        return new{D}(erings, attrs, kp)
    end
end

struct ERingIncluding{D} <: AbstractVector{Vector{Int}}
    eras::ERingAttributions{D}
    i::Int
    ofs::SVector{D,Int}
    buffer::Vector{Int}
    ERingIncluding(eras::ERingAttributions{D}, i, ofs) where {D} = new{D}(eras, i, ofs, Int[])
end

function Base.getindex(eri::ERingIncluding{D}, i::Int) where {D}
    newring_idx, idx = eri.eras.attrs[eri.i][i]
    newering = eri.eras.erings[newring_idx]
    k = newering[idx]
    n = length(newering)
    resize!(eri.buffer, n)
    if k == i
        unsafe_copyto!(eri.buffer, 1, newering, 1, n)
    else
        kp = eri.eras.kp
        (_, ofs_ref2), = kp[k]
        ofs_ref = ofs_ref2 .- eri.ofs
        for (j, e) in enumerate(newering)
            (x1, ofs1), (x2, ofs2) = kp[e]
            eri.buffer[j] = get!(kp, (PeriodicVertex{D}(x1, ofs1 .- ofs_ref), PeriodicVertex{D}(x2, ofs2 .- ofs_ref)))
        end
        sort!(eri.buffer)
    end
    return eri.buffer
end
Base.size(eri::ERingIncluding) = size(eri.eras.attrs[eri.i])
Base.IndexStyle(::Type{ERingIncluding{D}}) where {D} = Base.IndexLinear()

Base.@propagate_inbounds function Base.getindex(eras::ERingAttributions, i::Integer)
    j, ofs = ref_vertexpair(eras.kp, i)
    @boundscheck checkbounds(eras.attrs, j)
    ERingIncluding(eras, i, ofs)
end
Base.size(eras::ERingAttributions) = size(eras.attrs)
Base.IndexStyle(::Type{ERingAttributions{D}}) where {D} = Base.IndexLinear()


function canonical_ering!(ering::Vector{Int}, kp::EdgeDict{D}) where D
    lenc = length(ering)
    ((_, ofs), x2), fst = findmin(j -> kp[j], ering)
    if !iszero(ofs)
        for i in 1:lenc
            (v1, ofs1), (v2, ofs2) = ering[i]
            ering[i] = kp[(PeriodicVertex{D}(v1, ofs1 .- ofs), PeriodicVertex{D}(v2, ofs2 .- ofs))]
        end
        sort!(ering)
    end
    nothing
end

function find_common_vertex((xa, xb), (ya, yb))
    (xa == ya || xa == yb) && return xa
    return xb
end

function identify_junction(er1, er2, kp)
    len1 = length(er1)
    len2 = length(er2)
    start1 = 1
    _start2 = findfirst(==(er1[1]), er2)
    local start2::Int
    if _start2 isa Int
        start2 = _start2
    else
        start1 = 1 + length(er1)÷2
        start2 = findfirst(==(er1[start1]), er2)
    end
    after1 = er1[start1+1]
    after2 = er2[mod1(start2+1, len2)]
    before2 = er2[mod1(start2-1, len2)]
    direct = signbit(after1 == after2 || after1 != before2)
    ia = mod1(start1 + 1, len1)
    ja = mod1(start2 + direct, len2)
    while er1[ia] == er2[ja]
        ia = mod1(ia + 1, len1)
        ja = mod1(ja + direct, len2)
    end
    ib = mod1(start1 - 1, len1)
    jb = mod2(start2 - direct, len2)
    while er1[ib] == er2[jb]
        ib = mod1(ib - 1, len1)
        jb = mod1(jb - direct, len2)
    end
    _junctiona = find_common_vertex(er1[ia], er2[ja])
    _junctionb = find_common_vertex(er1[ib], er2[jb])
    (xa, ofsa), (xb, ofsb) = minmax(_junctiona, _junctionb)
    return PeriodicEdge{D}(xa, xb, ofsb .- ofsa), ofsa, ia, ib, ja, jb
end

function add_phantomedges!(erings::Vector{Vector{Int}}, kp::EdgeDict{D}) where D
    n = length(erings)
    eringdict = Dict{Vector{Int},Int}(x => i for (i,x) in enumerate(erings))
    junctions = PeriodicEdge{D}[] # the list of newly added edges
    junctions_dict = Dict{PeriodicEdge{D},Int}() # reverse map to junctions
    junctions_per_ring = [Tuple{PeriodicVertex{D},Int,Int}[] for _ in 1:n]
    # junctions_per_ring[i] is the list of tuples (x, j, k) indicating that junction x is
    # present between vertices j and k of rings[i].
    known_erings = Dict{Vector{Int},Int}()
    for (i,e) in enumerate(erings)
        e2 = copy(e)
        canonical_ering!(e2, kp)
        known_erings[e2] = i
    end
    buffer = Int[]
    eras = ERingAttributions(erings, kp)
    for (i1, er1) in enumerate(erings)
        len1 = length(er1)
        @show len1
        for x1 in er1
            @show x1
            @show eras[x1]
            
            for er2 in eras[x1]
                len2 = length(er2)
                abs(len2 - len1) > 1 && continue
                PeriodicGraphs.symdiff_cycles!(buffer, er1, er2)
                length(buffer) ≤ min(len1, len2) || continue
                canonical_ering!(buffer, kp)
                i = get(known_erings, buffer, 0)
                i == 0 && continue
                junction, ofs, ia, ib , ja, jb = identify_junction(er1, er2, kp)
                junction_idx = get!(junctions_dict, junction, length(junctions)+1)
                junction_idx > length(junctions) && push!(junctions, junction)

            end
        end
    end
end