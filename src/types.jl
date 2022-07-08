
include("TriangleIntersect.jl")
using .TriangleIntersect

export Tiling

function cycle_center(cycle, pge::PeriodicGraphEmbedding{D,T}) where {D,T}
    sum(pge[x] for x in cycle; init=zero(SVector{D,T})) / length(cycle)
end

struct Tiling{D,T}
    pge::PeriodicGraphEmbedding{D,T}
    g::PeriodicGraph{D} # pge.g possibly including new edges
    rings::Vector{Vector{PeriodicVertex{D}}} # each sublist is the list of vertices
    erings::Vector{Vector{Int}} # each sublist is the list of edge indices in kp of a ring
    tiles::Vector{Vector{PeriodicVertex{D}}} # each sublist is a list of indices of rings
    tileedges::Vector{Vector{Int}} # each sublist is the list of edge indices in kp of a tile
    tilevertices::Vector{Vector{PeriodicVertex{D}}} # each sublist is the list of vertices in a tile
    tiledict::Dict{Vector{Int},Int} # for each tile (represented as its unique edges), its index
    tilesofring::Vector{Tuple{PeriodicVertex{D},PeriodicVertex{D}}} # for each ring, the two tiles attached to it
    ringsofedge::Dict{PeriodicEdge{D},Vector{PeriodicVertex{D}}} # for each edge, the list of associated ring
    rgraph::PeriodicGraph{D} # The graph of rings: two rings are bonded if they share an edge
    ringcenters::Vector{SVector{3,T}} # The center of each ring, if D == 3
    triangles::Vector{Vector{Triangle{T}}} # Triangulation of each ring, if D == 3
    kp::EdgeDict{D} # correspondance between edges and their index
end

function Tiling(pge::PeriodicGraphEmbedding{D,T}, g, rings, erings, kp) where {D,T}
    n = length(erings)
    ringsofedge = Dict{PeriodicEdge{D},Vector{PeriodicVertex{D}}}()
    rgraph = PeriodicGraph{D}(n)
    multipleneighbors = PeriodicEdge{D}[]
    for (i, ering) in enumerate(erings)
        for _e in ering
            u1, u2 = kp[_e]
            @assert (u1, u2) == minmax(u1, u2)
            e = PeriodicEdge{D}(u1.v, u2.v, u2.ofs .- u1.ofs)
            ringsofe = get!(ringsofedge, e, PeriodicVertex{D}[])
            for idx in ringsofe
                newofs = u1.ofs - idx.ofs
                if idx.v != i || !iszero(newofs)
                    new_e = PeriodicEdge{D}(i, PeriodicVertex{D}(idx.v, newofs))
                    add_edge!(rgraph, new_e) || push!(multipleneighbors, new_e)
                end
            end
            push!(ringsofe, PeriodicVertex{D}(i, u1.ofs))
        end
    end
    # for e in unique!(sort!(multipleneighbors))
    #     rem_edge!(rgraph, e)
    # end

    tiles = Vector{PeriodicVertex{D}}[]
    tileedges = Vector{Int}[]
    tilevertices = Vector{PeriodicVertex{D}}[]
    tilesofring = [(PeriodicVertex{D}(0),PeriodicVertex{D}(0)) for _ in 1:n]
    tiledict = Dict{Vector{Int},Int}()
    ringcenters = D == 3 ? [cycle_center(r, pge) for r in rings] : SVector{3,T}[]
    triangles = D == 3 ? [collect(TriangleIterator(pge, r, center)) for (r, center) in zip(rings, ringcenters)] : Vector{Triangle{T}}[]
    return Tiling{D,T}(pge, g, rings, erings, tiles, tileedges, tilevertices, tiledict, tilesofring, ringsofedge, rgraph, ringcenters, triangles, kp)
end


"""
    countconsecutiveunique(list::Vector)

Return the list of pair `(count, item)` where `count` is the number of times `item` appears
consecutively at this point of the `list`.
"""
function countconsecutiveunique(l::Vector{T}) where T
    ret = Tuple{Int,T}[]
    isempty(l) && return ret
    count = 1
    lst = first(l)
    for x in Iterators.rest(l, 2)
        if x == lst
            count += 1
        else
            push!(ret, (count, lst))
            count = 1
            lst = x
        end
    end
    push!(ret, (count, lst))
    return ret
end

# const superscriptdigits = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')
# function tosuperscript(io::IO, x::Int)
#     ds = digits(x)
#     for d in Iterators.reverse(ds)
#         print(io, superscriptdigits[d+1])
#     end
#     nothing
# end

function Base.show(io::IO, ::MIME"text/plain", t::Tiling{D}) where D
    if isempty(t.tiles)
        print(io, "No tiling found")
        return
    end
    tiles = [sort!([length(t.rings[x.v]) for x in _t]) for _t in t.tiles]
    sort!(tiles)
    sort!(tiles; by=length)
    ctiles = countconsecutiveunique([countconsecutiveunique(t) for t in tiles])
    d = gcd([c for (c,_) in ctiles])
    for (i, (c, tile)) in enumerate(ctiles)
        i == 1 || print(io, " + ")
        c == d || print(io, div(c, d))
        print(io, '[')
        for (j, (c2, x)) in enumerate(tile)
            j == 1 || print(io, ", ")
            print(io, x)
            if c2 != 1
                print(io, '^', c2)
            end
        end
        print(io, ']')
    end
    if d != 1
        print(io, " (×", d, ')')
    end
    nothing
end


"""
    TilingAroundCycle{D}

Internal struct containing the information on the incremental tiling built around a given
cycle. When the final tiling is built, each essential cycle should separate two tiles.
"""
struct TilingAroundCycle{D}
    gauss::IterativeGaussianEliminationDecomposition
    encountered::Set{PeriodicVertex{D}}
    Q::Vector{Tuple{PeriodicVertex{D},Int}}
    tracks::Vector{Vector{Vector{Int32}}}
    restart::Base.RefValue{Int}
end
function TilingAroundCycle{D}(i) where D
    gauss = IterativeGaussianEliminationDecomposition()
    encountered = Set{PeriodicVertex{D}}((PeriodicVertex{D}(i),))
    tracks = Vector{Vector{Vector{Int32}}}(undef, 1)
    Q = [(PeriodicVertex{D}(i), 0)]
    return TilingAroundCycle{D}(gauss, encountered, Q, tracks, Ref(1))
end

"""
    ering_at_pos(tiling::Tiling{D}, (i, ofs)::PeriodicVertex{D})

Retrieve the ering in `tiling` corresponding to ring number `i` offset by `ofs`.
"""
function ering_at_pos(tiling::Tiling{D}, (i, ofs)::PeriodicVertex{D}) where D
    return convert_to_ering(tiling.rings[i], tiling.kp, ofs)
end

function minimal_track!(gauss, tracks)
    track = retrieve_track!(gauss)
    sort!(track)
    buffer = Int32[]
    @label restart
    for i in track
        if isassigned(tracks, i)
            for t in tracks[i]
                PeriodicGraphs.symdiff_cycles!(buffer, track, t)
                if length(buffer) < length(track) # TODO: handle ≤ ?
                    track, buffer = buffer, track
                    @goto restart
                end
            end
        end
    end
    return track
end

"""
    explore_around_cycle!(tac::TilingAroundCycle{D}, tiling::Tiling{D}, untilfirstfound=false) where D

Internal utility to build an incremental tiling around the given cycle. Each repeated call
will increment the exploration radius, which is the maximal distance between the starting
cycle and any other considered for tiling in the cycle graph.
"""
function explore_around_cycle!(tac::TilingAroundCycle{D}, tiling::Tiling{D}, untilfirstfound=false) where D
    Q = tac.Q
    restart = tac.restart[]
    restart > length(Q) && return Vector{PeriodicVertex{D}}[]
    maxdist = untilfirstfound ? typemax(Int) : isone(restart) ? 1 : last(Q[restart]) + 1
    tiles = Vector{PeriodicVertex{D}}[]
    tracks = tac.tracks
    for (u, dist) in Iterators.rest(Q, restart)
        dist > maxdist && break
        restart += 1
        if gaussian_elimination!(tac.gauss, ering_at_pos(tiling, u)) # a sum of previously encountered rings is empty
            track = minimal_track!(tac.gauss, tracks)
            if first(track) == 1 # otherwise the tile does not contain the current cycle.
                if untilfirstfound
                    maxdist = dist
                end
                push!(tiles, sort!([first(Q[x]) for x in track]))
            else
                resize!(tracks, restart-1)
                for i in track
                    if isassigned(tracks, i)
                        push!(tracks[i], track)
                    else
                        tracks[i] = [track]
                    end
                end
            end
        end
        for x in neighbors(tiling.rgraph, u)
            x ∈ tac.encountered && continue
            push!(tac.encountered, x)
            push!(Q, (x, dist+1))
        end
    end
    tac.restart[] = restart
    return tiles
end

function unique_vertices(tile::Vector{PeriodicVertex{D}}, tiling::Tiling{D}) where D
    doublevertices = reduce(vcat, PeriodicGraphs.OffsetVertexIterator{D}(ofs, tiling.rings[x]) for (x, ofs) in tile; init=PeriodicVertex{D}[])
    return unique!(sort!(doublevertices))
end
