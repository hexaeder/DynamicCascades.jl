using Test

using DynamicCascades
using MetaGraphs
using Missings
using Graphs
using Unitful

@testset "comparison between differen rts systems" begin
    g1 = import_system(:rtsgmlc)
    n1 = describe_nodes(g1)
    e1 = describe_edges(g1)

    g2 = import_system(:rts96raw)
    n2 = describe_nodes(g2)
    e2 = describe_edges(g2)

    g3 = import_system(:rts96prepared)
    n3 = describe_nodes(g3)
    e3 = describe_edges(g3)

    @test n1.n == n2.n == n3.n
    @test n1.id == n2.id == n3.id
    n1.type == n2.type == n3.type
    @test n1.Vbase == n2.Vbase == n3.Vbase

    #inertia
    inert1 = Missings.replace(n1.H, 0.0u"MJ/MW") |> collect
    inert2 = Missings.replace(n2.H, 0.0u"MJ/MW") |> collect
    inert3 = Missings.replace(n3.H, 0.0u"MJ/MW") |> collect
    # Anton sets interita for "Sync Cond" nodes to 1.0
    inert3[14] = 0.0u"MJ/MW" # 114
    inert3[38] = 0.0u"MJ/MW" # 214
    inert3[62] = 0.0u"MJ/MW" # 314

    @test inert2 ≈ inert3 atol=1e-8u"MJ/MW" # same with those 'fixed'

    Δinert = inert1 .- inert2
    idx = findall(x->!isapprox(x, 0u"MJ/MW"), Δinert)
    @info "Difference in inert to GMLC" idx .=> Δinert[idx]

    # there are some expected differences in power between 2/3 and 1
    @test n2.P ≈ n3.P
    ΔP = n1.P - n2.P
    idx = findall(x->!isapprox(x, 0.0), n1.P - n2.P)
    @info "Difference in P to GMLC" idx .=> ΔP[idx]

    # comparison of edges
    @test e1.src == e2.src == e3.src
    @test e1.dst == e2.dst == e3.dst

    # e1 and e2 is identical for edge data
    @test e1.R ≈ e2.R
    @test e1.X ≈ e2.X
    @test e1.rating ≈ e2.rating
    # idx = [27, 33, 34, 35, 63, 69, 70, 71, 98, 104, 105, 106]
    # a = findall(x->!isapprox(x, 0.0), e2.rating - e3.rating)
    # b = findall(x->!isapprox(x, 0.0), e2.X - e3.X)
    # c = findall(x->!isapprox(x, 0.0), e2.R - e3.R)
    # @test idx == a == b == c
end

@testset "graph distance" begin
    g = path_graph(7)
    dist = edge_distance_to(g, 1) do g, v1, v2
        d = length(a_star(g, v1, v2))
    end
    @test dist == [0,1,2,3,4,5]

    dist = edge_distance_to(g, 3) do g, v1, v2
        d = length(a_star(g, v1, v2))
    end
    @test dist == [2,1,0,1,2,3]

    network =  import_system(:kaiser2020, gen_γ=0.8, load_τ=0.1)
    src = edgeidx(network, 2, 17)
    src=2
    d1 = resistance_distances(network, src)
    d2 = naive_distances(network, src)

    # elabels = ["" for i in 1:ne(network)]
    # elabels[src] = "src"
    # graphplot(network; edge_color=d, elabels)
    # graphplot(network; edge_color=d, elabels)
end

@testset "time series stuff" begin
    using DynamicCascades: remove_zero_tail!
    x = [1, 2, 3, 4, 5]
    y = [2, 4, 1, 0, 0]
    @test remove_zero_tail!(x, y) == (x[1:3], y[1:3])
end
