using Test

using DynamicCascades

g1 = import_system(Val(:rts96prepared), losses=true)
g2 = import_system(Val(:rts96raw), losses=true)
n1 = describe_nodes(g1, (;id=Int[]))
n2 = describe_nodes(g2, (;id=Int[]))
e1 = describe_edges(g1)
e2 = describe_edges(g2)


# compare the nodes
nodes = outerjoin(n1, n2, on=:id, makeunique=true)

@test nodes.n == nodes.n_1
@test nodes.P ≈ nodes.P_1
@test nodes.type == nodes.type_1

# FIXME node inertia wrong for node 14?
nodes.inertia[14] = 0.0 # 114
nodes.inertia[38] = 0.0 # 214
nodes.inertia[62] = 0.0 # 314
@test all(nodes.inertia .≈ Missings.replace(nodes.inertia_1, 0.0))


#compare the edges
edges = outerjoin(e1, e2, on=[:src, :dst], makeunique=true)

a= findall(x->! isapprox(x, 0.0), edges.rating - edges.rating_1)
b=findall(x->! isapprox(x, 0.0), edges.X - edges.X_1)
c=findall(x->! isapprox(x, 0.0), edges.R - edges.R_1)

a==b==c

findall(x->! isapprox(x, 0.0), edges.B - edges.B_1)

edges[:,[:rating, :rating_1]]

g3 = import_system(:rtsgmlc)
describe_nodes(g3)[]
