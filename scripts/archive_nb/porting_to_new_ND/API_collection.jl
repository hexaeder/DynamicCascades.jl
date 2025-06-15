"""
Things I've tried out during porting to new ND.
"""
###
### API collection
### 
using SymbolicIndexingInterface

i = 1 # vertex index
vidxs(nw, i, "ω") # returns symbolic index of vertex IF it has a state ω
vidxs(nw, i, :ω)  # returns symbolic index of vertex REGARDLESS if it has a ω state or not

vidxs(1:10, :ω) 

#= Returns index of :ω internal STATE of vertex i. Here one needs to keep in mind that the states are NOT ordered by
the vertices but by the different vertex models. For the RTS this means first the SwingDynLoad nodes then
the DynLoad nodes. =#
SII.variable_index(nw, VIndex(i, :ω)) 

s0.v[i,:ω] # returns ω-value of vertex i


# Get index instead of symbolic index. Sollte man nicht in hot loops verwenden wg. performance
map(idx -> idx.compidx, vidxs(nw, :, "ω"))

vidxs(nw, :, "θ")
vidxs(nw, :, "ω")
vpidxs(nw, :, "Pmech")

vidxs(nw, 13, "ω")[1].compidx # 13
vidxs(nw, 13, "ω")[1].subidx # ω

p.v[map(idx -> idx.compidx, vidxs(nw, :, "ω")), :Pmech]


using SymbolicIndexingInterface: SymbolicIndexingInterface as SII
ω_idxs = Int[] 
for i in 1:nv(nw)
    try
        idx = SII.variable_index(nw, VIndex(i, :ω))
        push!(ω_idxs, idx)
    catch
        nothing
    end
end
ω_idxs