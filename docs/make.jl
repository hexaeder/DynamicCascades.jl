using GraphMakie
using DynamicCascades
using Documenter
using Literate
using CairoMakie

# preload the deps from the examples to supress precompilation output in docs
using NetworkDynamics
using Graphs

DocMeta.setdocmeta!(DynamicCascades, :DocTestSetup, :(using DynamicCascades); recursive=true)

# generate examples
script_dir = abspath(joinpath(@__DIR__, "..", "scripts"))
outdir = joinpath(@__DIR__, "src", "generated")
isdir(outdir) && rm(outdir, recursive=true)
mkpath(outdir)

@info "Create Mardown files from scripts"
scripts = ["rtsgmlc_plot.jl", "squaregrid.jl", "insulators.jl"]
for script in scripts
    path = joinpath(script_dir, script)
    Literate.markdown(path, outdir)
end

makedocs(; modules=[DynamicCascades], authors="Hans Würfel",
         repo="https://github.com/hexaeder/DynamicCascades.jl/blob/{commit}{path}#{line}",
         sitename="DynamicCascades.jl",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true",
                                canonical="https://wuerfel.io/DynamicCascades.jl", assets=String[]),
         pages=["Home" => "index.md",
                "Examples" => [
                    "RTS-GMLC" => "generated/rtsgmlc_plot.md",
                    "Disturbance in Grid" => "generated/squaregrid.md",
                    "Topological insulator" => "generated/insulators.md",
                ]
                ])

# if gh_pages branch gets to big, check out
# https://juliadocs.github.io/Documenter.jl/stable/man/hosting/#gh-pages-Branch
deploydocs(;repo="github.com/hexaeder/DynamicCascades.jl",
           push_preview=true)
