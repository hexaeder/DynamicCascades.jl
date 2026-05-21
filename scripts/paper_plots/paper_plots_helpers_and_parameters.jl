using MathTeXEngine
MathTeXEngine.set_texfont_family!(MathTeXEngine.FontFamily("TeXGyreHeros"))
set_theme!(fonts = (; regular = "TeX Gyre Heros Makie"))

fontsize = labelsize = 48
labellettersize = 4
linewidth = 4

# line_colors = [Makie.wong_colors()[3], Makie.wong_colors()[5], Makie.wong_colors()[6]]
line_colors = ["#5C3D99FF", "#56B4E9FF", "#D55E00FF"] # [deep indigo, Hellblau, Orange]